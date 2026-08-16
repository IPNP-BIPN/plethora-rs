//! `make_bed.sh`: from an alignment to a coverage figure per domain.
//!
//! ```text
//! samtools sort -n -@ 12 -m 2G -o ${output}_sorted.bam $bam
//! bedtools bamtobed -split -bedpe -i ${output}_sorted.bam > $output.bed
//! code/merge_pairs.pl $output.bed
//! sort -k 1,1 -k 2,2n -T ./ ${output}_edited.bed > ${output}_sorted.bed
//! bedtools intersect -wao -sorted -a $reference_bed -b ${output}_sorted.bed > ${output}_temp.bed
//! awk 'OFS="\t" {print $4,$2,$3,$1,$13}' ${output}_temp.bed | bedtools merge -c 5 -o sum -i - > ${output}_coverage.bed
//! awk 'OFS="\t" { print $1, $4 / ($3 - $2 + 1)}' ${output}_coverage.bed > ${output}_read_depth.bed
//! ```
//!
//! Every intermediate file keeps its upstream name, so a run can be compared
//! against the original pipeline step by step rather than only at the end.
//!
//! The last line carries a quirk worth naming: the divisor is
//! `end - start + 1`, one more than the domain's length. A 1000 bp domain is
//! divided by 1001. It is kept, because the published copy numbers were
//! produced with it and correcting it would silently rescale every result by
//! 0.1%. See `DIVERGENCES.md`.

use std::io::{self, BufRead, Write};
use std::path::{Path, PathBuf};

use plethora_compat::awk;

use crate::bam::bamtobed::{self, BedpeIter, is_pairable};
use crate::bam::namesort;
use crate::bam::reader::read_bam_streaming;
use crate::bed::{intersect, merge, sort};
use crate::io::Compress;

/// How the reads were sequenced, which decides how they become intervals.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Pairing {
    /// Mates are combined into fragments by `merge_pairs`.
    Paired,
    /// Each alignment becomes one interval; `merge_pairs` is not involved.
    Single,
}

impl std::str::FromStr for Pairing {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "paired" => Ok(Self::Paired),
            "single" => Ok(Self::Single),
            other => Err(format!("unknown pairing scheme: {other}")),
        }
    }
}

/// The files one run produces, by their upstream names.
#[derive(Debug, Clone)]
pub struct Outputs {
    /// `${output}.bed`, the BEDPE. Paired mode only.
    pub bedpe: Option<PathBuf>,
    /// `${output}_sorted.bed`, the sorted intervals.
    pub sorted: PathBuf,
    /// `${output}_coverage.bed`, overlapping bases summed per domain.
    pub coverage: PathBuf,
    /// `${output}_read_depth.bed`, the figure the pipeline is after.
    pub read_depth: PathBuf,
}

/// Runs the whole chain for one sample.
///
/// `output` is a prefix, not a file: upstream passes `results/HG00250` and gets
/// `results/HG00250_read_depth.bed` among others.
///
/// Secondary and supplementary alignments are dropped before pairing. bedtools
/// does not do that, and cannot: its reader consumes records two at a time, so
/// a third record under one name desynchronises everything after it. Under
/// bowtie2 the situation never arises; under BWA-MEM it always does. See
/// [`crate::bam::bamtobed::is_pairable`].
///
/// # Errors
/// Returns an error if any stage fails to read or write.
pub fn make_bed(
    bam: &Path,
    reference: &Path,
    pairing: Pairing,
    output: &Path,
    tmp_dir: &Path,
    compress: Compress,
) -> Result<Outputs, Box<dyn std::error::Error>> {
    let prefix = output.to_string_lossy().into_owned();
    // Every output takes the same form, so a run is all compressed or all not
    // and the stage that reads a file always looks for the name the stage
    // before it wrote.
    let path = |suffix: &str| compress.apply(&PathBuf::from(format!("{prefix}{suffix}")));

    let edited = path("_edited.bed");
    let mut bedpe_path = None;

    // Nothing here holds the alignment in memory. A 30x whole genome is over a
    // billion records, and at 181 bytes each a `Vec` of them would be 137 GB;
    // every stage reads them once, in order, so none of them needed one.
    match pairing {
        Pairing::Paired => {
            let mut failed: Option<std::io::Error> = None;
            let pairable = read_bam_streaming(bam)?
                .map_while(|r| capture(r, &mut failed))
                .filter(|a| is_pairable(a.flags));
            let sorted_records =
                namesort::sort_by_name_streaming(pairable, namesort::DEFAULT_RUN_RECORDS, tmp_dir)?;
            check(failed)?;

            let mut failed: Option<std::io::Error> = None;
            let bedpe = path(".bed");
            let mut writer = crate::io::create(&bedpe)?;
            let mut pairs = BedpeIter::new(sorted_records.map_while(|r| capture(r, &mut failed)));
            for record in pairs.by_ref() {
                writeln!(writer, "{record}")?;
            }
            let orphans = pairs.orphans();
            drop(pairs);
            // Dropped, not merely flushed, and before anything reads the file
            // back. A BGZF writer finalises on drop: that is when its workers
            // are joined and the end-of-file block is written, so a flush alone
            // leaves the last blocks unwritten and `run_to` below would open a
            // truncated file.
            writer.flush()?;
            drop(writer);
            check(failed)?;
            if orphans > 0 {
                eprintln!("warning: {orphans} record(s) had no adjacent mate and were skipped");
            }
            bedpe_path = Some(bedpe.clone());
        }
        Pairing::Single => {
            // No pairing, no fragment reconstruction: each mapped alignment is
            // its own interval and goes straight to the sort.
            let mut writer = crate::io::create(&edited)?;
            for record in read_bam_streaming(bam)? {
                if let Some(line) = bamtobed::bed(&record?) {
                    writeln!(writer, "{line}")?;
                }
            }
            writer.flush()?;
        }
    }

    // `_edited.bed` and `_temp.bed` never exist. Upstream writes them and
    // deletes them on the way out; here each is a pipe into the stage that
    // reads it, so nothing observable changes and the run's peak disk falls by
    // more than half. See `crate::io::piped_lines`.
    let mut failed: Option<std::io::Error> = None;
    let sorted = path("_sorted.bed");
    {
        // A BEDPE exists exactly when the paired arm ran, so it selects the
        // source rather than `pairing` being consulted a second time.
        let edited: Box<dyn Iterator<Item = std::io::Result<String>>> = match bedpe_path.clone() {
            Some(bedpe) => Box::new(crate::io::piped_lines(move |out| {
                crate::merge_pairs::emit_to(&bedpe, out).map(|_stats| ())
            })?),
            // The single-end arm has no fragments to reconstruct, so its
            // intervals were written straight out above.
            None => Box::new(lines_of(&edited)?.map(Ok)),
        };
        sort::sort_lines(
            edited.map_while(|l| capture(l, &mut failed)),
            sort::DEFAULT_RUN_LINES,
            tmp_dir,
            crate::io::create(&sorted)?,
        )?;
    }
    check(failed)?;

    let reference_lines: Vec<String> = lines_of(reference)?.collect();
    let sorted_for_intersect = sorted.clone();
    let temp = crate::io::piped_lines(move |out| {
        intersect::intersect_wao(
            reference_lines.into_iter(),
            lines_of(&sorted_for_intersect)?,
            out,
        )
        .map_err(|e| std::io::Error::other(e.to_string()))
    })?;

    // The awk that permutes thirteen columns into five, putting the domain name
    // where bedtools reads a chromosome so the merge groups by domain.
    let mut failed: Option<std::io::Error> = None;
    let coverage = path("_coverage.bed");
    merge::merge_sum(
        permute(temp.map_while(|l| capture(l, &mut failed))),
        5,
        crate::io::create(&coverage)?,
    )?;
    check(failed)?;

    let read_depth = path("_read_depth.bed");
    write_read_depth(lines_of(&coverage)?, crate::io::create(&read_depth)?)?;

    // The single-end arm is the only one that still writes this.
    let _ = std::fs::remove_file(&edited);

    Ok(Outputs {
        bedpe: bedpe_path,
        sorted,
        coverage,
        read_depth,
    })
}

/// Yields the value and stops the iterator on the first error, keeping it.
///
/// `map_while(Result::ok)` would do the same thing and lose the error, which
/// is how a gzipped BEDPE once read as zero proper pairs instead of failing.
fn capture<T>(result: std::io::Result<T>, failed: &mut Option<std::io::Error>) -> Option<T> {
    match result {
        Ok(value) => Some(value),
        Err(e) => {
            *failed = Some(e);
            None
        }
    }
}

/// Turns a captured error back into a failure.
fn check(failed: Option<std::io::Error>) -> Result<(), std::io::Error> {
    failed.map_or(Ok(()), Err)
}

/// `awk 'OFS="\t" {print $4,$2,$3,$1,$13}'`.
///
/// Lines with fewer than thirteen fields are dropped rather than emitted with
/// empty columns; awk would print empty strings, which the merge would then
/// refuse to parse.
fn permute<I: Iterator<Item = String>>(lines: I) -> impl Iterator<Item = String> {
    lines.filter_map(|line| {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 13 {
            return None;
        }
        Some(format!("{}\t{}\t{}\t{}\t{}", f[3], f[1], f[2], f[0], f[12]))
    })
}

/// `awk 'OFS="\t" { print $1, $4 / ($3 - $2 + 1)}'`.
///
/// The divisor is one more than the domain length; see the module note. Numbers
/// are spelled the way awk spells them, which is what gives the published
/// figures their six significant digits.
///
/// # Errors
/// Returns an error if writing fails.
pub fn write_read_depth<I, W>(lines: I, mut out: W) -> io::Result<()>
where
    I: Iterator<Item = String>,
    W: Write,
{
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 4 {
            continue;
        }
        let (Ok(start), Ok(end), Ok(sum)) = (
            f[1].parse::<i64>(),
            f[2].parse::<i64>(),
            f[3].parse::<f64>(),
        ) else {
            continue;
        };
        let depth = sum / (end - start + 1) as f64;
        writeln!(out, "{}\t{}", f[0], awk::print_number(depth))?;
    }
    out.flush()
}

/// Reads a file as lines, decompressing it if it is gzipped and skipping lines
/// that fail to decode.
fn lines_of(path: &Path) -> io::Result<impl Iterator<Item = String>> {
    Ok(crate::io::open(path)?.lines().map_while(Result::ok))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_permutation_moves_the_domain_name_first() {
        let line = "chr1\t100\t200\tdomA\t255\t+\tchr1\t150\t250\tr1\t60\t+\t50";
        let out: Vec<String> = permute(std::iter::once(line.to_string())).collect();
        assert_eq!(out, ["domA\t100\t200\tchr1\t50"]);
    }

    #[test]
    fn short_lines_are_dropped_by_the_permutation() {
        let out: Vec<String> = permute(std::iter::once("chr1\t100\t200".to_string())).collect();
        assert!(out.is_empty());
    }

    /// The divisor is the domain length plus one, which is upstream's.
    #[test]
    fn depth_divides_by_the_length_plus_one() {
        let mut out = Vec::new();
        write_read_depth(std::iter::once("domA\t100\t200\t101".to_string()), &mut out).unwrap();
        // 101 / (200 - 100 + 1) = 1 exactly, and awk prints an integer.
        assert_eq!(String::from_utf8(out).unwrap(), "domA\t1\n");
    }

    #[test]
    fn an_uncovered_domain_reports_zero() {
        let mut out = Vec::new();
        write_read_depth(std::iter::once("domC\t500\t600\t0".to_string()), &mut out).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "domC\t0\n");
    }

    /// Six significant digits, as awk's OFMT gives.
    #[test]
    fn depth_is_printed_the_way_awk_prints_it() {
        let mut out = Vec::new();
        write_read_depth(std::iter::once("domA\t100\t200\t60".to_string()), &mut out).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "domA\t0.594059\n");
    }

    #[test]
    fn pairing_parses_the_two_upstream_spellings() {
        assert_eq!("paired".parse::<Pairing>().unwrap(), Pairing::Paired);
        assert_eq!("single".parse::<Pairing>().unwrap(), Pairing::Single);
        assert!("both".parse::<Pairing>().is_err());
    }
}
