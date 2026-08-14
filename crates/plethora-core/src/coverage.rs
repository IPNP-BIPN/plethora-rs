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

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use plethora_compat::awk;

use crate::bam::bamtobed::{self, Aln, BedpeIter, is_pairable};
use crate::bam::namesort;
use crate::bam::reader::read_bam;
use crate::bed::{intersect, merge, sort};

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
) -> Result<Outputs, Box<dyn std::error::Error>> {
    let prefix = output.to_string_lossy().into_owned();
    let path = |suffix: &str| PathBuf::from(format!("{prefix}{suffix}"));

    let records: Vec<Aln> = read_bam(bam)?;
    let edited = path("_edited.bed");
    let mut bedpe_path = None;

    match pairing {
        Pairing::Paired => {
            let pairable: Vec<Aln> = records.into_iter().filter(|a| is_pairable(a.flags)).collect();
            let sorted_records = namesort::sort_by_name(pairable, namesort::DEFAULT_RUN_RECORDS, tmp_dir)?;

            let bedpe = path(".bed");
            let mut writer = BufWriter::new(File::create(&bedpe)?);
            let mut pairs = BedpeIter::new(sorted_records.into_iter());
            for record in pairs.by_ref() {
                writeln!(writer, "{record}")?;
            }
            writer.flush()?;
            if pairs.orphans() > 0 {
                eprintln!(
                    "warning: {} record(s) had no adjacent mate and were skipped",
                    pairs.orphans()
                );
            }
            bedpe_path = Some(bedpe.clone());

            crate::merge_pairs::run(&bedpe)?;
        }
        Pairing::Single => {
            // No pairing, no fragment reconstruction: each mapped alignment is
            // its own interval and goes straight to the sort.
            let mut writer = BufWriter::new(File::create(&edited)?);
            for record in &records {
                if let Some(line) = bamtobed::bed(record) {
                    writeln!(writer, "{line}")?;
                }
            }
            writer.flush()?;
        }
    }

    let sorted = path("_sorted.bed");
    sort::sort_lines(
        lines_of(&edited)?,
        sort::DEFAULT_RUN_LINES,
        tmp_dir,
        BufWriter::new(File::create(&sorted)?),
    )?;

    let temp = path("_temp.bed");
    intersect::intersect_wao(
        lines_of(reference)?,
        lines_of(&sorted)?,
        BufWriter::new(File::create(&temp)?),
    )?;

    // The awk that permutes thirteen columns into five, putting the domain name
    // where bedtools reads a chromosome so the merge groups by domain.
    let permuted = permute(lines_of(&temp)?);
    let coverage = path("_coverage.bed");
    merge::merge_sum(permuted, 5, BufWriter::new(File::create(&coverage)?))?;

    let read_depth = path("_read_depth.bed");
    write_read_depth(lines_of(&coverage)?, BufWriter::new(File::create(&read_depth)?))?;

    // Upstream removes these two and keeps the rest.
    let _ = std::fs::remove_file(&edited);
    let _ = std::fs::remove_file(&temp);

    Ok(Outputs {
        bedpe: bedpe_path,
        sorted,
        coverage,
        read_depth,
    })
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

/// Reads a file as lines, skipping ones that fail to decode.
fn lines_of(path: &Path) -> io::Result<impl Iterator<Item = String>> {
    Ok(BufReader::new(File::open(path)?).lines().map_while(Result::ok))
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
