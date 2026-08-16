//! `merge_pairs.pl`: fragments from proper pairs, extended reads from the rest.
//!
//! ```text
//! code/merge_pairs.pl $output.bed
//! ```
//!
//! Reads the BEDPE that `bamtobed` produced and writes `*_edited.bed`. A pair
//! it considers proper collapses into one interval spanning both reads. Anything
//! else is emitted as up to two single-end reads, each extended outwards by half
//! a fragment length drawn from the sample's own fragment-size distribution.
//!
//! The draw is seeded from an MD5 of the input line, so it is deterministic per
//! read rather than per run. That is what makes reproducing this script a
//! question of reproducing Perl's RANDLIB exactly; see
//! [`plethora_compat::randlib`].
//!
//! Two passes over the same file: the first measures the fragment-size
//! distribution, the second writes the output using it.

use std::io::{self, BufRead, Write};
use std::path::Path;

use md5::{Digest, Md5};
use plethora_compat::randlib::Randlib;

/// The inner-distance ceiling before the sample's own distribution is known.
///
/// `merge_pairs.pl` declares this once, at file scope, and overwrites it after
/// the first pass. So the two passes do not apply the same test: the first
/// classifies pairs against 800, the second against `mean + 5 * sd`. That is
/// not a bug to fix, it is what decides which reads contribute to the
/// distribution in the first place.
pub const INITIAL_MAX_INNER_DISTANCE: i64 = 800;

/// Where the distance sampling stops.
///
/// `last if ($#distance > $sufficient_number_of_reads)` compares the last
/// index, not the length, so collection stops once the vector holds
/// 50,000,002 entries. A whole-genome sample reaches this, which is why the
/// record order coming out of the name sort reaches the output.
pub const SUFFICIENT_NUMBER_OF_READS: usize = 50_000_000;

/// The fragment-size distribution measured in the first pass.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FragmentStats {
    /// Mean inner distance, rounded to an integer as `sprintf("%.0f", ...)`.
    pub mean: i64,
    /// Standard deviation, likewise rounded.
    ///
    /// A population standard deviation, dividing by n rather than n - 1, and
    /// computed against the already-rounded mean.
    pub sd: i64,
    /// `mean + 5 * sd`, the ceiling the second pass uses.
    pub max_inner_distance: i64,
    /// How many distances went into the estimate.
    pub n: usize,
}

/// One parsed BEDPE line.
///
/// Field names follow the BEDPE column order rather than read 1 and read 2,
/// because `bamtobed` orders the blocks by position: the first block is not
/// necessarily the first mate. See [`crate::bam::bamtobed::bedpe`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BedpeLine<'a> {
    pub chrom1: &'a str,
    pub start1: i64,
    pub end1: i64,
    pub chrom2: &'a str,
    pub start2: i64,
    pub end2: i64,
    pub name: &'a str,
    pub score: &'a str,
    pub strand1: &'a str,
    pub strand2: &'a str,
}

impl<'a> BedpeLine<'a> {
    /// Splits a line into its ten fields.
    ///
    /// Perl's `split(/\t/, $line)` on a short line simply yields fewer
    /// elements, and the comparisons that follow treat a missing field as
    /// undefined. Here a short line is rejected instead, since every line
    /// `bamtobed` writes has ten fields and a shorter one means the input is
    /// not what the stage expects.
    #[must_use]
    pub fn parse(line: &'a str) -> Option<Self> {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            return None;
        }
        // "." positions come through as -1; Perl numifies them to -1 too.
        let num = |s: &str| s.parse::<i64>().unwrap_or(0);
        Some(Self {
            chrom1: f[0],
            start1: num(f[1]),
            end1: num(f[2]),
            chrom2: f[3],
            start2: num(f[4]),
            end2: num(f[5]),
            name: f[6],
            score: f[7],
            strand1: f[8],
            strand2: f[9],
        })
    }

    /// `isAProperPair`, against the ceiling in force at the time.
    ///
    /// The ceiling is passed in rather than read from a global, so the two
    /// passes cannot accidentally share one. They must not: the first uses 800
    /// and the second uses `mean + 5 * sd`.
    #[must_use]
    pub fn is_a_proper_pair(&self, max_inner_distance: i64) -> bool {
        // An unmapped end is never part of a proper pair.
        if self.chrom1 == "." || self.chrom2 == "." {
            return false;
        }
        if self.chrom1 != self.chrom2 {
            return false;
        }
        // Same strand means the mates do not face each other.
        if self.strand1 == self.strand2 {
            return false;
        }
        if self.inner_distance() > max_inner_distance {
            return false;
        }
        // Overlapping reads count as a proper pair.
        if (self.start1 >= self.start2 && self.start1 <= self.end2)
            || (self.end1 >= self.start2 && self.end1 <= self.end2)
        {
            return true;
        }
        // Block 1 starting after block 2 would mean the fragment runs backwards.
        // bamtobed orders the blocks, so this is close to unreachable, and it is
        // kept because it is what the Perl tests.
        if self.start1 > self.start2 {
            return false;
        }
        true
    }

    /// The gap between the two reads: `start2 - end1`.
    #[must_use]
    pub const fn inner_distance(&self) -> i64 {
        self.start2 - self.end1
    }
}

/// No pair looked proper, so there is no distribution to draw from.
///
/// The Perl divides by zero here and dies with "Illegal division by zero".
/// Saying what went wrong is more use: on a real sample it means the alignment
/// is against the wrong reference, or that the BAM was not name-sorted so
/// `bamtobed -bedpe` paired nothing.
#[derive(Debug, thiserror::Error)]
#[error(
    "no proper pairs found, so there is no fragment size distribution. The \
     alignment may be against a different reference than the domains, or the \
     BAM may not have been name-sorted"
)]
pub struct NoProperPairs;

/// First pass: measure the fragment-size distribution.
///
/// Distances are collected only from pairs that look proper against
/// [`INITIAL_MAX_INNER_DISTANCE`], and collection stops after
/// [`SUFFICIENT_NUMBER_OF_READS`].
///
/// # Errors
/// Returns [`NoProperPairs`] when nothing qualified.
pub fn measure<I: Iterator<Item = String>>(lines: I) -> Result<FragmentStats, NoProperPairs> {
    let mut distance: Vec<i64> = Vec::new();

    for line in lines {
        if let Some(p) = BedpeLine::parse(&line)
            && p.is_a_proper_pair(INITIAL_MAX_INNER_DISTANCE)
        {
            distance.push(p.inner_distance());
        }
        // `$#distance > 50_000_000` compares the last index, so this triggers
        // at 50,000,002 entries, not 50,000,001.
        if distance.len() > SUFFICIENT_NUMBER_OF_READS + 1 {
            break;
        }
    }

    if distance.is_empty() {
        return Err(NoProperPairs);
    }

    let n = distance.len();
    let mut mean = 0.0_f64;
    for &d in &distance {
        mean += d as f64;
    }
    mean /= n as f64;
    // sprintf("%.0f") rounds half to even, which is what round_ties_even does.
    // Formatting through "{:.0}" and parsing back agrees on all four million
    // values we checked, halves included; this way is simply infallible.
    let mean = mean.round_ties_even() as i64;

    // Population variance, against the rounded mean.
    let mut sd = 0.0_f64;
    for &d in &distance {
        let centred = (d - mean) as f64;
        sd += centred * centred;
    }
    sd /= n as f64;
    let sd = sd.sqrt().round_ties_even() as i64;

    Ok(FragmentStats {
        mean,
        sd,
        max_inner_distance: mean + 5 * sd,
        n,
    })
}

/// How far to extend one unpaired read, drawn from the sample's distribution.
///
/// Seeded from an MD5 of the whole input line, so the same read always gets the
/// same extension however many times the pipeline runs. Half a fragment,
/// rounded, and never negative: reads that would overlap are simply not
/// trimmed.
///
/// # Panics
/// Panics if the rounded draw does not parse back as an integer, which cannot
/// happen for a finite mean and deviation.
#[must_use]
pub fn inner_distance_for(line: &str, stats: &FragmentStats) -> i64 {
    let digest = Md5::digest(line.as_bytes());
    let phrase = format!("{digest:x}");

    let mut rng = Randlib::new();
    rng.set_seed_from_phrase(phrase.as_bytes());
    let draw = rng.gennor(stats.mean as f64, stats.sd as f64) / 2.0;

    let rounded: i64 = format!("{draw:.0}")
        .parse()
        .expect("a rounded draw is an integer");
    rounded.max(0)
}

/// Second pass: write one line per fragment, or up to two per broken pair.
///
/// # Errors
/// Returns an error if the output cannot be written.
///
/// # Panics
/// Panics if a rounded draw does not parse back as an integer, which cannot
/// happen for a finite mean and deviation.
pub fn emit<I, W>(lines: I, stats: &FragmentStats, mut out: W) -> io::Result<()>
where
    I: Iterator<Item = String>,
    W: Write,
{
    for line in lines {
        let Some(p) = BedpeLine::parse(&line) else {
            continue;
        };

        if p.is_a_proper_pair(stats.max_inner_distance) {
            // One interval spanning both reads. The strand kept is block 1's.
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}",
                p.chrom1,
                p.start1.min(p.start2),
                p.end1.max(p.end2),
                p.name,
                p.score,
                p.strand1
            )?;
            continue;
        }

        let extension = inner_distance_for(&line, stats);

        // Each end that has an alignment is emitted alone, extended away from
        // where its mate would have been.
        for (chrom, start, end, strand) in [
            (p.chrom1, p.start1, p.end1, p.strand1),
            (p.chrom2, p.start2, p.end2, p.strand2),
        ] {
            if chrom == "." {
                continue;
            }
            let mut start = start;
            let mut end = end;
            if strand == "-" {
                start -= extension;
            }
            if strand == "+" {
                end += extension;
            }
            // Extending past the start of the chromosome clamps to zero, and
            // nothing clamps the far end: a read can be extended past the end
            // of the reference. That is upstream's behaviour and the intersect
            // stage drops the overhang anyway.
            start = start.max(0);
            writeln!(
                out,
                "{chrom}\t{start}\t{end}\t{}\t{}\t{strand}",
                p.name, p.score
            )?;
        }
    }
    out.flush()
}

/// Runs both passes over a BEDPE, writing the intervals to `output`.
///
/// The input is read twice, as the Perl does, rather than held in memory: a
/// whole-genome BEDPE runs to tens of gigabytes. Both paths go through
/// [`crate::io`], so a gzipped input reads and a `.gz` output compresses.
///
/// The output path is given rather than derived. Upstream derives it with a
/// substitution that cannot see a `.gz` suffix, and a caller that compresses
/// its intermediates already knows what it called them.
///
/// # Errors
/// Returns an error if the input cannot be read, the output cannot be written,
/// or no pair looked proper.
pub fn run_to(path: &Path, output: &Path) -> io::Result<FragmentStats> {
    let mut out = crate::io::create(output)?;
    let stats = emit_to(path, &mut out)?;
    // Dropped before returning, so whatever reads the file next finds a
    // finished one: a BGZF writer finalises on drop, not on flush.
    out.flush()?;
    drop(out);
    Ok(stats)
}

/// Both passes over the BEDPE, writing the edited intervals wherever `out`
/// goes.
///
/// Separated from [`run_to`] because the intervals usually go straight into the
/// sort rather than to a file: upstream writes `*_edited.bed` and deletes it a
/// few lines later, so there is nothing to keep.
///
/// # Errors
/// Returns an error if the input cannot be read, if it holds no proper pair, or
/// if the output cannot be written.
pub fn emit_to<W: Write>(path: &Path, out: W) -> io::Result<FragmentStats> {
    let read_lines =
        || -> io::Result<_> { Ok(crate::io::open(path)?.lines().map_while(Result::ok)) };

    let stats =
        measure(read_lines()?).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    emit(read_lines()?, &stats, out)?;
    Ok(stats)
}

/// Runs both passes, writing `*_edited.bed` beside the input as upstream does.
///
/// `$outfile =~ s/.bed$/_edited.bed/` is a regex where the dot matches any
/// character, so it also rewrites a file ending in "Xbed"; only the literal case
/// arises here. A `.gz` suffix is carried across, which the Perl cannot do.
///
/// # Errors
/// Returns an error if the name does not end in `.bed` or `.bed.gz`, or if
/// either pass fails.
pub fn run(path: &Path) -> io::Result<FragmentStats> {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "the input has no file name"))?;

    let (stem, suffix) = if let Some(stem) = name.strip_suffix(".bed.gz") {
        (stem, ".bed.gz")
    } else if let Some(stem) = name.strip_suffix(".bed") {
        (stem, ".bed")
    } else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "expected a .bed or .bed.gz file",
        ));
    };

    let output = path.with_file_name(format!("{stem}_edited{suffix}"));
    run_to(path, &output)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One BEDPE line. Grouped rather than passed as eight arguments so the
    /// call sites read as two blocks and a pair of strands.
    struct L<'a>(&'a str, i64, i64, &'a str, i64, i64, &'a str, &'a str);

    #[allow(clippy::many_single_char_names)]
    fn line_of(L(c1, s1, e1, c2, s2, e2, strand1, strand2): L<'_>) -> String {
        format!("{c1}\t{s1}\t{e1}\t{c2}\t{s2}\t{e2}\tread\t60\t{strand1}\t{strand2}")
    }

    #[test]
    fn a_facing_pair_within_range_is_proper() {
        let l = line_of(L("chr1", 100, 150, "chr1", 300, 350, "+", "-"));
        let p = BedpeLine::parse(&l).unwrap();
        assert!(p.is_a_proper_pair(800));
        assert_eq!(p.inner_distance(), 150);
    }

    #[test]
    fn an_unmapped_end_is_never_proper() {
        let l = line_of(L(".", -1, -1, "chr1", 300, 350, ".", "-"));
        assert!(!BedpeLine::parse(&l).unwrap().is_a_proper_pair(800));
    }

    #[test]
    fn different_chromosomes_are_never_proper() {
        let l = line_of(L("chr1", 100, 150, "chr2", 300, 350, "+", "-"));
        assert!(!BedpeLine::parse(&l).unwrap().is_a_proper_pair(800));
    }

    #[test]
    fn same_strand_is_never_proper() {
        let l = line_of(L("chr1", 100, 150, "chr1", 300, 350, "+", "+"));
        assert!(!BedpeLine::parse(&l).unwrap().is_a_proper_pair(800));
    }

    /// The ceiling is the only thing separating a proper pair from a pair of
    /// extended singletons, and it differs between the two passes.
    #[test]
    fn the_ceiling_decides() {
        let l = line_of(L("chr1", 100, 150, "chr1", 1000, 1050, "+", "-"));
        let p = BedpeLine::parse(&l).unwrap();
        assert_eq!(p.inner_distance(), 850);
        assert!(!p.is_a_proper_pair(800), "beyond the first-pass ceiling");
        assert!(
            p.is_a_proper_pair(900),
            "within a wider second-pass ceiling"
        );
    }

    /// Overlapping reads short-circuit to proper, even though the inner
    /// distance is negative and the start ordering would otherwise reject them.
    #[test]
    fn overlapping_reads_are_proper() {
        let l = line_of(L("chr1", 300, 360, "chr1", 320, 380, "+", "-"));
        let p = BedpeLine::parse(&l).unwrap();
        assert!(p.inner_distance() < 0);
        assert!(p.is_a_proper_pair(800));
    }

    #[test]
    fn a_proper_pair_collapses_to_one_span() {
        let l = line_of(L("chr1", 100, 150, "chr1", 300, 350, "+", "-"));
        let stats = FragmentStats {
            mean: 200,
            sd: 30,
            max_inner_distance: 350,
            n: 10,
        };
        let mut out = Vec::new();
        emit(std::iter::once(l), &stats, &mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "chr1\t100\t350\tread\t60\t+\n"
        );
    }

    /// A broken pair becomes two reads, each pushed outwards on its own strand.
    #[test]
    fn a_broken_pair_becomes_two_extended_reads() {
        let l = line_of(L("chr1", 100, 150, "chr1", 90000, 90050, "+", "-"));
        let stats = FragmentStats {
            mean: 200,
            sd: 30,
            max_inner_distance: 350,
            n: 10,
        };
        let extension = inner_distance_for(&l, &stats);
        assert!(extension > 0);

        let mut out = Vec::new();
        emit(std::iter::once(l), &stats, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 2);
        // Forward read grows at its end, reverse read grows at its start.
        assert_eq!(
            lines[0],
            format!("chr1\t100\t{}\tread\t60\t+", 150 + extension)
        );
        assert_eq!(
            lines[1],
            format!("chr1\t{}\t90050\tread\t60\t-", 90000 - extension)
        );
    }

    /// An unmapped end contributes no line at all, so a half-mapped pair yields
    /// exactly one read.
    #[test]
    fn an_unmapped_end_emits_nothing() {
        let l = line_of(L(".", -1, -1, "chr1", 300, 350, ".", "+"));
        let stats = FragmentStats {
            mean: 200,
            sd: 30,
            max_inner_distance: 350,
            n: 10,
        };
        let mut out = Vec::new();
        emit(std::iter::once(l), &stats, &mut out).unwrap();
        assert_eq!(String::from_utf8(out).unwrap().lines().count(), 1);
    }

    /// Extending past the start of a chromosome clamps to zero.
    #[test]
    fn extension_clamps_at_the_chromosome_start() {
        let l = line_of(L("chr1", 5, 20, "chr1", 90000, 90050, "-", "+"));
        let stats = FragmentStats {
            mean: 400,
            sd: 50,
            max_inner_distance: 650,
            n: 10,
        };
        let mut out = Vec::new();
        emit(std::iter::once(l.clone()), &stats, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(
            text.lines().next().unwrap().starts_with("chr1\t0\t"),
            "expected a clamped start, got {text}"
        );
    }

    /// The same line must always draw the same extension, whatever else ran
    /// before it. That is the whole point of seeding from the line's MD5.
    #[test]
    fn the_extension_is_deterministic_per_line() {
        let stats = FragmentStats {
            mean: 317,
            sd: 45,
            max_inner_distance: 542,
            n: 10,
        };
        let a = line_of(L("chr1", 100, 150, "chr1", 90000, 90050, "+", "-"));
        let b = line_of(L("chr1", 200, 250, "chr1", 90000, 90050, "+", "-"));

        assert_eq!(
            inner_distance_for(&a, &stats),
            inner_distance_for(&a, &stats)
        );
        // Interleaving another draw must not shift the answer.
        let first = inner_distance_for(&a, &stats);
        let _ = inner_distance_for(&b, &stats);
        assert_eq!(inner_distance_for(&a, &stats), first);
    }

    #[test]
    fn the_extension_is_never_negative() {
        // A distribution wide enough that half a draw often comes out below zero.
        let stats = FragmentStats {
            mean: 10,
            sd: 500,
            max_inner_distance: 2510,
            n: 10,
        };
        for i in 0..200 {
            let l = line_of(L("chr1", i, i + 50, "chr1", 90000, 90050, "+", "-"));
            assert!(inner_distance_for(&l, &stats) >= 0);
        }
    }

    #[test]
    fn measure_uses_the_population_deviation_and_a_rounded_mean() {
        // Inner distances of 100, 200 and 300: mean 200, population sd 81.65.
        let lines = vec![
            line_of(L("chr1", 0, 50, "chr1", 150, 200, "+", "-")),
            line_of(L("chr1", 0, 50, "chr1", 250, 300, "+", "-")),
            line_of(L("chr1", 0, 50, "chr1", 350, 400, "+", "-")),
        ];
        let stats = measure(lines.into_iter()).expect("the corpus has proper pairs");
        assert_eq!(stats.n, 3);
        assert_eq!(stats.mean, 200);
        assert_eq!(stats.sd, 82, "population deviation, rounded");
        assert_eq!(stats.max_inner_distance, 200 + 5 * 82);
    }

    /// Pairs beyond 800 do not contribute to the distribution, however many
    /// there are.
    #[test]
    fn the_first_pass_ceiling_filters_the_sample() {
        let lines = vec![
            line_of(L("chr1", 0, 50, "chr1", 150, 200, "+", "-")),
            // Inner distance 2000: excluded from the estimate.
            line_of(L("chr1", 0, 50, "chr1", 2050, 2100, "+", "-")),
        ];
        let stats = measure(lines.into_iter()).expect("the corpus has proper pairs");
        assert_eq!(stats.n, 1);
        assert_eq!(stats.mean, 100);
    }

    #[test]
    fn short_lines_are_skipped() {
        assert!(BedpeLine::parse("chr1\t1\t2").is_none());
    }
}
