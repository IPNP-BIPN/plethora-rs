//! `clean_files.pl`: remove an intermediate once the next one is known good.
//!
//! ```text
//! code/clean_files.pl --rm-fastq $sample
//! ```
//!
//! Whole-genome data is large enough that the pipeline deletes as it goes. The
//! rule is a chain of counts: a file may go once the file derived from it holds
//! the number of reads it should. Getting that wrong deletes data, so the
//! decision is separated from the deletion here, and is what the tests exercise.
//!
//! The expected count comes either from the sequence index or, absent that,
//! from the FASTQ itself, in which case the check on the FASTQ is vacuous and
//! the chain still guards everything downstream.

use std::fmt;
use std::io::BufRead;
use std::path::{Path, PathBuf};

/// How the reads were sequenced, which decides the BED comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Pairing {
    /// Two records per fragment in the BAM, one line per fragment in the BED.
    Paired,
    /// One record per read throughout.
    Single,
}

/// What each stage holds, as far as it exists.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Counts {
    /// Reads across the sample's FASTQ files.
    pub fastq: Option<u64>,
    /// Records in `alignments/<sample>.bam`.
    pub bam: Option<u64>,
    /// Records in `results/<sample>_sorted.bam`.
    pub sorted_bam: Option<u64>,
    /// Lines in `results/<sample>.bed`.
    pub bed: Option<u64>,
}

impl Counts {
    /// How many stages are present, which is what decides there is work to do.
    #[must_use]
    pub fn present(&self) -> usize {
        [self.fastq, self.bam, self.sorted_bam, self.bed]
            .iter()
            .filter(|v| v.is_some())
            .count()
    }
}

/// A file the plan says can go.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Removal {
    /// The sample's FASTQ files. Only ever proposed when explicitly allowed,
    /// since they may be the only copy.
    Fastq,
    /// `alignments/<sample>.bam`.
    Bam,
    /// `results/<sample>_sorted.bam`.
    SortedBam,
    /// `results/<sample>.bed`, once the sorted BED exists.
    Bed,
}

/// A count that did not match, which stops the plan.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Mismatch {
    /// Which file disagreed.
    pub stage: &'static str,
    pub expected: u64,
    pub found: u64,
}

impl fmt::Display for Mismatch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "something is wrong with the {} for this sample: expected {} reads and counted {}",
            self.stage, self.expected, self.found
        )
    }
}

impl std::error::Error for Mismatch {}

/// What can be deleted, given what each stage holds.
///
/// `expected` is the count from the sequence index when one was given. Without
/// it the first stage present anchors the chain, exactly as upstream does:
/// FASTQ, else BAM, else sorted BAM. The comparison against whichever anchored
/// it is vacuous; the ones after it are not.
///
/// Nothing is proposed unless at least two stages are present: with one there
/// is nothing to corroborate it against.
///
/// # Errors
/// Returns the first count that disagreed. Upstream exits at that point, having
/// possibly already deleted something; here nothing is deleted, because the
/// caller has not been told to delete anything yet.
pub fn plan(
    counts: &Counts,
    expected: Option<u64>,
    pairing: Pairing,
    remove_fastq: bool,
) -> Result<Vec<Removal>, Mismatch> {
    Ok(plan_with_reasons(counts, expected, pairing, remove_fastq)?
        .into_iter()
        .map(|step| step.removal)
        .collect())
}

/// One removal, and the stage whose count justified it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Step {
    pub removal: Removal,
    /// `"bam"`, `"sorted bam"` or `"bed"`, spelled as upstream logs it.
    pub verified: &'static str,
}

/// [`plan`], keeping the stage that justified each removal.
///
/// The log line upstream writes names both, and `trim_qc_report.R` reads the
/// verified stage back out of it to decide how far a sample got. Dropping it
/// would leave the QC report with nothing to read.
///
/// # Errors
/// As [`plan`].
pub fn plan_with_reasons(
    counts: &Counts,
    expected: Option<u64>,
    pairing: Pairing,
    remove_fastq: bool,
) -> Result<Vec<Step>, Mismatch> {
    let mut removals: Vec<Step> = Vec::new();

    if counts.present() < 2 {
        return Ok(removals);
    }

    // The expectation falls back through the chain, as upstream does: the
    // manifest if there is one, else the FASTQ, else the BAM, else the sorted
    // BAM. Whichever anchors it, the comparison against it is vacuous and the
    // ones after it are not.
    let Some(expected) = expected
        .or(counts.fastq)
        .or(counts.bam)
        .or(counts.sorted_bam)
    else {
        return Ok(removals);
    };

    if let Some(found) = counts.fastq
        && found != expected
    {
        return Err(Mismatch {
            stage: "fastq file(s)",
            expected,
            found,
        });
    }

    if let Some(found) = counts.bam {
        if found != expected {
            return Err(Mismatch {
                stage: "bam file",
                expected,
                found,
            });
        }
        // Only when there is one: upstream loops over a list that is empty in
        // that case, so it proposes nothing either.
        if remove_fastq && counts.fastq.is_some() {
            removals.push(Step {
                removal: Removal::Fastq,
                verified: "bam",
            });
        }
    }

    if let Some(found) = counts.sorted_bam {
        if found != expected {
            return Err(Mismatch {
                stage: "sorted bam file",
                expected,
                found,
            });
        }
        if counts.bam.is_some() {
            removals.push(Step {
                removal: Removal::Bam,
                verified: "sorted bam",
            });
        }
    }

    if let Some(found) = counts.bed {
        // A paired BED holds one line per fragment, so half the records.
        let expected_bed = match pairing {
            Pairing::Paired => expected / 2,
            Pairing::Single => expected,
        };
        if found != expected_bed {
            return Err(Mismatch {
                stage: "bed file",
                expected: expected_bed,
                found,
            });
        }
        let already = |r: Removal, done: &[Step]| done.iter().any(|s| s.removal == r);
        if remove_fastq && counts.fastq.is_some() && !already(Removal::Fastq, &removals) {
            removals.push(Step {
                removal: Removal::Fastq,
                verified: "bed",
            });
        }
        if counts.bam.is_some() && !already(Removal::Bam, &removals) {
            removals.push(Step {
                removal: Removal::Bam,
                verified: "bed",
            });
        }
        if counts.sorted_bam.is_some() {
            removals.push(Step {
                removal: Removal::SortedBam,
                verified: "bed",
            });
        }
        // The unsorted BED goes once the sorted one exists; the caller checks
        // for it, since the count says nothing about whether it was written.
        removals.push(Step {
            removal: Removal::Bed,
            verified: "bed",
        });
    }

    Ok(removals)
}

/// The line `clean_files.pl` prints for one removal, which
/// `trim_qc_report.R` greps back out of `logs/clean_*.out`.
///
/// ```text
/// correct number of reads for:\tHG00250\tbed\tremoving:\tbam
/// ```
#[must_use]
pub fn log_line(sample: &str, step: Step) -> String {
    let removed = match step.removal {
        Removal::Fastq => "fastq",
        Removal::Bam => "bam",
        Removal::SortedBam => "sorted bam",
        Removal::Bed => "bed",
    };
    format!(
        "correct number of reads for:\t{sample}\t{}\tremoving:\t{removed}",
        step.verified
    )
}

/// Where a sample's files live, which upstream takes as `-f`, `-a` and `-b`.
#[derive(Debug, Clone, Copy)]
pub struct Paths<'a> {
    /// Holds one directory per sample.
    pub fastq: &'a Path,
    /// Holds `<sample>.bam`.
    pub alignments: &'a Path,
    /// Holds `<sample>.bed` and `<sample>_sorted.bam`.
    pub results: &'a Path,
}

/// What was found on disk for one sample.
#[derive(Debug, Clone, Default)]
pub struct Gathered {
    pub counts: Counts,
    /// The FASTQ files that were counted, in the order upstream's glob returns
    /// them, since those are the files a removal would delete.
    pub fastq_files: Vec<PathBuf>,
    /// Whether the file names look paired.
    ///
    /// Upstream sets this from whichever file the loop happened to end on,
    /// rather than from whether any mate 2 exists. A `glob` sorts, so the last
    /// file is the mate 2 when there is one and the answer comes out right;
    /// it is an accident of ordering, not a decision. This reports what those
    /// names say so a caller can notice when its configuration disagrees.
    pub named_pairing: Option<Pairing>,
}

/// Counts what each stage holds, the way `clean_files.pl` does.
///
/// A missing file is `None` rather than zero: the plan distinguishes "this
/// stage is not there" from "this stage is empty", and deleting on the second
/// would be wrong.
///
/// # Errors
/// Returns an error if a file that exists cannot be read.
pub fn gather(sample: &str, paths: &Paths) -> std::io::Result<Gathered> {
    let mut out = Gathered::default();

    // `glob "fastq/<sample>/*.fastq.gz"`, falling back to the uncompressed
    // form, which is what upstream does when the first glob is empty.
    let dir = paths.fastq.join(sample);
    for suffix in [".fastq.gz", ".fastq"] {
        out.fastq_files = list_with_suffix(&dir, suffix)?;
        if !out.fastq_files.is_empty() {
            break;
        }
    }
    if !out.fastq_files.is_empty() {
        let mut total = 0;
        for file in &out.fastq_files {
            total += count_fastq_records(file)?;
            // As upstream: taken from each file in turn, so the last one wins.
            out.named_pairing = Some(if is_mate_two(file) {
                Pairing::Paired
            } else {
                Pairing::Single
            });
        }
        out.counts.fastq = Some(total);
    }

    out.counts.bam = count_bam(&paths.alignments.join(format!("{sample}.bam")))?;
    out.counts.sorted_bam = count_bam(&paths.results.join(format!("{sample}_sorted.bam")))?;
    out.counts.bed = count_lines(&paths.results.join(format!("{sample}.bed")))?;

    Ok(out)
}

/// The files one removal covers, which is what a caller deletes.
#[must_use]
pub fn files_for(
    removal: Removal,
    sample: &str,
    paths: &Paths,
    gathered: &Gathered,
) -> Vec<PathBuf> {
    match removal {
        Removal::Fastq => gathered.fastq_files.clone(),
        Removal::Bam => vec![paths.alignments.join(format!("{sample}.bam"))],
        Removal::SortedBam => vec![paths.results.join(format!("{sample}_sorted.bam"))],
        Removal::Bed => vec![paths.results.join(format!("{sample}.bed"))],
    }
}

/// `*_2.fastq.gz`, `*_2.filt.fastq.gz` and `*_2_filtered.fastq.gz`, which is
/// upstream's `/_2((.filt)|(_filtered))*.fastq.*/`.
fn is_mate_two(path: &Path) -> bool {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or_default();
    let Some(rest) = name.split("_2").nth(1) else {
        return false;
    };
    let rest = rest.strip_prefix(".filt").unwrap_or(rest);
    let rest = rest.strip_prefix("_filtered").unwrap_or(rest);
    rest.starts_with(".fastq")
}

/// Names ending in `suffix`, sorted, which is what a shell glob yields.
fn list_with_suffix(dir: &Path, suffix: &str) -> std::io::Result<Vec<PathBuf>> {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return Ok(Vec::new());
    };
    let mut out: Vec<PathBuf> = Vec::new();
    for entry in entries {
        let path = entry?.path();
        if path
            .file_name()
            .and_then(|n| n.to_str())
            .is_some_and(|n| n.ends_with(suffix))
        {
            out.push(path);
        }
    }
    out.sort();
    Ok(out)
}

/// `gunzip -c file | paste - - - - | wc -l`, which counts groups of four and
/// so counts a trailing partial record as one.
fn count_fastq_records(path: &Path) -> std::io::Result<u64> {
    let lines = count_lines(path)?.unwrap_or(0);
    Ok(lines.div_ceil(4))
}

fn count_lines(path: &Path) -> std::io::Result<Option<u64>> {
    if !path.exists() {
        return Ok(None);
    }
    let mut n = 0;
    let mut reader = crate::io::open(path)?;
    let mut buffer = Vec::new();
    loop {
        buffer.clear();
        if reader.read_until(b'\n', &mut buffer)? == 0 {
            break;
        }
        n += 1;
    }
    Ok(Some(n))
}

/// `samtools view -c`.
fn count_bam(path: &Path) -> std::io::Result<Option<u64>> {
    if !path.exists() {
        return Ok(None);
    }
    Ok(Some(crate::bam::reader::read_bam(path)?.len() as u64))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts(
        fastq: Option<u64>,
        bam: Option<u64>,
        sorted: Option<u64>,
        bed: Option<u64>,
    ) -> Counts {
        Counts {
            fastq,
            bam,
            sorted_bam: sorted,
            bed,
        }
    }

    /// The counts come off disk in the same units upstream's shell pipelines
    /// produce: records for a FASTQ, lines for a BED.
    #[test]
    fn gather_counts_what_is_there() {
        let dir = tempfile::tempdir().expect("tempdir");
        let root = dir.path();
        std::fs::create_dir_all(root.join("fastq/S1")).expect("fastq dir");
        std::fs::create_dir_all(root.join("results")).expect("results dir");
        std::fs::create_dir_all(root.join("alignments")).expect("alignments dir");

        // Three records per mate: `paste - - - - | wc -l` counts groups of four.
        let record = "@r\nACGT\n+\nIIII\n";
        std::fs::write(root.join("fastq/S1/S1_1.fastq"), record.repeat(3)).expect("mate 1");
        std::fs::write(root.join("fastq/S1/S1_2.fastq"), record.repeat(3)).expect("mate 2");
        std::fs::write(root.join("results/S1.bed"), "a\nb\nc\n").expect("bed");

        let fastq = root.join("fastq");
        let alignments = root.join("alignments");
        let results = root.join("results");
        let paths = Paths {
            fastq: &fastq,
            alignments: &alignments,
            results: &results,
        };
        let got = gather("S1", &paths).expect("gather");

        assert_eq!(
            got.counts.fastq,
            Some(6),
            "three records in each of two files"
        );
        assert_eq!(got.counts.bed, Some(3));
        assert_eq!(got.counts.bam, None, "absent, not zero");
        assert_eq!(got.counts.sorted_bam, None);
        assert_eq!(got.fastq_files.len(), 2);
        assert_eq!(
            got.named_pairing,
            Some(Pairing::Paired),
            "the last file is a mate 2"
        );
    }

    /// The gz form is read the same way, since upstream pipes it through gunzip.
    #[test]
    fn gather_reads_compressed_fastq() {
        use std::io::Write as _;
        let dir = tempfile::tempdir().expect("tempdir");
        let root = dir.path();
        std::fs::create_dir_all(root.join("fastq/S1")).expect("fastq dir");
        let file = std::fs::File::create(root.join("fastq/S1/S1_1.fastq.gz")).expect("create");
        let mut w = flate2::write::GzEncoder::new(file, flate2::Compression::fast());
        w.write_all("@r\nACGT\n+\nIIII\n".repeat(5).as_bytes())
            .expect("write");
        w.finish().expect("finish");
        // An uncompressed file that the glob must not reach, since the first
        // pattern matched.
        std::fs::write(root.join("fastq/S1/other.fastq"), "@r\nA\n+\nI\n").expect("other");

        let fastq = root.join("fastq");
        let empty = root.join("nowhere");
        let paths = Paths {
            fastq: &fastq,
            alignments: &empty,
            results: &empty,
        };
        let got = gather("S1", &paths).expect("gather");
        assert_eq!(got.counts.fastq, Some(5));
        assert_eq!(got.fastq_files.len(), 1, "only the .gz glob was used");
        assert_eq!(got.named_pairing, Some(Pairing::Single));
    }

    /// The line `trim_qc_report.R` greps back out of the clean-up logs.
    #[test]
    fn the_log_line_is_upstreams() {
        let step = Step {
            removal: Removal::Bam,
            verified: "bed",
        };
        assert_eq!(
            log_line("HG00250", step),
            "correct number of reads for:\tHG00250\tbed\tremoving:\tbam"
        );
    }

    /// Every removal carries the stage that justified it, which is what makes
    /// the log line reconstructible.
    #[test]
    fn each_removal_names_what_verified_it() {
        let c = counts(Some(100), Some(100), None, Some(50));
        let steps = plan_with_reasons(&c, None, Pairing::Paired, true).expect("plan");
        assert!(!steps.is_empty());
        for step in &steps {
            assert!(
                ["bam", "sorted bam", "bed"].contains(&step.verified),
                "unexpected justification {}",
                step.verified
            );
        }
        // And the wrapper still answers the older question.
        let plain = plan(&c, None, Pairing::Paired, true).expect("plan");
        assert_eq!(plain, steps.iter().map(|s| s.removal).collect::<Vec<_>>());
    }

    /// One stage on its own corroborates nothing.
    #[test]
    fn a_single_stage_proposes_nothing() {
        let c = counts(Some(100), None, None, None);
        assert!(plan(&c, None, Pairing::Paired, true).unwrap().is_empty());
    }

    /// The BAM matching the FASTQ is what releases the FASTQ, and only when
    /// asked.
    #[test]
    fn the_fastq_goes_only_when_allowed() {
        let c = counts(Some(100), Some(100), None, None);
        assert_eq!(
            plan(&c, None, Pairing::Paired, true).unwrap(),
            [Removal::Fastq]
        );
        assert!(plan(&c, None, Pairing::Paired, false).unwrap().is_empty());
    }

    #[test]
    fn the_sorted_bam_releases_the_bam() {
        let c = counts(None, Some(100), Some(100), None);
        assert_eq!(
            plan(&c, None, Pairing::Paired, false).unwrap(),
            [Removal::Bam]
        );
    }

    /// A paired BED holds half as many lines as the BAM holds records.
    #[test]
    fn a_paired_bed_is_compared_against_half_the_reads() {
        let c = counts(Some(100), Some(100), None, Some(50));
        let plan = plan(&c, None, Pairing::Paired, false).unwrap();
        assert!(plan.contains(&Removal::Bam));
        assert!(plan.contains(&Removal::Bed));

        // The same count under single-end reads is wrong.
        let err = super::plan(&c, None, Pairing::Single, false).unwrap_err();
        assert_eq!(err.stage, "bed file");
        assert_eq!(err.expected, 100);
        assert_eq!(err.found, 50);
    }

    #[test]
    fn a_single_end_bed_is_compared_against_all_the_reads() {
        let c = counts(Some(100), Some(100), None, Some(100));
        assert!(
            plan(&c, None, Pairing::Single, false)
                .unwrap()
                .contains(&Removal::Bed)
        );
    }

    /// The point of the whole thing: a short file stops the plan rather than
    /// releasing what produced it.
    #[test]
    fn a_short_bam_stops_everything() {
        let c = counts(Some(100), Some(99), Some(100), Some(50));
        let err = plan(&c, None, Pairing::Paired, true).unwrap_err();
        assert_eq!(err.stage, "bam file");
        assert_eq!(err.expected, 100);
        assert_eq!(err.found, 99);
        assert!(
            err.to_string()
                .contains("expected 100 reads and counted 99")
        );
    }

    /// With a manifest the FASTQ is checked too, and a truncated download is
    /// caught before anything is deleted.
    #[test]
    fn a_manifest_makes_the_fastq_check_meaningful() {
        let c = counts(Some(90), Some(90), None, None);
        // Without the manifest the FASTQ defines the truth, so this passes.
        assert_eq!(
            plan(&c, None, Pairing::Paired, true).unwrap(),
            [Removal::Fastq]
        );
        // With it, the shortfall is caught.
        let err = plan(&c, Some(100), Pairing::Paired, true).unwrap_err();
        assert_eq!(err.stage, "fastq file(s)");
        assert_eq!(err.found, 90);
    }

    /// The full chain, everything present and consistent.
    #[test]
    fn a_complete_chain_releases_everything_it_can() {
        let c = counts(Some(100), Some(100), Some(100), Some(50));
        let plan = plan(&c, Some(100), Pairing::Paired, true).unwrap();
        assert!(plan.contains(&Removal::Fastq));
        assert!(plan.contains(&Removal::Bam));
        assert!(plan.contains(&Removal::SortedBam));
        assert!(plan.contains(&Removal::Bed));
        // Nothing is proposed twice, which would double-delete.
        let mut seen = plan.clone();
        seen.sort_by_key(|r| format!("{r:?}"));
        seen.dedup();
        assert_eq!(seen.len(), plan.len(), "a file was proposed twice");
    }

    /// With no FASTQ and no manifest the BAM anchors the chain, so the sorted
    /// BAM is still checked against something real.
    #[test]
    fn the_chain_anchors_on_the_first_stage_present() {
        let c = counts(None, Some(100), Some(100), None);
        assert_eq!(
            plan(&c, None, Pairing::Paired, true).unwrap(),
            [Removal::Bam]
        );

        // And a sorted BAM that disagrees with it is caught.
        let c = counts(None, Some(100), Some(99), None);
        assert_eq!(
            plan(&c, None, Pairing::Paired, true).unwrap_err().stage,
            "sorted bam file"
        );
    }

    #[test]
    fn nothing_at_all_proposes_nothing() {
        assert!(
            plan(&Counts::default(), None, Pairing::Paired, true)
                .unwrap()
                .is_empty()
        );
    }
}
