//! `trim_qc_report.R`: how far each sample got, and which ones look wrong.
//!
//! Reads the trimming log, the sequence index, the sample list and the
//! alignment report, and answers three questions: did every file arrive, how
//! much was lost to trimming, and how far down the pipeline each sample is.
//!
//! **One deliberate divergence.** The script ends with a bare call to
//! `cleanup_old_files()`, so merely producing the report deletes FASTQ files:
//!
//! ```text
//! cleanup_old_files()
//!
//! write.table(sample.stage, file = "sample_stages.txt", ...)
//! ```
//!
//! Here the report is read-only and deletion is a separate, explicit step. A
//! report that deletes data as a side effect is not one you can run to find out
//! what happened.

use std::collections::{HashMap, HashSet};
use std::hash::BuildHasher;
use std::io::BufRead;

/// A line of `logs/trim_stats.txt`: a file, a kind of count, and the count.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrimStat {
    /// The path as logged, `fastq/<sample>/<file>`.
    pub file: String,
    /// `total` or `discarded`.
    pub kind: String,
    pub reads: u64,
}

impl TrimStat {
    /// Parses one line.
    #[must_use]
    pub fn parse(line: &str) -> Option<Self> {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            return None;
        }
        Some(Self {
            file: f[0].to_string(),
            kind: f[1].to_string(),
            reads: f[2].trim().parse().ok()?,
        })
    }

    /// The sample the file belongs to: the directory under `fastq/`.
    #[must_use]
    pub fn sample(&self) -> Option<&str> {
        self.file.strip_prefix("fastq/")?.split('/').next()
    }
}

/// Keeps the last entry for each file and kind.
///
/// Upstream explains why: "duplicate entries are possible if the trim script
/// failed and had to be rerun; for these duplicates, take the most recent stats
/// (largest row number)". It does that by numbering the rows, reversing,
/// dropping duplicates and reversing back. Keeping the last occurrence and
/// preserving first-seen order is the same thing, said once.
#[must_use]
pub fn deduplicate(stats: Vec<TrimStat>) -> Vec<TrimStat> {
    let mut latest: HashMap<(String, String), TrimStat> = HashMap::new();
    let mut order: Vec<(String, String)> = Vec::new();
    for stat in stats {
        let key = (stat.file.clone(), stat.kind.clone());
        if !latest.contains_key(&key) {
            order.push(key.clone());
        }
        latest.insert(key, stat);
    }
    order
        .into_iter()
        .filter_map(|key| latest.remove(&key))
        .collect()
}

/// What trimming did to one sample.
#[derive(Debug, Clone, PartialEq)]
pub struct TrimSummary {
    pub sample: String,
    pub total_reads: u64,
    pub filtered_reads: u64,
    /// Files with a `total` line, which is how many were trimmed.
    pub files: usize,
}

impl TrimSummary {
    /// Reads left after trimming.
    #[must_use]
    pub const fn remaining_reads(&self) -> u64 {
        self.total_reads.saturating_sub(self.filtered_reads)
    }

    /// The fraction lost, or zero when nothing was read.
    #[must_use]
    pub fn percent_filtered(&self) -> f64 {
        if self.total_reads == 0 {
            return 0.0;
        }
        self.filtered_reads as f64 / self.total_reads as f64
    }
}

/// Summarises the trimming log per sample.
#[must_use]
pub fn summarise(stats: &[TrimStat]) -> Vec<TrimSummary> {
    let mut order: Vec<String> = Vec::new();
    let mut by_sample: HashMap<String, TrimSummary> = HashMap::new();

    for stat in stats {
        let Some(sample) = stat.sample() else {
            continue;
        };
        let entry = by_sample.entry(sample.to_string()).or_insert_with(|| {
            order.push(sample.to_string());
            TrimSummary {
                sample: sample.to_string(),
                total_reads: 0,
                filtered_reads: 0,
                files: 0,
            }
        });
        match stat.kind.as_str() {
            "total" => {
                entry.total_reads += stat.reads;
                entry.files += 1;
            }
            "discarded" => entry.filtered_reads += stat.reads,
            _ => {}
        }
    }

    order
        .into_iter()
        .filter_map(|s| by_sample.remove(&s))
        .collect()
}

/// How far a sample got.
///
/// The stages are cumulative: a sample at `Bed` passed everything before it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Stage {
    /// Nothing yet, or not all files arrived.
    None = 0,
    /// Every expected FASTQ is present.
    Fastq = 1,
    /// Trimming produced the expected counts.
    Trimmed = 2,
    /// The BAM passed its check and was cleaned up.
    Bam = 3,
    /// The BED passed its check and was cleaned up.
    Bed = 4,
}

/// What is known about one sample at the end of the report.
#[derive(Debug, Clone, PartialEq)]
pub struct SampleReport {
    pub sample: String,
    /// Position in the sample list, which is what the job arrays index by.
    pub index: Option<usize>,
    pub stage: Stage,
    pub total_reads: u64,
    pub remaining_reads: u64,
    pub percent_filtered: f64,
    /// Distinct aligned fragments, from the alignment report.
    pub aligned_fragments: Option<usize>,
}

impl SampleReport {
    /// Aligned fragments over reads surviving trimming.
    #[must_use]
    pub fn percent_aligned(&self) -> Option<f64> {
        let aligned = self.aligned_fragments?;
        if self.remaining_reads == 0 {
            return None;
        }
        Some(aligned as f64 / self.remaining_reads as f64)
    }

    /// Whether the sample looks wrong enough to look at.
    ///
    /// Upstream's two tests: fewer than 100 million reads left, or more than
    /// a tenth lost to trimming.
    #[must_use]
    pub fn has_quality_problem(&self) -> bool {
        self.remaining_reads < 100_000_000 || self.percent_filtered > 0.1
    }
}

/// Reads the `SAMPLES=( ... )` array out of a `config.sh`.
///
/// Upstream slices the file between the line matching `SAMPLES` and the first
/// line that is exactly `)`. The index it derives is what `LSB_JOBINDEX` maps
/// to, so it is one-based.
#[must_use]
pub fn parse_sample_list(config: &str) -> Vec<String> {
    let mut lines = config.lines().skip_while(|l| !l.contains("SAMPLES"));
    // The line naming the array is not itself a sample.
    lines.next();
    lines
        .take_while(|l| l.trim() != ")")
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .map(String::from)
        .collect()
}

/// Builds the report.
///
/// `expected` gives each sample's expected read count and file count, halved
/// from the index because it lists both mates. A sample whose counts match is
/// at least [`Stage::Trimmed`]; the alignment and BED stages come from the
/// clean-up log, since that is the only record that they passed.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn build<E: BuildHasher, A: BuildHasher, B: BuildHasher, D: BuildHasher>(
    samples: &[String],
    summaries: &[TrimSummary],
    expected: &HashMap<String, (u64, usize), E>,
    aligned: &HashMap<String, usize, A>,
    finished_bams: &HashSet<String, B>,
    finished_beds: &HashSet<String, D>,
) -> Vec<SampleReport> {
    let by_sample: HashMap<&str, &TrimSummary> =
        summaries.iter().map(|s| (s.sample.as_str(), s)).collect();

    // Every sample in the list, plus any that only appear in the log.
    let mut names: Vec<String> = samples.to_vec();
    for summary in summaries {
        if !names.contains(&summary.sample) {
            names.push(summary.sample.clone());
        }
    }

    names
        .into_iter()
        .map(|sample| {
            let summary = by_sample.get(sample.as_str());
            let (total_reads, remaining_reads, percent_filtered, files) =
                summary.map_or((0, 0, 0.0, 0), |s| {
                    (
                        s.total_reads,
                        s.remaining_reads(),
                        s.percent_filtered(),
                        s.files,
                    )
                });

            let mut stage = Stage::None;
            if let Some((expected_reads, expected_files)) = expected.get(&sample) {
                if files == *expected_files && files > 0 {
                    stage = Stage::Fastq;
                }
                if stage == Stage::Fastq && total_reads == *expected_reads {
                    stage = Stage::Trimmed;
                }
            }
            if finished_bams.contains(&sample) {
                stage = Stage::Bam;
            }
            if finished_beds.contains(&sample) {
                stage = Stage::Bed;
            }

            SampleReport {
                index: samples.iter().position(|s| *s == sample).map(|i| i + 1),
                sample: sample.clone(),
                stage,
                total_reads,
                remaining_reads,
                percent_filtered,
                aligned_fragments: aligned.get(&sample).copied(),
            }
        })
        .collect()
}

/// Reads a trimming log.
///
/// # Errors
/// Returns an error if the input cannot be read.
pub fn read_trim_stats<R: BufRead>(input: R) -> std::io::Result<Vec<TrimStat>> {
    let mut out = Vec::new();
    for line in input.lines() {
        if let Some(stat) = TrimStat::parse(&line?) {
            out.push(stat);
        }
    }
    Ok(deduplicate(out))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stat(file: &str, kind: &str, reads: u64) -> TrimStat {
        TrimStat {
            file: file.into(),
            kind: kind.into(),
            reads,
        }
    }

    #[test]
    fn the_sample_is_the_directory_under_fastq() {
        let s = stat("fastq/HG00250/ERR001_1.fastq.gz", "total", 10);
        assert_eq!(s.sample(), Some("HG00250"));
        assert_eq!(stat("elsewhere/x.gz", "total", 1).sample(), None);
    }

    /// A rerun appends rather than replacing, so the last entry wins.
    #[test]
    fn a_rerun_supersedes_the_earlier_entry() {
        let stats = vec![
            stat("fastq/S/a.gz", "total", 100),
            stat("fastq/S/a.gz", "discarded", 5),
            stat("fastq/S/a.gz", "total", 200),
        ];
        let deduped = deduplicate(stats);
        assert_eq!(deduped.len(), 2, "one total and one discarded");
        let total = deduped.iter().find(|s| s.kind == "total").unwrap();
        assert_eq!(total.reads, 200, "the later run wins");
        // And first-seen order is kept, so the file reads as it was written.
        assert_eq!(deduped[0].kind, "total");
    }

    #[test]
    fn summarising_adds_up_totals_and_counts_files() {
        let stats = vec![
            stat("fastq/S1/a_1.gz", "total", 100),
            stat("fastq/S1/a_1.gz", "discarded", 10),
            stat("fastq/S1/a_2.gz", "total", 100),
            stat("fastq/S1/a_2.gz", "discarded", 10),
            stat("fastq/S2/b_1.gz", "total", 50),
        ];
        let summaries = summarise(&stats);
        assert_eq!(summaries.len(), 2);

        let s1 = &summaries[0];
        assert_eq!(s1.sample, "S1");
        assert_eq!(s1.total_reads, 200);
        assert_eq!(s1.filtered_reads, 20);
        assert_eq!(s1.files, 2, "two files had a total line");
        assert_eq!(s1.remaining_reads(), 180);
        assert!((s1.percent_filtered() - 0.1).abs() < 1e-12);
    }

    #[test]
    fn a_sample_with_no_reads_has_no_percentage_rather_than_a_division_by_zero() {
        let summary = TrimSummary {
            sample: "S".into(),
            total_reads: 0,
            filtered_reads: 0,
            files: 0,
        };
        assert_eq!(summary.percent_filtered(), 0.0);
        assert_eq!(summary.remaining_reads(), 0);
    }

    /// The array in config.sh, sliced between the declaration and the closing
    /// parenthesis.
    #[test]
    fn the_sample_list_is_read_out_of_the_config() {
        let config = "#!/usr/bin/env bash\ngenome=x\n\nSAMPLES=(\nNA19914\nHG00623\nHG01139\n)\n";
        assert_eq!(parse_sample_list(config), ["NA19914", "HG00623", "HG01139"]);
    }

    #[test]
    fn an_absent_sample_list_reads_as_empty() {
        assert!(parse_sample_list("nothing here\n").is_empty());
    }

    /// The stage is cumulative, and the index is one-based because that is what
    /// the job arrays use.
    #[test]
    fn the_stage_follows_how_far_the_sample_got() {
        let samples = vec!["S1".to_string(), "S2".to_string()];
        let summaries = vec![
            TrimSummary {
                sample: "S1".into(),
                total_reads: 200,
                filtered_reads: 10,
                files: 2,
            },
            TrimSummary {
                sample: "S2".into(),
                total_reads: 100,
                filtered_reads: 5,
                files: 1,
            },
        ];
        let expected: HashMap<String, (u64, usize)> =
            [("S1".to_string(), (200, 2)), ("S2".to_string(), (200, 2))]
                .into_iter()
                .collect();

        let report = build(
            &samples,
            &summaries,
            &expected,
            &HashMap::new(),
            &HashSet::new(),
            &HashSet::new(),
        );

        assert_eq!(report[0].stage, Stage::Trimmed, "counts match");
        assert_eq!(report[0].index, Some(1), "one-based");
        assert_eq!(report[1].stage, Stage::None, "a file is missing");
    }

    #[test]
    fn the_clean_up_log_carries_the_later_stages() {
        let samples = vec!["S1".to_string()];
        let summaries = vec![TrimSummary {
            sample: "S1".into(),
            total_reads: 200,
            filtered_reads: 0,
            files: 2,
        }];
        let expected: HashMap<String, (u64, usize)> =
            [("S1".to_string(), (200, 2))].into_iter().collect();

        let bams: HashSet<String> = ["S1".to_string()].into_iter().collect();
        let report = build(
            &samples,
            &summaries,
            &expected,
            &HashMap::new(),
            &bams,
            &HashSet::new(),
        );
        assert_eq!(report[0].stage, Stage::Bam);

        let beds: HashSet<String> = ["S1".to_string()].into_iter().collect();
        let report = build(
            &samples,
            &summaries,
            &expected,
            &HashMap::new(),
            &bams,
            &beds,
        );
        assert_eq!(report[0].stage, Stage::Bed);
    }

    #[test]
    fn percent_aligned_divides_by_the_surviving_reads() {
        let mut r = SampleReport {
            sample: "S".into(),
            index: Some(1),
            stage: Stage::Bed,
            total_reads: 200,
            remaining_reads: 180,
            percent_filtered: 0.1,
            aligned_fragments: Some(90),
        };
        assert!((r.percent_aligned().unwrap() - 0.5).abs() < 1e-12);

        r.aligned_fragments = None;
        assert_eq!(r.percent_aligned(), None, "no report, no figure");

        r.aligned_fragments = Some(10);
        r.remaining_reads = 0;
        assert_eq!(r.percent_aligned(), None, "no division by zero");
    }

    /// The two tests upstream flags a sample on.
    #[test]
    fn quality_problems_are_shallow_or_over_trimmed_samples() {
        let base = SampleReport {
            sample: "S".into(),
            index: Some(1),
            stage: Stage::Bed,
            total_reads: 200_000_000,
            remaining_reads: 190_000_000,
            percent_filtered: 0.05,
            aligned_fragments: None,
        };
        assert!(!base.has_quality_problem());

        let shallow = SampleReport {
            remaining_reads: 99_000_000,
            ..base.clone()
        };
        assert!(shallow.has_quality_problem(), "too few reads left");

        let trimmed = SampleReport {
            percent_filtered: 0.11,
            ..base
        };
        assert!(trimmed.has_quality_problem(), "too much lost to trimming");
    }

    /// A sample present in the log but absent from the list still appears, and
    /// has no index because the job arrays cannot reach it.
    #[test]
    fn a_sample_outside_the_list_still_appears() {
        let summaries = vec![TrimSummary {
            sample: "stray".into(),
            total_reads: 10,
            filtered_reads: 0,
            files: 1,
        }];
        let report = build(
            &[],
            &summaries,
            &HashMap::new(),
            &HashMap::new(),
            &HashSet::new(),
            &HashSet::new(),
        );
        assert_eq!(report.len(), 1);
        assert_eq!(report[0].sample, "stray");
        assert_eq!(report[0].index, None);
    }

    #[test]
    fn reading_a_log_deduplicates_as_it_goes() {
        let log = "fastq/S/a.gz\ttotal\t100\nfastq/S/a.gz\ttotal\t200\nbroken line\n";
        let stats = read_trim_stats(log.as_bytes()).unwrap();
        assert_eq!(stats.len(), 1);
        assert_eq!(stats[0].reads, 200);
    }
}
