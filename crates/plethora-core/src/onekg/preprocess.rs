//! `preprocessing_1000genomes.R`: choose which public samples to process.
//!
//! The full 1000 Genomes index lists more than a thousand usable samples, which
//! is more data than the pipeline can hold at once. This picks a few hundred:
//! everything already studied or already on hand, then a quota from each
//! population so the mean and variance stay representative.
//!
//! The upstream script does not run. Its final `write.table` closes `paste0`
//! after the last argument instead of after the filename, so `sep`, `quote` and
//! `row.names` are swallowed and `write.table(` is never closed; the file fails
//! to parse. Reported as dpastling/plethora#15. The selection below is what the
//! rest of the script says it does.

use std::collections::{HashMap, HashSet};

use super::sample_index::Record;

/// Samples wanted from each population.
pub const POPULATION_QUOTA: usize = 25;

/// The genome size the coverage estimate divides by.
pub const GENOME_SIZE: f64 = 3.235e9;

/// Sequencing centres in preference order.
///
/// The comment upstream explains why: BCM's libraries had a short insert size
/// and BI's had low overall qualities, so both are ranked below the rest. The
/// first five are deliberately not ordered among themselves.
pub const CENTRE_PREFERENCE: [&str; 7] = ["BGI", "SC", "ILLUMINA", "MPIMG", "WUGSC", "BCM", "BI"];

/// The populations the study keeps.
pub const POPULATIONS: [&str; 14] = [
    "MXL", "CLM", "PUR", "ASW", "LWK", "YRI", "JPT", "CHB", "CHS", "TSI", "CEU", "IBS", "FIN",
    "GBR",
];

/// One candidate: a sample as sequenced by one centre from one library.
#[derive(Debug, Clone, PartialEq)]
pub struct Candidate {
    pub sample_name: String,
    pub centre_name: String,
    pub library_name: String,
    pub population: String,
    /// Files contributing, counting the second mate only.
    pub files: usize,
    pub reads: u64,
    pub mean_insert_size: f64,
}

impl Candidate {
    /// Estimated coverage: `reads * insert size / genome size`.
    ///
    /// The comment upstream reads "10x coverage is 100 million reads for an
    /// insert size of 300bp", which is what this arithmetic encodes. It counts
    /// fragment length rather than sequenced bases, so it is not coverage in
    /// the usual sense; it is the number the thresholds were chosen against.
    #[must_use]
    pub fn coverage(&self) -> f64 {
        (self.reads as f64 * self.mean_insert_size) / GENOME_SIZE
    }

    /// Where the centre sits in the preference order.
    ///
    /// `as.numeric(factor(CENTER_NAME, levels = ...))` gives `NA` for a centre
    /// that is not listed, and `ifelse(center.rank < 6, 1, center.rank)` leaves
    /// that `NA` alone, so an unlisted centre sorts last rather than first.
    /// The first five collapse to one because their order does not matter.
    #[must_use]
    pub fn centre_rank(&self) -> Option<u8> {
        let position = CENTRE_PREFERENCE
            .iter()
            .position(|c| *c == self.centre_name)?;
        let rank = u8::try_from(position + 1).ok()?;
        Some(if rank < 6 { 1 } else { rank })
    }
}

/// Which rows survive the first filter.
///
/// `HiSeq`, Illumina, paired, not withdrawn, and not an exome study. The exome
/// test is case-insensitive upstream, and so is this.
#[must_use]
pub fn is_usable(record: &Record) -> bool {
    record.instrument_model.contains("HiSeq")
        && record.instrument_platform == "ILLUMINA"
        && record.library_layout == "PAIRED"
        && !record.withdrawn
        && !record.study_name.to_lowercase().contains("exome")
}

/// Groups usable rows into candidates.
///
/// Only the second mate of each pair is counted, which is how upstream avoids
/// double counting: `filter(grepl("_2.(filt.)*fastq.gz", FASTQ_FILE))`. The
/// read count it then sums is therefore per fragment, not per read.
#[must_use]
pub fn candidates(records: &[Record]) -> Vec<Candidate> {
    // Insertion order is kept so the result does not depend on hashing.
    let mut order: Vec<(String, String, String)> = Vec::new();
    let mut groups: HashMap<(String, String, String), Vec<&Record>> = HashMap::new();

    for record in records {
        if !is_usable(record) || record.mate() != Some(2) {
            continue;
        }
        let key = (
            record.sample_name.clone(),
            record.centre_name(),
            record.library_name.clone(),
        );
        if !groups.contains_key(&key) {
            order.push(key.clone());
        }
        groups.entry(key).or_default().push(record);
    }

    order
        .into_iter()
        .map(|key| {
            let rows = &groups[&key];
            let inserts: Vec<f64> = rows
                .iter()
                .filter_map(|r| r.insert_size)
                .map(|v| v as f64)
                .collect();
            Candidate {
                sample_name: key.0,
                centre_name: key.1,
                library_name: key.2,
                population: rows[0].population.clone(),
                files: rows.len(),
                reads: rows.iter().filter_map(|r| r.read_count).sum(),
                mean_insert_size: if inserts.is_empty() {
                    // `mean` of nothing is NaN in R too, and the coverage
                    // filter then rejects the candidate.
                    f64::NAN
                } else {
                    inserts.iter().sum::<f64>() / inserts.len() as f64
                },
            }
        })
        .collect()
}

impl Record {
    /// The sequencing centre, named so the grouping reads as it does upstream.
    #[must_use]
    pub fn centre_name(&self) -> String {
        self.center_name.clone()
    }
}

/// Keeps one candidate per sample: the best-ranked centre, then the deepest.
///
/// `arrange(center.rank, desc(n.reads))` followed by taking the first is a
/// stable sort in dplyr, so candidates tying on both keep their input order.
#[must_use]
pub fn best_per_sample(mut candidates: Vec<Candidate>) -> Vec<Candidate> {
    // Sorted stably, exactly as arrange does, then the first of each sample.
    candidates.sort_by(|a, b| {
        // An unranked centre sorts after every ranked one, as NA does.
        let rank = |c: &Candidate| c.centre_rank().unwrap_or(u8::MAX);
        rank(a).cmp(&rank(b)).then_with(|| b.reads.cmp(&a.reads))
    });

    let mut seen: HashSet<String> = HashSet::new();
    candidates
        .into_iter()
        .filter(|c| seen.insert(c.sample_name.clone()))
        .collect()
}

/// The lists that steer the selection.
#[derive(Debug, Clone, Default)]
pub struct Lists {
    /// Samples that failed QC before and are never chosen.
    pub failed: HashSet<String>,
    /// Samples from Sudmant et al. 2010, always chosen.
    pub sudmant: HashSet<String>,
    /// Samples the lab holds DNA for, always chosen.
    pub dna: HashSet<String>,
    /// Samples with Irys data, always chosen.
    pub irys: HashSet<String>,
}

impl Lists {
    /// Samples taken regardless of quota.
    #[must_use]
    pub fn of_interest(&self) -> HashSet<String> {
        self.sudmant
            .union(&self.dna)
            .chain(self.irys.iter())
            .cloned()
            .collect()
    }
}

/// The thresholds the second pass applies.
///
/// Upstream filters twice with different numbers, and only the second set
/// decides anything: the first pass's `coverage > 10` is subsumed by the
/// second's `coverage >= 12`. Both are kept because the first also drops
/// samples before the per-sample ranking, which changes which library wins.
#[derive(Debug, Clone, Copy)]
pub struct Thresholds {
    pub min_coverage: f64,
    pub min_reads: u64,
    pub max_reads: u64,
}

impl Default for Thresholds {
    fn default() -> Self {
        // "we want the unfiltered coverage to be 12 so allow for up to 10% of
        // the reads to be trimmed"
        Self {
            min_coverage: 12.0,
            min_reads: 90_000_000,
            max_reads: 400_000_000,
        }
    }
}

impl Thresholds {
    /// The first pass's looser set.
    #[must_use]
    pub fn first_pass() -> Self {
        Self {
            min_coverage: 10.0,
            min_reads: 0,
            max_reads: 400_000_000,
        }
    }

    /// True when a candidate clears them.
    #[must_use]
    pub fn admits(&self, candidate: &Candidate) -> bool {
        let coverage = candidate.coverage();
        coverage.is_finite()
            && coverage >= self.min_coverage
            && candidate.reads >= self.min_reads
            && candidate.reads < self.max_reads
    }
}

/// Runs the selection, returning the chosen sample names in order.
///
/// Samples of interest come first, then each population is topped up to the
/// quota from the deepest remaining candidates, best-ranked centre first.
#[must_use]
pub fn select(records: &[Record], lists: &Lists) -> Vec<String> {
    // First pass: the looser thresholds, applied before the per-sample choice,
    // because dropping a library here can change which one wins.
    let first = Thresholds::first_pass();
    let surviving: Vec<Candidate> = candidates(records)
        .into_iter()
        .filter(|c| first.admits(c) && !lists.failed.contains(&c.sample_name))
        .collect();

    let second = Thresholds::default();
    let mut pool: Vec<Candidate> = best_per_sample(surviving)
        .into_iter()
        .filter(|c| second.admits(c) && !lists.failed.contains(&c.sample_name))
        .collect();

    let of_interest = lists.of_interest();
    let mut chosen: Vec<String> = Vec::new();
    let mut taken: HashSet<String> = HashSet::new();

    // Everything already studied or already on hand.
    for candidate in &pool {
        if of_interest.contains(&candidate.sample_name)
            && taken.insert(candidate.sample_name.clone())
        {
            chosen.push(candidate.sample_name.clone());
        }
    }

    // How much of each population that already covers.
    let mut filled: HashMap<String, usize> = HashMap::new();
    for candidate in &pool {
        if taken.contains(&candidate.sample_name) {
            *filled.entry(candidate.population.clone()).or_default() += 1;
        }
    }

    // Then top up, in the order the ranking already put them.
    pool.retain(|c| {
        !taken.contains(&c.sample_name) && POPULATIONS.contains(&c.population.as_str())
    });
    for candidate in &pool {
        let already = filled.entry(candidate.population.clone()).or_default();
        if *already >= POPULATION_QUOTA {
            continue;
        }
        if taken.insert(candidate.sample_name.clone()) {
            *already += 1;
            chosen.push(candidate.sample_name.clone());
        }
    }

    chosen
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(
        sample: &str,
        centre: &str,
        library: &str,
        population: &str,
        reads: u64,
        insert: u64,
        mate: u8,
    ) -> Record {
        let mut f: Vec<String> = vec![String::new(); 26];
        f[0] = format!("ftp://x/{sample}_{library}_{mate}.fastq.gz");
        f[4] = "1000 Genomes GBR population sequencing".into();
        f[5] = centre.into();
        f[9] = sample.into();
        f[10] = population.into();
        f[12] = "ILLUMINA".into();
        f[13] = "Illumina HiSeq 2000".into();
        f[14] = library.into();
        f[17] = insert.to_string();
        f[18] = "PAIRED".into();
        f[20] = "0".into();
        f[23] = reads.to_string();
        Record::parse(&f.join("\t")).expect("a well-formed row")
    }

    #[test]
    fn the_first_filter_keeps_only_paired_hiseq_illumina() {
        let good = record("S1", "BGI", "L1", "GBR", 100, 300, 2);
        assert!(is_usable(&good));

        let mut f: Vec<String> = vec![String::new(); 26];
        f[12] = "ILLUMINA".into();
        f[13] = "Illumina GAIIx".into();
        f[18] = "PAIRED".into();
        f[20] = "0".into();
        assert!(
            !is_usable(&Record::parse(&f.join("\t")).unwrap()),
            "not HiSeq"
        );
    }

    /// An exome study is dropped, and the test is case-insensitive.
    #[test]
    fn exome_studies_are_dropped_whatever_the_case() {
        let mut r = record("S1", "BGI", "L1", "GBR", 100, 300, 2);
        r.study_name = "1000 Genomes Exome Project".into();
        assert!(!is_usable(&r));
        r.study_name = "1000 genomes exome".into();
        assert!(!is_usable(&r));
    }

    #[test]
    fn a_withdrawn_run_is_dropped() {
        let mut r = record("S1", "BGI", "L1", "GBR", 100, 300, 2);
        r.withdrawn = true;
        assert!(!is_usable(&r));
    }

    /// Only the second mate is counted, so a pair contributes its reads once.
    #[test]
    fn only_the_second_mate_is_counted() {
        let records = vec![
            record("S1", "BGI", "L1", "GBR", 100_000_000, 300, 1),
            record("S1", "BGI", "L1", "GBR", 100_000_000, 300, 2),
        ];
        let c = candidates(&records);
        assert_eq!(c.len(), 1);
        assert_eq!(c[0].files, 1, "one file, not two");
        assert_eq!(c[0].reads, 100_000_000);
    }

    #[test]
    fn candidates_group_by_sample_centre_and_library() {
        let records = vec![
            record("S1", "BGI", "L1", "GBR", 50_000_000, 300, 2),
            record("S1", "BGI", "L1", "GBR", 50_000_000, 300, 2),
            record("S1", "BGI", "L2", "GBR", 90_000_000, 300, 2),
            record("S1", "BI", "L1", "GBR", 90_000_000, 300, 2),
        ];
        let c = candidates(&records);
        assert_eq!(c.len(), 3, "three distinct (sample, centre, library)");
        assert_eq!(c[0].reads, 100_000_000, "the two rows of L1 summed");
        assert_eq!(c[0].files, 2);
    }

    /// The five preferred centres collapse to one rank; BCM and BI keep theirs;
    /// an unlisted centre has none and sorts last.
    #[test]
    fn centre_ranks_follow_the_stated_preference() {
        let rank = |centre: &str| record("S", centre, "L", "GBR", 1, 1, 2);
        let c = candidates(&[rank("BGI")]);
        assert_eq!(c[0].centre_rank(), Some(1));
        assert_eq!(candidates(&[rank("WUGSC")])[0].centre_rank(), Some(1));
        assert_eq!(candidates(&[rank("BCM")])[0].centre_rank(), Some(6));
        assert_eq!(candidates(&[rank("BI")])[0].centre_rank(), Some(7));
        assert_eq!(candidates(&[rank("SOMEWHERE")])[0].centre_rank(), None);
    }

    #[test]
    fn coverage_is_reads_times_insert_over_the_genome() {
        let c = &candidates(&[record("S1", "BGI", "L1", "GBR", 100_000_000, 300, 2)])[0];
        // The comment upstream: 100 million reads at 300 bp is about 10x.
        assert!((c.coverage() - 9.27).abs() < 0.01, "got {}", c.coverage());
    }

    /// The best library per sample: preferred centre first, then depth.
    #[test]
    fn the_best_library_per_sample_wins() {
        let pool = candidates(&[
            record("S1", "BI", "deep", "GBR", 300_000_000, 400, 2),
            record("S1", "BGI", "shallow", "GBR", 100_000_000, 400, 2),
            record("S1", "BGI", "deeper", "GBR", 200_000_000, 400, 2),
        ]);
        let best = best_per_sample(pool);
        assert_eq!(best.len(), 1);
        assert_eq!(
            best[0].library_name, "deeper",
            "BGI beats BI, then depth decides"
        );
    }

    #[test]
    fn an_unranked_centre_loses_to_a_ranked_one() {
        let pool = candidates(&[
            record("S1", "SOMEWHERE", "deep", "GBR", 400_000_000, 400, 2),
            record("S1", "BI", "shallow", "GBR", 100_000_000, 400, 2),
        ]);
        assert_eq!(best_per_sample(pool)[0].library_name, "shallow");
    }

    /// Samples on the lists are taken whatever their population's quota.
    #[test]
    fn samples_of_interest_are_taken_first() {
        let records: Vec<Record> = (0..5)
            .map(|i| record(&format!("S{i}"), "BGI", "L", "GBR", 150_000_000, 400, 2))
            .collect();
        let mut lists = Lists::default();
        lists.sudmant.insert("S3".into());

        let chosen = select(&records, &lists);
        assert_eq!(chosen[0], "S3", "the sample of interest leads");
        assert_eq!(chosen.len(), 5);
    }

    #[test]
    fn a_failed_sample_is_never_chosen() {
        let records: Vec<Record> = (0..3)
            .map(|i| record(&format!("S{i}"), "BGI", "L", "GBR", 150_000_000, 400, 2))
            .collect();
        let mut lists = Lists::default();
        lists.failed.insert("S1".into());
        // Even being of interest does not override it.
        lists.sudmant.insert("S1".into());

        let chosen = select(&records, &lists);
        assert!(!chosen.contains(&"S1".to_string()));
        assert_eq!(chosen.len(), 2);
    }

    /// The quota caps each population, and the deepest are kept.
    #[test]
    fn each_population_is_capped_at_the_quota() {
        let records: Vec<Record> = (0..40)
            .map(|i| {
                record(
                    &format!("S{i:02}"),
                    "BGI",
                    "L",
                    "GBR",
                    100_000_000 + u64::from(i as u32) * 1_000_000,
                    400,
                    2,
                )
            })
            .collect();
        let chosen = select(&records, &Lists::default());
        assert_eq!(chosen.len(), POPULATION_QUOTA);
        // Deepest first, so the highest-numbered samples are the ones kept.
        assert_eq!(chosen[0], "S39");
    }

    #[test]
    fn a_population_outside_the_study_is_not_topped_up() {
        let records: Vec<Record> = (0..5)
            .map(|i| record(&format!("S{i}"), "BGI", "L", "PEL", 150_000_000, 400, 2))
            .collect();
        assert!(
            select(&records, &Lists::default()).is_empty(),
            "PEL is not among the fourteen populations"
        );
    }

    /// A sample below the coverage floor is dropped even if nothing else
    /// competes for its quota.
    #[test]
    fn the_thresholds_drop_shallow_samples() {
        let shallow = record("S1", "BGI", "L", "GBR", 90_000_000, 100, 2);
        assert!(shallow_coverage(&shallow) < 12.0);
        assert!(select(&[shallow], &Lists::default()).is_empty());
    }

    fn shallow_coverage(r: &Record) -> f64 {
        candidates(std::slice::from_ref(r))[0].coverage()
    }

    #[test]
    fn a_candidate_with_no_insert_size_is_rejected_rather_than_guessed() {
        let mut r = record("S1", "BGI", "L", "GBR", 150_000_000, 400, 2);
        r.insert_size = None;
        let c = &candidates(std::slice::from_ref(&r))[0];
        assert!(c.coverage().is_nan());
        assert!(!Thresholds::default().admits(c));
    }
}
