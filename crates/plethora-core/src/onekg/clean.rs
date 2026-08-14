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
    let mut removals = Vec::new();

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
            removals.push(Removal::Fastq);
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
            removals.push(Removal::Bam);
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
        if remove_fastq && counts.fastq.is_some() && !removals.contains(&Removal::Fastq) {
            removals.push(Removal::Fastq);
        }
        if counts.bam.is_some() && !removals.contains(&Removal::Bam) {
            removals.push(Removal::Bam);
        }
        if counts.sorted_bam.is_some() {
            removals.push(Removal::SortedBam);
        }
        // The unsorted BED goes once the sorted one exists; the caller checks
        // for it, since the count says nothing about whether it was written.
        removals.push(Removal::Bed);
    }

    Ok(removals)
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
