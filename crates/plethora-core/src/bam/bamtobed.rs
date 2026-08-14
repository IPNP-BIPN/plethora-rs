//! `bedtools bamtobed`, in the two forms `make_bed.sh` uses.
//!
//! ```text
//! bedtools bamtobed -split -bedpe -i ${output}_sorted.bam > $output.bed   # paired
//! bedtools bamtobed -i $bam > ${output}_edited.bed                        # single
//! ```
//!
//! Transcribed from bedtools' `src/bamToBed/bamToBed.cpp`: [`bedpe`] from
//! `PrintBedPE`, [`bed`] from `PrintBed`, and [`BedpeIter`] from
//! `ConvertBamToBedpe`.
//!
//! Three things about that source are worth stating, because none of them are
//! guessable from the command line:
//!
//! - `-split` is silently ignored when `-bedpe` is given. A read with an N in
//!   its CIGAR spans the gap either way. Upstream passes both flags; only
//!   `-bedpe` has any effect.
//! - The two blocks are ordered by comparing the chromosome *name as a string*,
//!   not the reference id, and only then by position. Since `.` (0x2E) sorts
//!   below every real chromosome name, an unmapped mate always lands in the
//!   first block. So the first block is not read 1: for a pair whose read 1 is
//!   downstream, the blocks come out swapped.
//! - `ConvertBamToBedpe` consumes records strictly two at a time and pairs them
//!   only if their names match. A third record under the same name, which is
//!   what BWA-MEM produces for a supplementary alignment, pushes every later
//!   record out of phase and makes bedtools skip pairs in a cascade. See
//!   [`is_pairable`].

use std::fmt;

/// SAM flag: the segment is part of a pair.
pub const PAIRED: u16 = 0x1;
/// SAM flag: the segment is unmapped.
pub const UNMAPPED: u16 = 0x4;
/// SAM flag: the segment is on the reverse strand.
pub const REVERSE: u16 = 0x10;
/// SAM flag: the first segment of the pair.
pub const FIRST_MATE: u16 = 0x40;
/// SAM flag: the last segment of the pair.
pub const SECOND_MATE: u16 = 0x80;
/// SAM flag: a secondary alignment.
pub const SECONDARY: u16 = 0x100;
/// SAM flag: a supplementary alignment.
pub const SUPPLEMENTARY: u16 = 0x800;

/// The fields of an alignment that `bamtobed` reads.
///
/// A narrow view rather than a BAM record, so the bedtools rules can be tested
/// without constructing BAM files.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Aln {
    /// QNAME.
    pub name: Vec<u8>,
    /// FLAG.
    pub flags: u16,
    /// Reference name, or `None` when unmapped.
    pub chrom: Option<String>,
    /// Zero-based leftmost position; meaningless when unmapped.
    pub start: i64,
    /// Half-open end, from the reference-consuming CIGAR operations only.
    /// Soft and hard clips are excluded, so a soft-clipped read spans less.
    pub end: i64,
    /// MAPQ.
    pub mapq: u8,
}

impl Aln {
    /// True when the record has a position on a reference.
    #[must_use]
    pub const fn is_mapped(&self) -> bool {
        self.flags & UNMAPPED == 0
    }

    #[must_use]
    const fn is_paired(&self) -> bool {
        self.flags & PAIRED != 0
    }

    #[must_use]
    const fn strand(&self) -> char {
        if self.flags & REVERSE == 0 { '+' } else { '-' }
    }
}

/// True when a record can take part in the BEDPE pairing.
///
/// bedtools has no such filter, which is the problem: it assumes exactly two
/// records per name. Under bowtie2 that holds, so the paper's pipeline never
/// met the failure. BWA-MEM emits supplementary and secondary records, and a
/// third record under one name desynchronises the two-at-a-time reader for the
/// whole rest of the file. Filtering here is a deliberate divergence, recorded
/// in `DIVERGENCES.md`, and it is what makes the BWA-MEM path usable at all.
#[must_use]
pub const fn is_pairable(flags: u16) -> bool {
    flags & (SECONDARY | SUPPLEMENTARY) == 0
}

/// A BED6 line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bed6 {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub name: String,
    pub score: u8,
    pub strand: char,
}

impl fmt::Display for Bed6 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.chrom, self.start, self.end, self.name, self.score, self.strand
        )
    }
}

/// A BEDPE line, with `.` and -1 standing in for an unmapped end.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BedPe {
    pub chrom1: String,
    pub start1: i64,
    pub end1: i64,
    pub chrom2: String,
    pub start2: i64,
    pub end2: i64,
    pub name: String,
    pub score: u8,
    pub strand1: char,
    pub strand2: char,
}

impl fmt::Display for BedPe {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.chrom1,
            self.start1,
            self.end1,
            self.chrom2,
            self.start2,
            self.end2,
            self.name,
            self.score,
            self.strand1,
            self.strand2
        )
    }
}

/// `PrintBed`: one BED6 line per mapped record.
///
/// Unmapped records produce nothing, which is why the single-end file has fewer
/// lines than the BAM has records. The `/1` and `/2` suffixes come from the
/// pair flags; a record with both set gets `/1/2`, as in the C++.
#[must_use]
pub fn bed(a: &Aln) -> Option<Bed6> {
    if !a.is_mapped() {
        return None;
    }

    let mut name = String::from_utf8_lossy(&a.name).into_owned();
    if a.flags & FIRST_MATE != 0 {
        name.push_str("/1");
    }
    if a.flags & SECOND_MATE != 0 {
        name.push_str("/2");
    }

    Some(Bed6 {
        chrom: a.chrom.clone().unwrap_or_else(|| ".".to_string()),
        start: a.start,
        end: a.end,
        name,
        score: a.mapq,
        strand: a.strand(),
    })
}

/// `PrintBedPE`: one BEDPE line per pair.
///
/// The name is taken from `bam1` before any swap, so it does not follow the
/// block that ends up first.
#[must_use]
pub fn bedpe(bam1: &Aln, bam2: &Aln) -> BedPe {
    let unset = |a: &Aln| -> (String, i64, i64, char) {
        if a.is_mapped() {
            (
                a.chrom.clone().unwrap_or_else(|| ".".to_string()),
                a.start,
                a.end,
                a.strand(),
            )
        } else {
            (".".to_string(), -1, -1, '.')
        }
    };

    let (mut chrom1, mut start1, mut end1, mut strand1) = unset(bam1);
    let (mut chrom2, mut start2, mut end2, mut strand2) = unset(bam2);

    // Ordered by chromosome name as a string, then by position. Not by
    // reference id: "chr10" sorts before "chr9" here, and "." before both.
    if chrom1 > chrom2 || (chrom1 == chrom2 && start1 > start2) {
        std::mem::swap(&mut chrom1, &mut chrom2);
        std::mem::swap(&mut start1, &mut start2);
        std::mem::swap(&mut end1, &mut end2);
        std::mem::swap(&mut strand1, &mut strand2);
    }

    // Zero unless both ends are mapped, even if one end has a high MAPQ.
    let score = if bam1.is_mapped() && bam2.is_mapped() {
        bam1.mapq.min(bam2.mapq)
    } else {
        0
    };

    BedPe {
        chrom1,
        start1,
        end1,
        chrom2,
        start2,
        end2,
        name: String::from_utf8_lossy(&bam1.name).into_owned(),
        score,
        strand1,
        strand2,
    }
}

/// `ConvertBamToBedpe`: pair up a name-grouped stream, two records at a time.
///
/// Reproduces the skip-forward recovery of the original, including the fact
/// that it emits the pair it recovers on. Records whose mate does not sit next
/// to them are reported through `on_orphan` rather than printed to stderr, so
/// the caller decides what a warning means.
pub struct BedpeIter<I: Iterator<Item = Aln>> {
    inner: I,
    orphans: usize,
}

impl<I: Iterator<Item = Aln>> BedpeIter<I> {
    /// Wraps a name-grouped record stream.
    pub fn new(inner: I) -> Self {
        Self { inner, orphans: 0 }
    }

    /// How many paired records were skipped for want of an adjacent mate.
    ///
    /// bedtools prints one warning line per skipped record. A non-zero count
    /// here on a bowtie2 alignment means the file was not name-sorted; on a
    /// BWA-MEM alignment it means supplementary records were not filtered.
    #[must_use]
    pub const fn orphans(&self) -> usize {
        self.orphans
    }
}

impl<I: Iterator<Item = Aln>> Iterator for BedpeIter<I> {
    type Item = BedPe;

    fn next(&mut self) -> Option<BedPe> {
        let mut bam1 = self.inner.next()?;
        let Some(mut bam2) = self.inner.next() else {
            // The C++ leaves bam2 holding its previous value at end of file and
            // pairs against it. Stopping is the one place this deliberately
            // does not reproduce bedtools, because reproducing it would emit a
            // record built from a stale buffer.
            return None;
        };

        if bam1.name != bam2.name {
            while bam1.name != bam2.name {
                if bam1.is_paired() {
                    self.orphans += 1;
                }
                bam1 = bam2;
                bam2 = self.inner.next()?;
            }
            return Some(bedpe(&bam1, &bam2));
        }

        if bam1.is_paired() && bam2.is_paired() {
            return Some(bedpe(&bam1, &bam2));
        }

        // Two records sharing a name but not flagged as paired produce nothing,
        // and both are consumed.
        self.next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln(name: &str, flags: u16, chrom: Option<&str>, start: i64, end: i64, mapq: u8) -> Aln {
        Aln {
            name: name.as_bytes().to_vec(),
            flags,
            chrom: chrom.map(String::from),
            start,
            end,
            mapq,
        }
    }

    /// The plain case, checked against `bedtools bamtobed -bedpe`.
    #[test]
    fn a_proper_pair_keeps_its_order() {
        let r1 = aln("r1", PAIRED | FIRST_MATE, Some("chr1"), 99, 149, 60);
        let r2 = aln("r1", PAIRED | SECOND_MATE | REVERSE, Some("chr1"), 299, 349, 60);
        assert_eq!(
            bedpe(&r1, &r2).to_string(),
            "chr1\t99\t149\tchr1\t299\t349\tr1\t60\t+\t-"
        );
    }

    /// Read 1 downstream of read 2: the blocks swap, so the first block is
    /// read 2 and carries read 2's strand.
    #[test]
    fn blocks_are_ordered_by_position_not_by_mate() {
        let r1 = aln("c", PAIRED | FIRST_MATE | REVERSE, Some("chr1"), 899, 949, 60);
        let r2 = aln("c", PAIRED | SECOND_MATE, Some("chr1"), 699, 749, 60);
        assert_eq!(
            bedpe(&r1, &r2).to_string(),
            "chr1\t699\t749\tchr1\t899\t949\tc\t60\t+\t-"
        );
    }

    /// An unmapped mate sorts first, because "." precedes every chromosome name.
    #[test]
    fn an_unmapped_mate_lands_in_the_first_block() {
        let r1 = aln("a", PAIRED | UNMAPPED | FIRST_MATE, None, 0, 0, 0);
        let r2 = aln("a", PAIRED | SECOND_MATE, Some("chr1"), 499, 549, 44);
        assert_eq!(
            bedpe(&r1, &r2).to_string(),
            ".\t-1\t-1\tchr1\t499\t549\ta\t0\t.\t+"
        );

        // And the other way round: still first, even though it is read 2.
        let r1 = aln("a", PAIRED | FIRST_MATE, Some("chr1"), 499, 549, 44);
        let r2 = aln("a", PAIRED | UNMAPPED | SECOND_MATE, None, 0, 0, 0);
        assert_eq!(
            bedpe(&r1, &r2).to_string(),
            ".\t-1\t-1\tchr1\t499\t549\ta\t0\t.\t+"
        );
    }

    /// The score is the smaller MAPQ, and zero as soon as one end is unmapped.
    #[test]
    fn the_score_is_the_minimum_mapping_quality() {
        let r1 = aln("b", PAIRED | FIRST_MATE, Some("chr1"), 99, 149, 60);
        let r2 = aln("b", PAIRED | SECOND_MATE | REVERSE, Some("chr1"), 299, 349, 10);
        assert_eq!(bedpe(&r1, &r2).score, 10);

        let r2 = aln("b", PAIRED | UNMAPPED | SECOND_MATE, None, 0, 0, 0);
        assert_eq!(bedpe(&r1, &r2).score, 0);
    }

    /// Chromosomes are compared as strings, so "chr10" precedes "chr2".
    #[test]
    fn chromosomes_compare_as_strings() {
        let r1 = aln("d", PAIRED | FIRST_MATE, Some("chr2"), 100, 150, 60);
        let r2 = aln("d", PAIRED | SECOND_MATE, Some("chr10"), 100, 150, 60);
        let out = bedpe(&r1, &r2);
        assert_eq!(out.chrom1, "chr10", "string order, not numeric or reference id");
        assert_eq!(out.chrom2, "chr2");
    }

    #[test]
    fn both_ends_unmapped_gives_an_empty_line() {
        let r1 = aln("e", PAIRED | UNMAPPED | FIRST_MATE, None, 0, 0, 0);
        let r2 = aln("e", PAIRED | UNMAPPED | SECOND_MATE, None, 0, 0, 0);
        assert_eq!(bedpe(&r1, &r2).to_string(), ".\t-1\t-1\t.\t-1\t-1\te\t0\t.\t.");
    }

    #[test]
    fn single_end_mode_suffixes_the_name_and_drops_unmapped() {
        let r1 = aln("r1", PAIRED | FIRST_MATE, Some("chr1"), 99, 149, 60);
        assert_eq!(bed(&r1).unwrap().to_string(), "chr1\t99\t149\tr1/1\t60\t+");

        let r2 = aln("r1", PAIRED | SECOND_MATE | REVERSE, Some("chr1"), 299, 349, 60);
        assert_eq!(bed(&r2).unwrap().to_string(), "chr1\t299\t349\tr1/2\t60\t-");

        let u = aln("r3", PAIRED | UNMAPPED | SECOND_MATE, None, 0, 0, 0);
        assert!(bed(&u).is_none());
    }

    /// Both pair bits set is not a state an aligner should produce, but the C++
    /// appends both suffixes rather than choosing, so this does too.
    #[test]
    fn both_pair_bits_append_both_suffixes() {
        let r = aln("x", PAIRED | FIRST_MATE | SECOND_MATE, Some("chr1"), 0, 50, 60);
        assert_eq!(bed(&r).unwrap().name, "x/1/2");
    }

    #[test]
    fn the_iterator_pairs_adjacent_records() {
        let records = vec![
            aln("r1", PAIRED | FIRST_MATE, Some("chr1"), 99, 149, 60),
            aln("r1", PAIRED | SECOND_MATE | REVERSE, Some("chr1"), 299, 349, 60),
            aln("r2", PAIRED | FIRST_MATE, Some("chr1"), 499, 549, 60),
            aln("r2", PAIRED | SECOND_MATE | REVERSE, Some("chr1"), 899, 949, 60),
        ];
        let mut it = BedpeIter::new(records.into_iter());
        assert_eq!(it.next().unwrap().name, "r1");
        assert_eq!(it.next().unwrap().name, "r2");
        assert!(it.next().is_none());
        assert_eq!(it.orphans(), 0);
    }

    /// A third record under one name desynchronises the reader: this is the
    /// bedtools behaviour that makes filtering supplementary records necessary.
    #[test]
    fn an_extra_record_desynchronises_the_pairing() {
        let records = vec![
            aln("r1", PAIRED | FIRST_MATE, Some("chr1"), 99, 149, 60),
            aln("r1", PAIRED | SECOND_MATE, Some("chr1"), 299, 349, 60),
            aln("r1", PAIRED | SUPPLEMENTARY | FIRST_MATE, Some("chr1"), 400, 420, 60),
            aln("r2", PAIRED | FIRST_MATE, Some("chr1"), 499, 549, 60),
            aln("r2", PAIRED | SECOND_MATE, Some("chr1"), 899, 949, 60),
        ];
        let mut it = BedpeIter::new(records.clone().into_iter());
        assert_eq!(it.next().unwrap().name, "r1");
        // The supplementary record is now bam1 and r2's first read is bam2, so
        // the names differ, the supplementary is skipped, and the pair emitted
        // is r2's two reads. One real pair survives out of two.
        let second = it.next().unwrap();
        assert_eq!(second.name, "r2");
        assert_eq!(it.orphans(), 1, "the supplementary record counts as an orphan");

        // Filtering first restores the expected pairing.
        let filtered: Vec<Aln> = records.into_iter().filter(|a| is_pairable(a.flags)).collect();
        let mut it = BedpeIter::new(filtered.into_iter());
        assert_eq!(it.next().unwrap().name, "r1");
        assert_eq!(it.next().unwrap().name, "r2");
        assert_eq!(it.orphans(), 0);
    }

    #[test]
    fn is_pairable_rejects_secondary_and_supplementary() {
        assert!(is_pairable(PAIRED | FIRST_MATE));
        assert!(!is_pairable(PAIRED | FIRST_MATE | SECONDARY));
        assert!(!is_pairable(PAIRED | FIRST_MATE | SUPPLEMENTARY));
    }
}
