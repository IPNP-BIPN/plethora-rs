//! `samtools sort -n` ordering.
//!
//! `make_bed.sh` name-sorts the alignment before converting it to BEDPE:
//!
//! ```text
//! samtools sort -n -@ 12 -m 2G -o ${output}_sorted.bam $bam
//! bedtools bamtobed -split -bedpe -i ${output}_sorted.bam > $output.bed
//! ```
//!
//! Two things ride on this order. Mates must come out adjacent, or `bamtobed
//! -bedpe` cannot pair them; that only needs *some* grouping by name. But
//! `merge_pairs.pl` also stops accumulating fragment lengths after 50 million
//! records, and a whole-genome sample has more pairs than that, so the mean and
//! standard deviation it derives, and therefore every extended single-end read
//! in the output, depend on *which* 50 million came first. Reproducing the
//! order exactly is the difference between matching upstream and merely
//! resembling it.
//!
//! Transcribed from samtools 1.24, `bam_sort.c`: [`strnum_cmp`] from the
//! function of the same name at line 172, and [`pair_rank`] from the `QueryName`
//! arm of `heap_lt` at line 242.

use std::cmp::Ordering;

/// samtools' own digit test, which is ASCII-only by construction.
#[inline]
const fn is_digit(c: u8) -> bool {
    c <= b'9' && c >= b'0'
}

/// Reads a C string one past its end as the NUL terminator.
///
/// The C code walks pointers and relies on the terminator to stop; QNAMEs
/// arrive here as slices without one, so the terminator is synthesised. BAM
/// forbids interior NULs in QNAME, so this is faithful.
#[inline]
fn at(s: &[u8], i: usize) -> u8 {
    if i < s.len() { s[i] } else { 0 }
}

/// Natural alphanumeric comparison: `a7b` sorts before `a12b`.
///
/// This is samtools' `strnum_cmp` with `natural_sort` left at its default of 1.
/// Passing `-N` to `samtools sort` would switch it to plain `strcmp`, which
/// plethora never does.
///
/// Leading zeros are skipped, so `a1`, `a01` and `a001` all compare equal;
/// runs of digits are ranked by length first and by the earliest differing
/// digit second, which is what lets it order numbers wider than any integer
/// type.
#[must_use]
pub fn strnum_cmp(a: &[u8], b: &[u8]) -> Ordering {
    let mut pa = 0_usize;
    let mut pb = 0_usize;

    while at(a, pa) != 0 && at(b, pb) != 0 {
        if !is_digit(at(a, pa)) || !is_digit(at(b, pb)) {
            if at(a, pa) != at(b, pb) {
                return at(a, pa).cmp(&at(b, pb));
            }
            pa += 1;
            pb += 1;
        } else {
            // Skip leading zeros.
            while at(a, pa) == b'0' {
                pa += 1;
            }
            while at(b, pb) == b'0' {
                pb += 1;
            }

            // Skip matching digits.
            while is_digit(at(a, pa)) && at(a, pa) == at(b, pb) {
                pa += 1;
                pb += 1;
            }

            // Now mismatching, so see which ends the number sooner.
            let diff = i32::from(at(a, pa)) - i32::from(at(b, pb));
            while is_digit(at(a, pa)) && is_digit(at(b, pb)) {
                pa += 1;
                pb += 1;
            }

            if is_digit(at(a, pa)) {
                return Ordering::Greater; // pa still going, so larger
            } else if is_digit(at(b, pb)) {
                return Ordering::Less; // pb still going, so larger
            } else if diff != 0 {
                return diff.cmp(&0); // same length, so earlier diff
            }
        }
    }

    if at(a, pa) != 0 {
        Ordering::Greater
    } else if at(b, pb) != 0 {
        Ordering::Less
    } else {
        Ordering::Equal
    }
}

/// Which samtools the tie-break should follow.
///
/// The rule changed in 1.20. Up to 1.19 a tie on QNAME was broken by the pair
/// bits alone, so a secondary and a supplementary record sharing those bits
/// kept their input order; from 1.20 the secondary and supplementary bits enter
/// the key as well.
///
/// This is unreachable from the plethora pipeline either way: bowtie2 emits
/// neither kind of record, and the BWA-MEM path filters both before sorting.
/// It is carried because the differential test runs against whatever samtools
/// is installed, and a runner with an older one would otherwise report a
/// disagreement that is a version difference rather than an error.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TieBreak {
    /// samtools 1.20 and later: READ1, READ2, then primary, supplementary,
    /// secondary.
    #[default]
    Since1_20,
    /// samtools 1.19 and earlier: the pair bits only.
    Before1_20,
}

impl TieBreak {
    /// Reads a `samtools --version` line and picks the matching rule.
    ///
    /// Anything unparseable is taken as current, since that is what a fresh
    /// install gives.
    #[must_use]
    pub fn from_version(version: &str) -> Self {
        let Some(rest) = version.split_whitespace().nth(1) else {
            return Self::Since1_20;
        };
        let mut parts = rest.split('.');
        let (Some(Ok(major)), Some(Ok(minor))) = (
            parts.next().map(str::parse::<u32>),
            parts.next().map(str::parse::<u32>),
        ) else {
            return Self::Since1_20;
        };
        if (major, minor) < (1, 20) {
            Self::Before1_20
        } else {
            Self::Since1_20
        }
    }
}

/// The secondary key: where a record sits among others sharing its QNAME.
///
/// samtools packs three flag bits into one integer so the comparison stays a
/// plain subtraction, giving the order READ1, READ2, then primary,
/// supplementary, secondary.
#[must_use]
pub const fn pair_rank(flag: u16) -> u32 {
    pair_rank_with(flag, TieBreak::Since1_20)
}

/// [`pair_rank`] under a named samtools rule.
#[must_use]
pub const fn pair_rank_with(flag: u16, rule: TieBreak) -> u32 {
    let f = flag as u32;
    match rule {
        TieBreak::Since1_20 => ((f & 0xc0) << 8) | ((f & 0x100) << 3) | ((f & 0x800) >> 3),
        TieBreak::Before1_20 => f & 0xc0,
    }
}

/// The complete `samtools sort -n` record comparison.
///
/// Records tying on both keys keep their input order: samtools' merge is
/// stable, so callers must pair this with a stable sort.
#[must_use]
pub fn cmp_by_qname(a_qname: &[u8], a_flag: u16, b_qname: &[u8], b_flag: u16) -> Ordering {
    cmp_by_qname_with(a_qname, a_flag, b_qname, b_flag, TieBreak::Since1_20)
}

/// [`cmp_by_qname`] under a named samtools rule.
#[must_use]
pub fn cmp_by_qname_with(
    a_qname: &[u8],
    a_flag: u16,
    b_qname: &[u8],
    b_flag: u16,
    rule: TieBreak,
) -> Ordering {
    strnum_cmp(a_qname, b_qname)
        .then_with(|| pair_rank_with(a_flag, rule).cmp(&pair_rank_with(b_flag, rule)))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Leading zeros are invisible to the comparison.
    #[test]
    fn leading_zeros_compare_equal() {
        assert_eq!(strnum_cmp(b"a1", b"a01"), Ordering::Equal);
        assert_eq!(strnum_cmp(b"a1", b"a001"), Ordering::Equal);
        assert_eq!(strnum_cmp(b"a0", b"a00"), Ordering::Equal);
    }

    /// Digit runs rank by length before value, which is what makes numbers
    /// wider than u64 sort correctly.
    #[test]
    fn longer_digit_runs_are_larger() {
        assert_eq!(strnum_cmp(b"a2", b"a10"), Ordering::Less);
        assert_eq!(
            strnum_cmp(b"9999999999999999999", b"10000000000000000000"),
            Ordering::Less
        );
    }

    #[test]
    fn digits_outrank_the_end_of_a_name() {
        assert_eq!(strnum_cmp(b"a0", b"a1"), Ordering::Less);
        assert_eq!(strnum_cmp(b"read", b"read0"), Ordering::Less);
        assert_eq!(strnum_cmp(b"a2", b"a2b1"), Ordering::Less);
    }

    /// Outside digit runs the comparison is plain ASCII.
    #[test]
    fn non_digits_compare_as_bytes() {
        assert_eq!(strnum_cmp(b"A", b"a"), Ordering::Less);
        assert_eq!(strnum_cmp(b"b1", b"b1c"), Ordering::Less);
        assert_eq!(strnum_cmp(b"", b""), Ordering::Equal);
        assert_eq!(strnum_cmp(b"", b"a"), Ordering::Less);
    }

    /// Embedded numbers are compared run by run.
    #[test]
    fn multiple_digit_runs() {
        assert_eq!(strnum_cmp(b"a1b2", b"a1b10"), Ordering::Less);
        assert_eq!(strnum_cmp(b"x2y", b"x10y"), Ordering::Less);
    }

    /// READ1 before READ2, primary before supplementary before secondary.
    #[test]
    fn pair_rank_orders_the_flag_bits() {
        let read1 = pair_rank(0x40);
        let read2 = pair_rank(0x80);
        let both = pair_rank(0xc0);
        assert!(pair_rank(0) < read1 && read1 < read2 && read2 < both);

        let primary = pair_rank(0x40);
        let supplementary = pair_rank(0x40 | 0x800);
        let secondary = pair_rank(0x40 | 0x100);
        assert!(primary < supplementary && supplementary < secondary);
    }

    /// The rule changed in samtools 1.20, and the parser has to place a version
    /// on the right side of that.
    #[test]
    fn the_tie_break_follows_the_samtools_version() {
        assert_eq!(TieBreak::from_version("samtools 1.24"), TieBreak::Since1_20);
        assert_eq!(TieBreak::from_version("samtools 1.20"), TieBreak::Since1_20);
        assert_eq!(
            TieBreak::from_version("samtools 1.19"),
            TieBreak::Before1_20
        );
        assert_eq!(
            TieBreak::from_version("samtools 1.13"),
            TieBreak::Before1_20
        );
        assert_eq!(TieBreak::from_version("samtools 2.0"), TieBreak::Since1_20);
        // Unparseable reads as current, which is what a fresh install gives.
        assert_eq!(TieBreak::from_version("nonsense"), TieBreak::Since1_20);
    }

    /// Under the older rule a secondary and a supplementary record sharing
    /// their pair bits compare equal, so the sort leaves them in input order.
    #[test]
    fn the_older_rule_ignores_the_secondary_and_supplementary_bits() {
        let supplementary = 0x80 | 0x800;
        let secondary = 0x80 | 0x100;
        assert_ne!(
            pair_rank_with(supplementary, TieBreak::Since1_20),
            pair_rank_with(secondary, TieBreak::Since1_20)
        );
        assert_eq!(
            pair_rank_with(supplementary, TieBreak::Before1_20),
            pair_rank_with(secondary, TieBreak::Before1_20)
        );
    }

    #[test]
    fn qname_beats_the_flag() {
        assert_eq!(cmp_by_qname(b"a1", 0x80, b"a2", 0x40), Ordering::Less);
        assert_eq!(cmp_by_qname(b"a1", 0x40, b"a1", 0x80), Ordering::Less);
        assert_eq!(cmp_by_qname(b"a1", 0x40, b"a01", 0x40), Ordering::Equal);
    }
}
