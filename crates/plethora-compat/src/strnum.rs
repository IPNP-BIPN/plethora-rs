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

/// The secondary key: where a record sits among others sharing its QNAME.
///
/// samtools packs three flag bits into one integer so the comparison stays a
/// plain subtraction, giving the order READ1, READ2, then primary,
/// supplementary, secondary. Bowtie2 emits neither supplementary nor secondary
/// records, so under the paper's pipeline only the 0xc0 bits matter; BWA-MEM
/// emits both, which is exactly why the full key is carried here.
#[must_use]
pub const fn pair_rank(flag: u16) -> u32 {
    let f = flag as u32;
    ((f & 0xc0) << 8) | ((f & 0x100) << 3) | ((f & 0x800) >> 3)
}

/// The complete `samtools sort -n` record comparison.
///
/// Records tying on both keys keep their input order: samtools' merge is
/// stable, so callers must pair this with a stable sort.
#[must_use]
pub fn cmp_by_qname(a_qname: &[u8], a_flag: u16, b_qname: &[u8], b_flag: u16) -> Ordering {
    strnum_cmp(a_qname, b_qname).then_with(|| pair_rank(a_flag).cmp(&pair_rank(b_flag)))
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

    #[test]
    fn qname_beats_the_flag() {
        assert_eq!(cmp_by_qname(b"a1", 0x80, b"a2", 0x40), Ordering::Less);
        assert_eq!(cmp_by_qname(b"a1", 0x40, b"a1", 0x80), Ordering::Less);
        assert_eq!(cmp_by_qname(b"a1", 0x40, b"a01", 0x40), Ordering::Equal);
    }
}
