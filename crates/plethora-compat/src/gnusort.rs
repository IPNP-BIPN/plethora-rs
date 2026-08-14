//! GNU `sort -k1,1 -k2,2n` ordering.
//!
//! `make_bed.sh` sorts the read intervals before intersecting them:
//!
//! ```text
//! sort -k 1,1 -k 2,2n -T ./ ${output}_edited.bed > ${output}_sorted.bed
//! ```
//!
//! The two named keys are the obvious part. The part that decides actual output
//! is the third, unnamed key: GNU `sort` is not stable unless asked, and it
//! breaks a tie on every named key by comparing the whole line byte for byte.
//! Lines sharing a chromosome and a start position are common in a BED file of
//! read intervals, so this last-resort comparison orders a large fraction of
//! the file. Reproducing only the two named keys would produce a differently
//! ordered file that still looks sorted.
//!
//! Scope: this emulates the exact invocation above, not GNU `sort` in general.
//! That matters for one detail. With no `-t`, a key for field 2 formally
//! includes the blanks in front of it, because the `b` modifier was not given.
//! Here that is provably inert: field 1 has no blanks before it, and field 2 is
//! compared numerically, which skips leading blanks anyway. A different `-k`
//! specification would need the general rule.
//!
//! Locale: comparison is byte-wise, as in `LC_ALL=C`. Under a UTF-8 locale GNU
//! `sort` compares through `strcoll`, which folds case and punctuation
//! differently and would order the file differently. Upstream never pins the
//! locale, so its own output is locale-dependent; see `DIVERGENCES.md`.

use std::cmp::Ordering;

/// GNU `sort`'s notion of a blank in the C locale.
#[inline]
const fn is_blank(c: u8) -> bool {
    c == b' ' || c == b'\t'
}

/// The `n`-th blank-separated field, 1-based, or an empty slice if absent.
///
/// Fields are maximal runs of non-blanks. Runs of blanks separate them, so a
/// line with repeated tabs has empty fields collapsed, exactly as GNU `sort`
/// treats them by default.
#[must_use]
pub fn field(line: &[u8], n: usize) -> &[u8] {
    let mut start = 0;
    let mut remaining = n;

    while start < line.len() {
        while start < line.len() && is_blank(line[start]) {
            start += 1;
        }
        let mut end = start;
        while end < line.len() && !is_blank(line[end]) {
            end += 1;
        }
        remaining -= 1;
        if remaining == 0 {
            return &line[start..end];
        }
        start = end;
    }

    &[]
}

/// A decimal number as GNU `sort -n` reads it.
struct Number<'a> {
    negative: bool,
    /// Integer digits with leading zeros removed; empty means zero.
    integer: &'a [u8],
    /// Fraction digits with trailing zeros left in place.
    fraction: &'a [u8],
}

impl<'a> Number<'a> {
    /// Parses as much of a field as looks numeric, treating the rest as absent.
    ///
    /// Anything unparseable reads as zero, which is what GNU `sort -n` does
    /// with a non-numeric field rather than erroring.
    fn parse(s: &'a [u8]) -> Self {
        let mut i = 0;
        while i < s.len() && is_blank(s[i]) {
            i += 1;
        }

        let mut negative = false;
        if i < s.len() && (s[i] == b'-' || s[i] == b'+') {
            negative = s[i] == b'-';
            i += 1;
        }

        let int_start = i;
        while i < s.len() && s[i].is_ascii_digit() {
            i += 1;
        }
        let mut integer = &s[int_start..i];
        while !integer.is_empty() && integer[0] == b'0' {
            integer = &integer[1..];
        }

        let mut fraction: &[u8] = &[];
        if i < s.len() && s[i] == b'.' {
            i += 1;
            let frac_start = i;
            while i < s.len() && s[i].is_ascii_digit() {
                i += 1;
            }
            fraction = &s[frac_start..i];
        }

        Self {
            negative,
            integer,
            fraction,
        }
    }

    /// True when every digit is zero, so that `-0` and `0` compare equal.
    fn is_zero(&self) -> bool {
        self.integer.is_empty() && self.fraction.iter().all(|&d| d == b'0')
    }

    /// Compares magnitudes only.
    fn cmp_magnitude(&self, other: &Self) -> Ordering {
        // More integer digits means a larger number, since leading zeros are gone.
        self.integer
            .len()
            .cmp(&other.integer.len())
            .then_with(|| self.integer.cmp(other.integer))
            .then_with(|| {
                // Compare fractions digit by digit, treating a missing digit as zero.
                let n = self.fraction.len().max(other.fraction.len());
                for i in 0..n {
                    let a = self.fraction.get(i).copied().unwrap_or(b'0');
                    let b = other.fraction.get(i).copied().unwrap_or(b'0');
                    match a.cmp(&b) {
                        Ordering::Equal => {}
                        other_order => return other_order,
                    }
                }
                Ordering::Equal
            })
    }
}

/// Numeric comparison of two fields, as `sort -n` performs it.
///
/// Comparison is on the decimal text rather than through `f64`, so it stays
/// exact for values wider than a double can hold. BED coordinates never get
/// that large, but the sum column downstream can.
#[must_use]
pub fn numeric_cmp(a: &[u8], b: &[u8]) -> Ordering {
    let x = Number::parse(a);
    let y = Number::parse(b);

    match (x.is_zero(), y.is_zero()) {
        (true, true) => return Ordering::Equal,
        (true, false) => {
            return if y.negative {
                Ordering::Greater
            } else {
                Ordering::Less
            };
        }
        (false, true) => {
            return if x.negative {
                Ordering::Less
            } else {
                Ordering::Greater
            };
        }
        (false, false) => {}
    }

    match (x.negative, y.negative) {
        (false, true) => Ordering::Greater,
        (true, false) => Ordering::Less,
        (true, true) => y.cmp_magnitude(&x),
        (false, false) => x.cmp_magnitude(&y),
    }
}

/// `sort -k1,1 -k2,2n`, including the whole-line last resort.
///
/// Lines are compared without their trailing newline; pass them already
/// stripped.
#[must_use]
pub fn cmp_k1_k2n(a: &[u8], b: &[u8]) -> Ordering {
    field(a, 1)
        .cmp(field(b, 1))
        .then_with(|| numeric_cmp(field(a, 2), field(b, 2)))
        // Last resort: GNU sort is not stable unless given -s, and falls back
        // to comparing entire lines.
        .then_with(|| a.cmp(b))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fields_are_runs_of_non_blanks() {
        assert_eq!(field(b"chr1\t100\tname", 1), b"chr1");
        assert_eq!(field(b"chr1\t100\tname", 2), b"100");
        assert_eq!(field(b"chr1\t100\tname", 3), b"name");
        assert_eq!(field(b"chr1\t100\tname", 4), b"");
        assert_eq!(field(b"", 1), b"");
    }

    #[test]
    fn numeric_key_beats_byte_order() {
        // "10" sorts after "2" numerically but before it byte-wise.
        assert_eq!(cmp_k1_k2n(b"a\t2\tx", b"a\t10\tx"), Ordering::Less);
    }

    /// The behaviour that pins the last-resort rule: numerically equal keys,
    /// ordered by the whole line. Confirmed against `LC_ALL=C gsort`.
    #[test]
    fn equal_keys_fall_back_to_the_whole_line() {
        assert_eq!(cmp_k1_k2n(b"a\t02\tsame", b"a\t2\tbbb"), Ordering::Less);
        assert_eq!(cmp_k1_k2n(b"a\t2\tbbb", b"a\t2\tqqq"), Ordering::Less);
        assert_eq!(cmp_k1_k2n(b"b\t10\taaa", b"b\t10\tzzz"), Ordering::Less);
    }

    #[test]
    fn first_key_is_a_byte_comparison() {
        assert_eq!(cmp_k1_k2n(b"A\t5\tx", b"a\t2\tx"), Ordering::Less);
        assert_eq!(cmp_k1_k2n(b"chr10\t1\tx", b"chr9\t1\tx"), Ordering::Less);
    }

    #[test]
    fn numeric_handles_signs_and_fractions() {
        assert_eq!(numeric_cmp(b"-1", b"1"), Ordering::Less);
        assert_eq!(numeric_cmp(b"-2", b"-1"), Ordering::Less);
        assert_eq!(numeric_cmp(b"0", b"-0"), Ordering::Equal);
        assert_eq!(numeric_cmp(b"1.5", b"1.25"), Ordering::Greater);
        assert_eq!(numeric_cmp(b"1.5", b"1.50"), Ordering::Equal);
        assert_eq!(numeric_cmp(b"007", b"7"), Ordering::Equal);
        // Wider than f64 can represent exactly.
        assert_eq!(
            numeric_cmp(b"9007199254740993", b"9007199254740992"),
            Ordering::Greater
        );
    }

    #[test]
    fn non_numeric_fields_read_as_zero() {
        assert_eq!(numeric_cmp(b"abc", b"0"), Ordering::Equal);
        assert_eq!(numeric_cmp(b"", b"0"), Ordering::Equal);
        assert_eq!(numeric_cmp(b"abc", b"1"), Ordering::Less);
    }
}
