//! R's `round()` and R's spelling of a double.
//!
//! `gc_correction.R` leans on both:
//!
//! ```text
//! X <- mutate(X, percent.gc = round(percent.gc, 2))
//! ...
//! write.table(X, output.file, sep = "\t", row.names = FALSE, quote = FALSE)
//! ```
//!
//! The first decides which GC bin a domain falls into, and therefore which
//! correction factor multiplies its coverage. The GC file is full of
//! three-decimal values because a 1000 bp domain's GC content is a count over
//! 1000, so exact halves like `0.345` are common, and halves are exactly where
//! a naive round disagrees with R's.
//!
//! The second decides the bytes of the output file. Matching upstream means
//! matching the spelling, not just the value.

/// R's `round(x, digits)`.
///
/// Transcribed from R's `src/nmath/fround.c` as it stands in R 4.x. The
/// algorithm changed in R 4.0.0, so this does not match R 3.x.
///
/// It is not simply "round half to even on the decimal value". R rounds the
/// scaled value both ways, converts each back, and picks whichever lands closer
/// to the original double, falling back to the even rule only on an exact tie.
/// That is why `round(0.345, 2)` and `round(0.001 * 345, 2)` disagree: those are
/// two different doubles, and each is genuinely nearer a different answer.
#[must_use]
pub fn fround(x: f64, digits: f64) -> f64 {
    /// `DBL_MAX_10_EXP`.
    const MAX10E: i32 = 308;
    /// `DBL_MAX_10_EXP + DBL_DIG`.
    const MAX_DIGITS: f64 = 323.0;
    /// `DBL_DIG`.
    const DBL_DIG: f64 = 15.0;
    /// R's `M_LOG10_2`, which rounds to the same double as this constant.
    use std::f64::consts::LOG10_2;

    if x.is_nan() || digits.is_nan() {
        return x + digits;
    }
    if !x.is_finite() {
        return x;
    }
    if digits > MAX_DIGITS || x == 0.0 {
        return x;
    }
    if digits < -f64::from(MAX10E) {
        return 0.0;
    }
    if digits == 0.0 {
        return x.round_ties_even();
    }

    let dig = (digits + 0.5).floor() as i32;

    let mut sgn = 1.0_f64;
    let mut x = x;
    if x < 0.0 {
        sgn = -1.0;
        x = -x;
    }

    // ~= log10(x), the way R computes it.
    let l10x = LOG10_2 * (0.5 + logb(x));
    if l10x + f64::from(dig) > DBL_DIG {
        // Rounding to so many digits that no rounding is needed.
        return sgn * x;
    }

    let (pow10, i10, xd, xu);
    if dig <= MAX10E {
        pow10 = pow_di(10.0, dig);
        let x10 = x * pow10;
        i10 = x10.floor();
        xd = i10 / pow10;
        xu = x10.ceil() / pow10;
    } else {
        // |x| << 1, around 10^-305: scale in two steps so neither overflows.
        let e10 = dig - MAX10E;
        let p10 = pow_di(10.0, e10);
        pow10 = pow_di(10.0, MAX10E);
        let x10 = (x * pow10) * p10;
        i10 = x10.floor();
        xd = i10 / pow10 / p10;
        xu = x10.ceil() / pow10 / p10;
    }

    let du = xu - x;
    let dd = x - xd;
    sgn * if du < dd || (i10 % 2.0 == 1.0 && du == dd) {
        xu
    } else {
        xd
    }
}

/// R's `R_pow_di`: `x^n` by binary exponentiation.
///
/// Kept rather than replaced by `powi` because the multiplication order decides
/// the last bit, and `fround` divides by the result.
#[must_use]
fn pow_di(x: f64, n: i32) -> f64 {
    let mut pow = 1.0_f64;
    let mut x = x;
    let mut n = n;

    if n != 0 {
        if !x.is_finite() {
            return x.powi(n);
        }
        let is_neg = n < 0;
        if is_neg {
            n = -n;
        }
        loop {
            if n & 1 != 0 {
                pow *= x;
            }
            n >>= 1;
            if n != 0 {
                x *= x;
            } else {
                break;
            }
        }
        if is_neg {
            pow = 1.0 / pow;
        }
    }

    pow
}

/// C's `logb`: the unbiased binary exponent, as a double.
///
/// Read off the bit pattern rather than computed as `log2(x).floor()`, which
/// disagrees at exact powers of two.
#[must_use]
fn logb(x: f64) -> f64 {
    let bits = x.abs().to_bits();
    let exp = ((bits >> 52) & 0x7ff) as i32;
    if exp != 0 {
        return f64::from(exp - 1023);
    }
    let mantissa = bits & ((1_u64 << 52) - 1);
    if mantissa == 0 {
        return f64::NEG_INFINITY;
    }
    // Highest set bit, as a power of two, for a subnormal.
    let highest = 63 - mantissa.leading_zeros() as i32;
    f64::from(highest - 1074)
}

/// How R spells a double, as `as.character()` and `write.table()` do.
///
/// This reproduces the *decision* made by R's `formatReal` and `scientific` in
/// `src/main/format.c` rather than transcribing them: take the value to 15
/// significant digits, drop trailing zeros to get the significant-digit count,
/// then measure how wide the fixed and the scientific spelling would be and
/// take the narrower, preferring fixed on a tie. That tie rule is what makes
/// `0.000123` come out fixed (8 characters either way) while `0.0001` comes out
/// as `1e-04` (5 against 6).
///
/// The equivalence is asserted, not assumed: `tests/rmath_parity.rs` checks it
/// against R for every vector in `tests/data/rmath_vectors.tsv`.
///
/// `scipen` is taken as 0, its default. Upstream never changes it.
///
/// # Panics
/// Panics if the standard library stops emitting `{:.14e}` in the form
/// `d.dddddddddddddde[+-]dd`, which the parsing below relies on.
#[must_use]
pub fn format_as_r(x: f64) -> String {
    if x.is_nan() {
        return "NaN".to_string();
    }
    if x.is_infinite() {
        return if x > 0.0 {
            "Inf".to_string()
        } else {
            "-Inf".to_string()
        };
    }
    if x == 0.0 {
        // Covers -0.0, which R spells without the sign.
        return "0".to_string();
    }

    let neg = usize::from(x < 0.0);

    // 15 significant digits, which is what R's `digits` is worth here.
    let sci15 = format!("{:.14e}", x.abs());
    let (mantissa, exponent) = sci15
        .split_once('e')
        .expect("scientific form has an exponent");
    let kpower: i32 = exponent.parse().expect("exponent is an integer");

    // Significant digits actually needed, after dropping trailing zeros.
    let digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
    let nsig = digits.trim_end_matches('0').len().max(1);

    // Width of "1.234e+05" style.
    let exp_digits = exponent
        .trim_start_matches(['+', '-'])
        .trim_start_matches('0')
        .len()
        .max(2);
    let sci_width = neg + nsig + usize::from(nsig > 1) + 2 + exp_digits;

    // Width of "123.4" style.
    let left = usize::try_from(kpower + 1).unwrap_or(1).max(1);
    let rgt = usize::try_from(nsig as i32 - kpower - 1).unwrap_or(0);
    let fixed_width = neg + left + if rgt > 0 { rgt + 1 } else { 0 };

    if fixed_width <= sci_width {
        format!("{x:.rgt$}")
    } else {
        // Take the mantissa from the 15-digit form rather than recomputing it.
        // Dividing by `10^kpower` would be a lossy round trip near the ends of
        // the exponent range: for 9.99999999999999e-301 it lands on 10.0.
        // Truncating is exact here, since `nsig` was defined as the digit count
        // left after stripping trailing zeros from these very digits.
        let kept = &digits[..nsig];
        let m = if nsig > 1 {
            format!("{}.{}", &kept[..1], &kept[1..])
        } else {
            kept.to_string()
        };
        let sign = if x < 0.0 { "-" } else { "" };
        let esign = if kpower < 0 { '-' } else { '+' };
        format!("{sign}{m}e{esign}{:0>2}", kpower.abs())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The pair that shows R's rounding is nearest-value, not nearest-decimal:
    /// two doubles that both print as "0.345" round to different answers.
    #[test]
    fn neighbouring_doubles_round_apart() {
        assert_eq!(fround(0.345, 2.0), 0.34);
        assert_eq!(fround(0.001 * 345.0, 2.0), 0.35);
    }

    #[test]
    fn exact_halves_go_to_even() {
        assert_eq!(fround(0.5, 0.0), 0.0);
        assert_eq!(fround(1.5, 0.0), 2.0);
        assert_eq!(fround(2.5, 0.0), 2.0);
        assert_eq!(fround(-0.5, 0.0), 0.0);
        assert_eq!(fround(-1.5, 0.0), -2.0);
    }

    #[test]
    fn zero_and_non_finite_pass_through() {
        assert_eq!(fround(0.0, 2.0), 0.0);
        assert!(fround(f64::NAN, 2.0).is_nan());
        assert_eq!(fround(f64::INFINITY, 2.0), f64::INFINITY);
    }

    #[test]
    fn fixed_wins_ties_against_scientific() {
        // Eight characters either way, so fixed.
        assert_eq!(format_as_r(0.000123), "0.000123");
        // Six against five, so scientific.
        assert_eq!(format_as_r(0.0001), "1e-04");
        assert_eq!(format_as_r(100000.0), "1e+05");
    }

    #[test]
    fn integers_keep_all_their_digits() {
        assert_eq!(format_as_r(1234567890123456.0), "1234567890123456");
        assert_eq!(format_as_r(1.0), "1");
        assert_eq!(format_as_r(-1.0), "-1");
        assert_eq!(format_as_r(0.0), "0");
        assert_eq!(format_as_r(-0.0), "0");
    }

    #[test]
    // The literal deliberately carries more precision than a double holds:
    // being truncated is the point of the assertion.
    #[allow(clippy::excessive_precision)]
    fn fifteen_significant_digits() {
        assert_eq!(format_as_r(1.0 / 3.0), "0.333333333333333");
        assert_eq!(format_as_r(2.0 / 3.0), "0.666666666666667");
        assert_eq!(format_as_r(123456789.123456789), "123456789.123457");
    }

    /// The three numbers the upstream README quotes as expected output.
    #[test]
    fn readme_values_round_trip() {
        assert_eq!(format_as_r(28.2794), "28.2794");
        assert_eq!(format_as_r(2.55338), "2.55338");
        assert_eq!(format_as_r(0.66548), "0.66548");
    }

    #[test]
    fn three_digit_exponents() {
        assert_eq!(format_as_r(1e300), "1e+300");
        assert_eq!(format_as_r(1e-300), "1e-300");
    }
}
