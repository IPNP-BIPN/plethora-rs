//! How `awk` spells a number.
//!
//! The last line of `make_bed.sh` is where the coverage numbers acquire their
//! printed form:
//!
//! ```text
//! awk 'OFS="\t" { print $1, $4 / ($3 - $2 + 1)}' ${output}_coverage.bed > ${output}_read_depth.bed
//! ```
//!
//! `_read_depth.bed` is the file the paper's expected output is quoted from,
//! and those numbers carry exactly six significant digits:
//!
//! ```text
//! NBPF1_CON1_1   28.2794
//! NBPF1L_CON1_1  2.55338
//! NBPF1_CON1_2   0.66548
//! ```
//!
//! That is `OFMT`, awk's output format for numbers, at its default of `%.6g`.
//! Formatting a coverage of 28.279401 as `28.279401` would not be a rounding
//! difference, it would be a different file.
//!
//! One rule sits on top of `OFMT` and is easy to miss: a value that is exactly
//! an integer is printed as an integer, not through `OFMT`. Uncovered domains
//! are common, and every one of them prints `0` rather than `0.00000`.

/// awk's default `OFMT`.
pub const DEFAULT_OFMT_PRECISION: usize = 6;

/// Formats a number the way `awk`'s `print` does.
///
/// Integral values go out as integers; everything else goes through `%.6g`.
#[must_use]
pub fn print_number(value: f64) -> String {
    if value.is_nan() {
        return "nan".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "inf" } else { "-inf" }.to_string();
    }
    // An exactly integral value prints as an integer. The bound keeps this to
    // the range where a double still represents consecutive integers exactly;
    // beyond it awk's own behaviour depends on its integer width, and coverage
    // values never come close.
    if value == value.trunc() && value.abs() < 1e15 {
        return format!("{}", value as i64);
    }
    format_g(value, DEFAULT_OFMT_PRECISION)
}

/// C's `%.*g`, which is what `OFMT` names.
///
/// Chooses between `%e` and `%f` by the decimal exponent, then strips trailing
/// zeros and a bare decimal point. Written out rather than delegated because
/// Rust has no `%g`, and the choice of form is the whole point.
///
/// # Panics
/// Panics if the standard library stops emitting `{:e}` in the form
/// `d.ddde[+-]dd`, which the exponent parsing relies on.
#[must_use]
pub fn format_g(value: f64, precision: usize) -> String {
    // C treats a precision of 0 as 1.
    let p = precision.max(1);

    if value == 0.0 {
        return "0".to_string();
    }

    // The decimal exponent X, as %e would report it.
    let exponent: i32 = {
        let e = format!("{:.*e}", p - 1, value);
        e.split_once('e')
            .map(|(_, x)| x.parse().expect("scientific form has an integer exponent"))
            .expect("scientific form has an exponent")
    };

    let mut s = if exponent < -4 || exponent >= i32::try_from(p).unwrap_or(i32::MAX) {
        // %e with precision p - 1, then a two-digit exponent as C writes it.
        let raw = format!("{:.*e}", p - 1, value);
        let (mantissa, exp) = raw
            .split_once('e')
            .expect("scientific form has an exponent");
        let mantissa = strip_trailing_zeros(mantissa);
        let exp: i32 = exp.parse().expect("integer exponent");
        format!(
            "{mantissa}e{}{:02}",
            if exp < 0 { '-' } else { '+' },
            exp.abs()
        )
    } else {
        // %f with precision p - 1 - X.
        let decimals =
            usize::try_from(i32::try_from(p).unwrap_or(i32::MAX) - 1 - exponent).unwrap_or(0);
        strip_trailing_zeros(&format!("{value:.decimals$}"))
    };

    // "-0" is not a form awk emits for a value that rounded to zero.
    if s == "-0" {
        s = "0".to_string();
    }
    s
}

/// Removes trailing zeros in the fractional part, and the point if it is left
/// bare. Only applies when the text has a decimal point at all.
fn strip_trailing_zeros(s: &str) -> String {
    if !s.contains('.') {
        return s.to_string();
    }
    let trimmed = s.trim_end_matches('0');
    trimmed.strip_suffix('.').unwrap_or(trimmed).to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The three numbers the upstream README quotes as the expected result.
    #[test]
    fn readme_values_round_trip() {
        assert_eq!(print_number(28.279_401_23), "28.2794");
        assert_eq!(print_number(2.553_380_456), "2.55338");
        assert_eq!(print_number(0.665_48), "0.66548");
    }

    /// The rule that sits on top of OFMT: integers print as integers. Every
    /// uncovered domain goes through this branch.
    #[test]
    fn integral_values_print_as_integers() {
        assert_eq!(print_number(0.0), "0");
        assert_eq!(print_number(-0.0), "0");
        assert_eq!(print_number(1.0), "1");
        assert_eq!(print_number(2.0), "2");
        assert_eq!(print_number(25.0), "25");
        assert_eq!(print_number(1_000_000.0), "1000000");
    }

    #[test]
    fn six_significant_digits() {
        assert_eq!(print_number(1.0 / 3.0), "0.333333");
        assert_eq!(print_number(60.0 / 101.0), "0.594059");
        assert_eq!(print_number(28564.0 / 1001.0), "28.5355");
        assert_eq!(print_number(1.0 / 101.0), "0.00990099");
    }

    #[test]
    fn short_values_keep_no_trailing_zeros() {
        assert_eq!(print_number(0.5), "0.5");
        assert_eq!(print_number(0.125), "0.125");
        assert_eq!(print_number(1.5), "1.5");
    }

    /// Beyond the %f window, %g switches to exponential with a two-digit
    /// exponent.
    #[test]
    fn very_small_and_very_large_go_exponential() {
        assert_eq!(format_g(0.000_012_345_6, 6), "1.23456e-05");
        assert_eq!(format_g(1_234_567.0, 6), "1.23457e+06");
        assert_eq!(format_g(1e-300, 6), "1e-300");
    }

    #[test]
    fn negative_values_keep_their_sign() {
        assert_eq!(print_number(-0.5), "-0.5");
        assert_eq!(print_number(-1.0), "-1");
        assert_eq!(print_number(-1.0 / 3.0), "-0.333333");
    }
}
