//! Pins `rmath` to the local R, bit for bit.
//!
//! Vectors come from `tests/oracle/gen_rmath_vectors.R`. Regenerate with:
//!
//! ```text
//! Rscript crates/plethora-compat/tests/oracle/gen_rmath_vectors.R \
//!     > crates/plethora-compat/tests/data/rmath_vectors.tsv
//! ```
//!
//! The rounding corpus covers every three-decimal value from 0 to 1, which is
//! the exact domain of the GC file for the 1000 bp baseline domains, so for that
//! input the check is exhaustive rather than a sample.

use plethora_compat::rmath::{format_as_r, fround};

const VECTORS: &str = include_str!("data/rmath_vectors.tsv");

fn rows(section: &str) -> impl Iterator<Item = Vec<&'static str>> {
    VECTORS
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .filter(move |f| f[0] == section)
}

#[test]
fn round_to_two_digits_matches_r() {
    let mut checked = 0;
    for f in rows("round") {
        let input: f64 = f[1].parse().expect("input is a double");
        let want: f64 = f[2].parse().expect("expected is a double");
        let got = fround(input, 2.0);
        assert_eq!(
            got.to_bits(),
            want.to_bits(),
            "round({input:.17e}, 2): got {got:.17e}, want {want:.17e}"
        );
        checked += 1;
    }
    assert!(
        checked > 1000,
        "expected the full corpus, saw {checked} rows"
    );
}

#[test]
fn formatting_matches_write_table() {
    let mut checked = 0;
    for f in rows("write") {
        let input: f64 = f[1].parse().expect("input is a double");
        let want = f[3];
        assert_eq!(
            format_as_r(input),
            want,
            "format_as_r({input:.17e}) disagrees with write.table"
        );
        checked += 1;
    }
    assert!(
        checked > 20,
        "expected the write corpus, saw {checked} rows"
    );
}
