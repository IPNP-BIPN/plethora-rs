//! Pins the RANDLIB port to Perl's `Math::Random`, bit for bit.
//!
//! The vectors in `tests/data/randlib_vectors.tsv` are produced by
//! `tests/oracle/gen_randlib_vectors.pl`. Regenerate them with:
//!
//! ```text
//! env -u PERL5LIB perl -I.oracle/perl5/lib/perl5 \
//!     crates/plethora-compat/tests/oracle/gen_randlib_vectors.pl \
//!     > crates/plethora-compat/tests/data/randlib_vectors.tsv
//! ```
//!
//! Comparison is on the bit pattern, not a tolerance. `%.17g` round-trips an
//! IEEE double exactly and Rust's float parser is correctly rounded, so any
//! surviving difference is a real difference in the arithmetic.

use plethora_compat::randlib::{Phrtsd, Randlib};

const VECTORS: &str = include_str!("data/randlib_vectors.tsv");

struct Vector {
    phrase: Vec<u8>,
    seed1: i64,
    seed2: i64,
    uniform: Vec<f64>,
    normal: Vec<f64>,
    gennor: f64,
}

fn unhex(s: &str) -> Vec<u8> {
    (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).expect("phrase is valid hex"))
        .collect()
}

fn parse_doubles(field: &str) -> Vec<f64> {
    field
        .split(',')
        .map(|v| v.parse::<f64>().expect("oracle emitted a parseable double"))
        .collect()
}

fn vectors() -> Vec<Vector> {
    VECTORS
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            assert_eq!(f.len(), 6, "malformed vector row: {line}");
            Vector {
                phrase: unhex(f[0]),
                seed1: f[1].parse().expect("seed1 is an integer"),
                seed2: f[2].parse().expect("seed2 is an integer"),
                uniform: parse_doubles(f[3]),
                normal: parse_doubles(f[4]),
                gennor: f[5].parse().expect("gennor is a double"),
            }
        })
        .collect()
}

/// Which `phrtsd` the oracle's `Math::Random` was compiled with.
///
/// `Makefile.PL` only passes `-DPHRTSD_ORIG` when the build was invoked with a
/// matching argument, so the answer depends on how the module was installed
/// rather than on anything in the source. Rather than assume, the variant is
/// resolved from the vectors and every other test uses whatever won.
fn detected_variant() -> Phrtsd {
    let vectors = vectors();
    let matches = |variant: Phrtsd| {
        let rng = Randlib::with_phrtsd(variant);
        vectors
            .iter()
            .all(|v| rng.phrtsd_seeds(&v.phrase) == (v.seed1, v.seed2))
    };

    match (matches(Phrtsd::New), matches(Phrtsd::Orig)) {
        (true, false) => Phrtsd::New,
        (false, true) => Phrtsd::Orig,
        (true, true) => panic!(
            "both phrtsd variants reproduce every vector, so the vectors do not \
             discriminate between them; add a phrase that does"
        ),
        (false, false) => panic!(
            "neither phrtsd variant reproduces the oracle seeds; the port is wrong \
             or Math::Random changed"
        ),
    }
}

#[test]
fn phrtsd_matches_the_oracle_seeds() {
    let variant = detected_variant();
    let rng = Randlib::with_phrtsd(variant);

    for v in vectors() {
        let got = rng.phrtsd_seeds(&v.phrase);
        assert_eq!(
            got,
            (v.seed1, v.seed2),
            "phrtsd({variant:?}) disagrees for phrase {:?}",
            String::from_utf8_lossy(&v.phrase)
        );
    }
}

#[test]
fn ranf_matches_the_oracle_uniforms() {
    let variant = detected_variant();

    for v in vectors() {
        let mut rng = Randlib::with_phrtsd(variant);
        rng.set_seed_from_phrase(&v.phrase);
        for (draw, &want) in v.uniform.iter().enumerate() {
            let got = rng.ranf();
            assert_eq!(
                got.to_bits(),
                want.to_bits(),
                "ranf draw {draw} for phrase {:?}: got {got:.17e}, want {want:.17e}",
                String::from_utf8_lossy(&v.phrase)
            );
        }
    }
}

#[test]
fn snorm_matches_the_oracle_normals() {
    let variant = detected_variant();

    for v in vectors() {
        let mut rng = Randlib::with_phrtsd(variant);
        rng.set_seed_from_phrase(&v.phrase);
        for (draw, &want) in v.normal.iter().enumerate() {
            let got = rng.snorm();
            assert_eq!(
                got.to_bits(),
                want.to_bits(),
                "snorm draw {draw} for phrase {:?}: got {got:.17e}, want {want:.17e}",
                String::from_utf8_lossy(&v.phrase)
            );
        }
    }
}

/// The exact call shape `merge_pairs.pl` makes.
#[test]
fn gennor_matches_the_oracle() {
    let variant = detected_variant();

    for v in vectors() {
        let mut rng = Randlib::with_phrtsd(variant);
        rng.set_seed_from_phrase(&v.phrase);
        let got = rng.gennor(317.0, 45.0);
        assert_eq!(
            got.to_bits(),
            v.gennor.to_bits(),
            "gennor(317, 45) for phrase {:?}: got {got:.17e}, want {:.17e}",
            String::from_utf8_lossy(&v.phrase),
            v.gennor
        );
    }
}

/// A regression guard on the one vector quoted in the module documentation.
#[test]
fn documented_vector_holds() {
    let mut rng = Randlib::with_phrtsd(detected_variant());
    rng.set_seed_from_phrase(b"hello");
    assert_eq!(rng.gennor(0.0, 1.0), -1.117_463_622_136_280_4);
}
