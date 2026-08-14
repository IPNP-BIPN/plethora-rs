//! Pins the loess port to the local R, stage by stage.
//!
//! Vectors come from `tests/oracle/gen_loess_vectors.R`, which exposes R's own
//! internal k-d tree through `fit$kd`. Regenerate with:
//!
//! ```text
//! Rscript crates/plethora-compat/tests/oracle/gen_loess_vectors.R \
//!     > crates/plethora-compat/tests/data/loess_vectors.tsv
//! ```
//!
//! Checking the tree separately from the fitted values is the point. An
//! end-to-end mismatch says only that something is wrong; a tree that matches
//! and vertex values that do not says the geometry is right and the regression
//! is not.

use plethora_compat::loess::Sizing;
use plethora_compat::loess::tree::KdTree;

const VECTORS: &str = include_str!("data/loess_vectors.tsv");

struct Case {
    name: &'static str,
    nc: usize,
    nv: usize,
    vert: Vec<f64>,
    a: Vec<usize>,
    xi: Vec<f64>,
    x: Vec<f64>,
}

fn doubles(field: &str) -> Vec<f64> {
    field
        .split(',')
        .map(|v| v.parse::<f64>().expect("oracle emitted a parseable double"))
        .collect()
}

fn cases() -> Vec<Case> {
    VECTORS
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert_eq!(f.len(), 11, "malformed vector row for {}", f[0]);
            Case {
                name: f[0],
                nc: f[2].parse().expect("nc is an integer"),
                nv: f[3].parse().expect("nv is an integer"),
                vert: doubles(f[4]),
                a: f[5]
                    .split(',')
                    .map(|v| v.parse::<usize>().expect("a is an integer"))
                    .collect(),
                xi: doubles(f[6]),
                x: doubles(f[8]),
            }
        })
        .collect()
}

/// R's own sizing: `nvmax = max(200, n)`, and `ncmax` follows it.
fn build(case: &Case) -> KdTree {
    let n = case.x.len();
    let sizing = Sizing::new(n, 0.75, 0.2, 2);
    let nvmax = n.max(200);
    KdTree::build(&case.x, sizing.fc, 0.0, nvmax, nvmax)
}

#[test]
fn bounding_box_matches_r() {
    for case in cases() {
        let tree = build(&case);
        assert_eq!(
            tree.v[1].to_bits(),
            case.vert[0].to_bits(),
            "{}: lower vertex, got {}, want {}",
            case.name,
            tree.v[1],
            case.vert[0]
        );
        assert_eq!(
            tree.v[tree.vc].to_bits(),
            case.vert[1].to_bits(),
            "{}: upper vertex, got {}, want {}",
            case.name,
            tree.v[tree.vc],
            case.vert[1]
        );
    }
}

#[test]
fn cell_and_vertex_counts_match_r() {
    for case in cases() {
        let tree = build(&case);
        assert_eq!(tree.nc, case.nc, "{}: cell count", case.name);
        assert_eq!(tree.nv, case.nv, "{}: vertex count", case.name);
    }
}

/// `a` is the split dimension per cell, 0 for a leaf. Matching it means the
/// tree has the same shape and the same cell numbering, which is what decides
/// the interpolation cells.
#[test]
fn split_dimensions_match_r() {
    for case in cases() {
        let tree = build(&case);
        let got: Vec<usize> = (1..=tree.nc).map(|p| tree.a[p]).collect();
        assert_eq!(got, case.a, "{}: split dimensions", case.name);
    }
}

/// `xi` is the split coordinate. Only internal cells have one; R leaves the
/// leaf slots at whatever they were, so only the internal ones are compared.
#[test]
fn split_coordinates_match_r() {
    for case in cases() {
        let tree = build(&case);
        for p in 1..=tree.nc {
            if tree.a[p] != 0 {
                assert_eq!(
                    tree.xi[p].to_bits(),
                    case.xi[p - 1].to_bits(),
                    "{}: split coordinate for cell {p}, got {}, want {}",
                    case.name,
                    tree.xi[p],
                    case.xi[p - 1]
                );
            }
        }
    }
}

/// Every data point must land in a leaf whose vertex interval brackets it, or
/// the Hermite interpolation would extrapolate.
#[test]
fn every_point_lands_in_a_bracketing_leaf() {
    for case in cases() {
        let tree = build(&case);
        for &z in &case.x {
            let cell = tree.locate(z);
            assert_eq!(
                tree.a[cell], 0,
                "{}: locate returned an internal cell",
                case.name
            );
            let lo = tree.v[tree.c[cell][1]];
            let hi = tree.v[tree.c[cell][tree.vc]];
            assert!(
                lo <= z && z <= hi,
                "{}: {z} outside [{lo}, {hi}]",
                case.name
            );
        }
    }
}

// ---------------------------------------------------------------------------
// The fitting stages, once the geometry is known to be right.
// ---------------------------------------------------------------------------

/// Relative difference against a scale, since several vertex values and
/// fitted values sit near zero after cancellation. Counting ULP there is
/// misleading: a relative gap of 1e-13 spans thousands of representable steps.
fn relative_error(got: f64, want: f64, scale: f64) -> f64 {
    if got == want {
        return 0.0;
    }
    (got - want).abs() / scale.max(1e-300)
}

/// A case together with R's `vval`, response and fitted values.
type FullCase = (Case, Vec<f64>, Vec<f64>, Vec<f64>);

fn full_case(line: &'static str) -> FullCase {
    let f: Vec<&str> = line.split('\t').collect();
    let case = Case {
        name: f[0],
        nc: f[2].parse().expect("nc"),
        nv: f[3].parse().expect("nv"),
        vert: doubles(f[4]),
        a: f[5].split(',').map(|v| v.parse().expect("a")).collect(),
        xi: doubles(f[6]),
        x: doubles(f[8]),
    };
    (case, doubles(f[7]), doubles(f[9]), doubles(f[10]))
}

fn full_cases() -> Vec<FullCase> {
    VECTORS
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(full_case)
        .collect()
}

/// `fit$kd$vval` interleaves value and derivative per vertex.
///
/// Checked before the interpolation, so that a failure here localises to the
/// regression rather than to the Hermite basis.
#[test]
fn vertex_fits_match_r() {
    let mut worst = 0.0_f64;
    let mut worst_where = String::new();

    for (case, vval, y, _fitted) in full_cases() {
        let fit = plethora_compat::loess::eval::Loess::fit(&case.x, &y);
        assert_eq!(fit.vval.len(), case.nv, "{}: vertex count", case.name);

        // Scale values against the spread of the response, derivatives against
        // the spread of the derivatives; the two have different units.
        let y_scale = y.iter().copied().fold(f64::NEG_INFINITY, f64::max)
            - y.iter().copied().fold(f64::INFINITY, f64::min);
        let d_scale = (0..case.nv)
            .map(|l| vval[2 * l + 1].abs())
            .fold(0.0_f64, f64::max);

        for l in 1..=case.nv {
            for (offset, got, scale) in [
                (0, fit.vval[l - 1].value, y_scale),
                (1, fit.vval[l - 1].derivative, d_scale),
            ] {
                let want = vval[2 * (l - 1) + offset];
                let rel = relative_error(got, want, scale);
                if rel > worst {
                    worst = rel;
                    worst_where = format!(
                        "{} vertex {l} {}: got {got:.17e}, want {want:.17e}",
                        case.name,
                        if offset == 0 { "value" } else { "derivative" }
                    );
                }
            }
        }
    }

    println!("worst vertex-fit relative error: {worst:.3e} ({worst_where})");
    assert!(
        worst < 1e-9,
        "vertex fits drifted {worst:.3e} relative from R, far beyond a BLAS \
         summation-order difference: {worst_where}"
    );
}

/// The end of the chain: what `gc_correction.R` divides by.
///
/// Reported as a relative error against the spread of the fitted values, which
/// is the meaningful scale. Counting ULP is misleading here: several fitted
/// values sit near zero after cancellation, where a relative difference of
/// 1e-13 spans thousands of representable steps.
#[test]
fn fitted_values_match_r() {
    let mut worst_rel = 0.0_f64;
    let mut worst_where = String::new();
    let mut exact = 0_usize;
    let mut total = 0_usize;

    for (case, _vval, y, fitted) in full_cases() {
        let fit = plethora_compat::loess::eval::Loess::fit(&case.x, &y);
        let got = fit.fitted(&case.x);
        assert_eq!(got.len(), fitted.len(), "{}: length", case.name);

        let lo = fitted.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = fitted.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let scale = (hi - lo).max(hi.abs()).max(1e-300);

        let mut case_worst = 0.0_f64;
        for (i, (&g, &w)) in got.iter().zip(&fitted).enumerate() {
            total += 1;
            if g == w {
                exact += 1;
            }
            let rel = (g - w).abs() / scale;
            if rel > case_worst {
                case_worst = rel;
            }
            if rel > worst_rel {
                worst_rel = rel;
                worst_where = format!("{} point {i}: got {g:.17e}, want {w:.17e}", case.name);
            }
        }
        println!("  {:<16} worst relative error {case_worst:.3e}", case.name);
    }

    println!(
        "fitted values: {exact}/{total} bit-exact, worst relative {worst_rel:.3e} ({worst_where})"
    );
    assert!(
        worst_rel < 1e-9,
        "fitted values drifted {worst_rel:.3e} relative from R, far beyond a \
         BLAS summation-order difference: {worst_where}"
    );
}
