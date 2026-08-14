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
            assert_eq!(tree.a[cell], 0, "{}: locate returned an internal cell", case.name);
            let lo = tree.v[tree.c[cell][1]];
            let hi = tree.v[tree.c[cell][tree.vc]];
            assert!(lo <= z && z <= hi, "{}: {z} outside [{lo}, {hi}]", case.name);
        }
    }
}
