//! R's `loess`, reproduced closely enough to match `predict()` bit for bit.
//!
//! `gc_correction.R` fits one:
//!
//! ```text
//! loess_fit <- function(x, y) { fit <- loess(y ~ x); predict(fit) }
//! ```
//!
//! and divides by the result to get the per-GC-bin correction factor `k.gc`.
//! Everything downstream, including the copy numbers the paper reports, is a
//! function of those fitted values.
//!
//! The defaults it inherits are `span = 0.75`, `degree = 2`,
//! `family = "gaussian"`, `surface = "interpolate"` and `cell = 0.2`, with one
//! predictor and one iteration. Two of those matter more than they look:
//!
//! - `surface = "interpolate"` means `predict()` reads off a cubic interpolant
//!   built over a k-d tree, not off the local regression. See [`tree`].
//! - `family = "gaussian"` means a single pass, with no robustness reweighting,
//!   which removes the iteration loop entirely.
//!
//! `normalize = TRUE` is also a default, but R only applies it for more than one
//! predictor, so it is inert here.

pub mod blas;
pub mod eval;
pub mod fit;
pub mod linpack;
pub mod tree;

/// The parameters R derives from `loess.control` before any fitting happens.
///
/// Split out because they are pure arithmetic on `n`, `span` and `cell`, and
/// getting one of them wrong shifts the whole tree without any obvious symptom.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Sizing {
    /// Points in the local neighbourhood, `min(n, floor(n * span))`.
    pub nf: usize,
    /// Bucket size: a cell holding this many points or fewer is a leaf.
    pub fc: usize,
    /// Coefficients in the local polynomial, `(d + 2)(d + 1) / 2` at degree 2.
    pub k: usize,
}

impl Sizing {
    /// Derives the sizes for `n` points at the given span and cell.
    ///
    /// # Panics
    /// Panics unless `degree` is 0, 1 or 2, which is all `loess` accepts.
    #[must_use]
    pub fn new(n: usize, span: f64, cell: f64, degree: usize) -> Self {
        assert!(degree <= 2, "loess accepts degree 0, 1 or 2, got {degree}");
        let d = 1_usize;

        let nf = n.min((n as f64 * span).floor() as usize);
        let fc = (n as f64 * cell * span).floor() as usize;
        let k = match degree {
            0 => 1,
            1 => d + 1,
            _ => (d + 2) * (d + 1) / 2,
        };

        Self { nf, fc, k }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The sizes R reports for the 53-point input used throughout these tests:
    /// `fit$kd$parameter` gives nc = 15 and nv = 9, which is what a bucket size
    /// of 7 produces.
    #[test]
    fn sizing_matches_r_for_the_reference_input() {
        let s = Sizing::new(53, 0.75, 0.2, 2);
        assert_eq!(s.nf, 39, "floor(53 * 0.75)");
        assert_eq!(s.fc, 7, "floor(53 * 0.2 * 0.75)");
        assert_eq!(
            s.k, 3,
            "a quadratic in one predictor has three coefficients"
        );
    }

    #[test]
    fn neighbourhood_never_exceeds_the_sample() {
        let s = Sizing::new(10, 2.0, 0.2, 2);
        assert_eq!(s.nf, 10);
    }

    #[test]
    fn degree_sets_the_coefficient_count() {
        assert_eq!(Sizing::new(53, 0.75, 0.2, 0).k, 1);
        assert_eq!(Sizing::new(53, 0.75, 0.2, 1).k, 2);
        assert_eq!(Sizing::new(53, 0.75, 0.2, 2).k, 3);
    }
}
