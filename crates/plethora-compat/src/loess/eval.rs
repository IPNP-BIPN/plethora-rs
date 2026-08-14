//! Cubic interpolation over the fitted vertices, and the whole `loess` in one
//! call.
//!
//! Transcribed from `ehg128` and `ehg133` in R's
//! `src/library/stats/src/loessf.f`, restricted to one predictor.
//!
//! `ehg128` in full also blends across cell faces for two and three predictors,
//! which is most of its length. In one dimension a cell is an interval, there
//! are no faces, and what remains is a Hermite cubic through the two endpoint
//! values and their two derivatives.

use super::Sizing;
use super::fit::{VertexFit, fit_vertices};
use super::tree::KdTree;

/// `ehg128` for one predictor: evaluate the interpolant at `z`.
///
/// The two vertex values and the two derivatives determine a unique cubic on
/// the cell. Scaling the derivative terms by the cell width is what makes the
/// interpolant continuous across cells of different widths.
///
/// # Panics
/// Panics if `z` falls outside the tree's bounding box by more than the
/// tolerance `ehg128` itself allows, which would mean the tree was built over
/// different data.
#[must_use]
pub fn interpolate(tree: &KdTree, vval: &[VertexFit], z: f64) -> f64 {
    let cell = tree.locate(z);
    let ll = tree.c[cell][1];
    let ur = tree.c[cell][tree.vc];

    let v_ll = tree.v[ll];
    let v_ur = tree.v[ur];
    let width = v_ur - v_ll;
    let h = (z - v_ll) / width;

    assert!(
        (-0.001..=1.001).contains(&h),
        "evaluation point {z} lies outside cell [{v_ll}, {v_ur}]"
    );

    let g0 = vval[ll - 1];
    let g1 = vval[ur - 1];

    // Hermite basis.
    let phi0 = (1.0 - h) * (1.0 - h) * (1.0 + 2.0 * h);
    let phi1 = h * h * (3.0 - 2.0 * h);
    let psi0 = h * (1.0 - h) * (1.0 - h);
    let psi1 = h * h * (h - 1.0);

    phi0 * g0.value + phi1 * g1.value + (psi0 * g0.derivative + psi1 * g1.derivative) * width
}

/// A fitted `loess`, holding everything `predict()` needs.
#[derive(Debug, Clone)]
pub struct Loess {
    /// The interpolation tree.
    pub tree: KdTree,
    /// Value and derivative at each vertex, indexed by vertex number minus one.
    pub vval: Vec<VertexFit>,
    /// The sizes derived from the control parameters.
    pub sizing: Sizing,
}

impl Loess {
    /// Fits `y ~ x` with R's defaults: span 0.75, degree 2, cell 0.2, one
    /// predictor, gaussian family, interpolated surface.
    ///
    /// # Panics
    /// Panics if `x` and `y` differ in length, or if `x` is empty.
    #[must_use]
    pub fn fit(x: &[f64], y: &[f64]) -> Self {
        Self::fit_with(x, y, 0.75, 2, 0.2)
    }

    /// Fits with explicit control parameters.
    ///
    /// # Panics
    /// Panics if `x` and `y` differ in length, or if `x` is empty.
    #[must_use]
    pub fn fit_with(x: &[f64], y: &[f64], span: f64, degree: usize, cell: f64) -> Self {
        assert_eq!(x.len(), y.len(), "x and y must have the same length");
        assert!(!x.is_empty(), "loess needs at least one point");

        let n = x.len();
        let sizing = Sizing::new(n, span, cell, degree);
        // R's loess_workspace uses nvmax = max(200, n), and ncmax follows it.
        let nvmax = n.max(200);

        let tree = KdTree::build(x, sizing.fc, 0.0, nvmax, nvmax);
        let vertices: Vec<f64> = (1..=tree.nv).map(|i| tree.v[i]).collect();
        let vval = fit_vertices(&vertices, x, y, sizing.nf, span, degree);

        Self { tree, vval, sizing }
    }

    /// `predict(fit)` with no new data: the interpolant at each original point.
    #[must_use]
    pub fn fitted(&self, x: &[f64]) -> Vec<f64> {
        x.iter().map(|&z| self.predict(z)).collect()
    }

    /// The interpolant at one point.
    #[must_use]
    pub fn predict(&self, z: f64) -> f64 {
        interpolate(&self.tree, &self.vval, z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Exactly linear data is in the model space at every vertex, so the whole
    /// pipeline must return the line.
    #[test]
    fn reproduces_a_line() {
        let x: Vec<f64> = (0..60).map(|i| f64::from(i) / 60.0).collect();
        let y: Vec<f64> = x.iter().map(|v| 3.0 + 2.0 * v).collect();
        let fit = Loess::fit(&x, &y);
        for (i, got) in fit.fitted(&x).iter().enumerate() {
            assert!(
                (got - y[i]).abs() < 1e-9,
                "point {i}: got {got}, want {}",
                y[i]
            );
        }
    }

    /// The interpolant must agree with the vertex fit exactly at a vertex,
    /// where the Hermite basis collapses to (1, 0, 0, 0).
    #[test]
    fn passes_through_the_vertex_values() {
        let x: Vec<f64> = (0..80).map(|i| f64::from(i) / 80.0).collect();
        let y: Vec<f64> = x.iter().map(|v| (6.0 * v).sin()).collect();
        let fit = Loess::fit(&x, &y);

        for i in 1..=fit.tree.nv {
            let v = fit.tree.v[i];
            let got = fit.predict(v);
            let want = fit.vval[i - 1].value;
            assert!(
                (got - want).abs() < 1e-12,
                "vertex {i} at {v}: interpolated {got}, fitted {want}"
            );
        }
    }

    /// A cubic interpolant with matched derivatives is continuous, so the fit
    /// must not jump across a cell boundary.
    #[test]
    fn is_continuous_across_cells() {
        let x: Vec<f64> = (0..100).map(|i| f64::from(i) / 100.0).collect();
        let y: Vec<f64> = x.iter().map(|v| (6.0 * v).sin()).collect();
        let fit = Loess::fit(&x, &y);

        for i in 1..=fit.tree.nv {
            let v = fit.tree.v[i];
            if v <= fit.tree.v[1] || v >= fit.tree.v[fit.tree.vc] {
                continue;
            }
            let eps = 1e-9;
            let left = fit.predict(v - eps);
            let right = fit.predict(v + eps);
            assert!(
                (left - right).abs() < 1e-6,
                "jump at vertex {i} ({v}): {left} against {right}"
            );
        }
    }

    /// Smoothing must actually smooth: the fit of noisy data varies less than
    /// the data.
    #[test]
    fn smooths_noise() {
        let x: Vec<f64> = (0..100).map(|i| f64::from(i) / 100.0).collect();
        // Deterministic alternating noise on a straight line.
        let y: Vec<f64> = x
            .iter()
            .enumerate()
            .map(|(i, v)| 2.0 * v + if i % 2 == 0 { 0.5 } else { -0.5 })
            .collect();
        let fit = Loess::fit(&x, &y);
        let fitted = fit.fitted(&x);

        let rough: f64 = y.windows(2).map(|w| (w[1] - w[0]).abs()).sum();
        let smooth: f64 = fitted.windows(2).map(|w| (w[1] - w[0]).abs()).sum();
        assert!(
            smooth < rough / 10.0,
            "smoothed {smooth} against raw {rough}"
        );
    }

    #[test]
    fn a_constant_stays_constant() {
        let x: Vec<f64> = (0..30).map(f64::from).collect();
        let y = vec![7.0; 30];
        let fit = Loess::fit(&x, &y);
        for got in fit.fitted(&x) {
            assert!((got - 7.0).abs() < 1e-10, "got {got}");
        }
    }
}
