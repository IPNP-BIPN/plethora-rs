//! The local weighted regression, evaluated at one point.
//!
//! Transcribed from `ehg127` in R's `src/library/stats/src/loessf.f`, with the
//! driver loop from `ehg139`.
//!
//! For each cell vertex, `loess` takes the `nf` nearest data points, weights
//! them by a tricube kernel, fits a weighted quadratic, and keeps both the
//! fitted value and its first derivative. The derivative is not a diagnostic:
//! it is what makes the interpolation between vertices cubic rather than
//! linear, so it has to come out of the same fit.

use super::blas::ddot;
use super::linpack::{dqrdc, dqrsl_qty, dsvdc};
use super::tree::select_kth;

/// `d1mach(4)`, the largest relative spacing.
///
/// `loessf.f` comments this as "1 / `DBL_EPSILON` === 2^52", which is wrong:
/// R's `d1mach.c` returns `DBL_EPSILON` for case 4. Taking the comment at face
/// value would make `tol` enormous and report every fit as singular.
const MACHEP: f64 = f64::EPSILON;

/// Leading dimension of the workspace matrices, fixed at 15 in the Fortran and
/// passed to `dsvdc` as such.
const LD: usize = 15;

/// A vertex fit: the value and the first derivative.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VertexFit {
    /// `s(0)`, the fitted value.
    pub value: f64,
    /// `s(1)`, the first derivative with respect to the predictor.
    pub derivative: f64,
    /// True when the design was rank deficient and a pseudo-inverse was used.
    pub singular: bool,
}

/// `ehg127`: fit the local regression at `q`.
///
/// `psi` carries the point ordering between calls and is permuted in place.
/// That is not an optimisation: `ehg139` initialises it once before the vertex
/// loop and lets each vertex inherit the previous vertex's permutation, and the
/// selection's result depends on its starting order whenever distances tie.
/// Reinitialising it per vertex gives a different, equally valid neighbourhood,
/// and different numbers.
///
/// # Panics
/// Panics if `nf` exceeds the sample size, or if the neighbourhood weights are
/// all zero, matching `ehg127`'s own abort.
#[must_use]
pub fn local_fit(
    q: f64,
    x: &[f64],
    y: &[f64],
    nf: usize,
    span: f64,
    degree: usize,
    psi: &mut [usize],
) -> VertexFit {
    let n = x.len();
    assert!(nf <= n && nf > 0, "neighbourhood size {nf} out of range for {n} points");

    let k = match degree {
        0 => 1,
        1 => 2,
        _ => 3,
    };

    // Squared distance to every point, indexed by the original point number.
    let mut dist = vec![0.0; n + 1];
    for i in 1..=n {
        dist[i] = (x[i - 1] - q) * (x[i - 1] - q);
    }

    select_kth(1, n, nf, &dist, psi);

    let rho = dist[psi[nf]] * span.max(1.0);
    assert!(rho > 0.0, "degenerate neighbourhood radius at q = {q}");

    // Tricube weights. The robustness weights are all 1 for family "gaussian",
    // so `sqrt(rw * ...)` collapses to `sqrt(...)`; the multiplication is kept
    // so the arithmetic matches operation for operation.
    let mut w = vec![0.0; nf + 1];
    for i in 1..=nf {
        w[i] = (dist[psi[i]] / rho).sqrt();
    }
    for i in 1..=nf {
        let u = 1.0 - w[i] * w[i] * w[i];
        w[i] = (1.0 * (u * u * u)).sqrt();
    }
    assert!(w[1..=nf].iter().any(|v| *v != 0.0), "all neighbourhood weights vanished at q = {q}");

    // Design matrix, nf by k, column-major.
    let mut b = vec![0.0; nf * k];
    b[..nf].copy_from_slice(&w[1..=nf]);
    if degree >= 1 {
        for i in 1..=nf {
            b[nf + i - 1] = w[i] * (x[psi[i] - 1] - q);
        }
    }
    if degree >= 2 {
        for i in 1..=nf {
            let dx = x[psi[i] - 1] - q;
            b[2 * nf + i - 1] = w[i] * dx * dx;
        }
    }

    let mut eta = vec![0.0; nf];
    for i in 1..=nf {
        eta[i - 1] = w[i] * y[psi[i] - 1];
    }

    // Equilibrate the columns. Written as a plain accumulation and a sqrt, as
    // in the Fortran, rather than through dnrm2, which would scale differently.
    let mut colnor = vec![1.0; k + 1];
    for j in 1..=k {
        let mut scal = 0.0;
        for i in 0..nf {
            scal += b[(j - 1) * nf + i] * b[(j - 1) * nf + i];
        }
        scal = scal.sqrt();
        if scal > 0.0 {
            for i in 0..nf {
                b[(j - 1) * nf + i] /= scal;
            }
            colnor[j] = scal;
        }
    }

    let mut qraux = vec![0.0; k];
    dqrdc(&mut b, nf, nf, k, &mut qraux, 0);

    let mut qty = vec![0.0; nf];
    dqrsl_qty(&mut b, nf, nf, k, &qraux, &eta, &mut qty);
    eta.copy_from_slice(&qty);

    // The upper triangle of R, in a 15 by 15 workspace as dsvdc expects.
    let mut u = vec![0.0; LD * LD];
    for i in 1..=k {
        for j in i..=k {
            u[(j - 1) * LD + (i - 1)] = b[(j - 1) * nf + (i - 1)];
        }
    }

    let mut sigma = vec![0.0; LD + 1];
    let mut g = vec![0.0; LD];
    let mut e = vec![0.0; LD * LD];
    let mut work = vec![0.0; LD];
    let info = dsvdc(
        &mut u.clone(),
        LD,
        k,
        k,
        &mut sigma,
        &mut g,
        &mut u,
        LD,
        &mut e,
        LD,
        &mut work,
        21,
    );
    assert_eq!(info, 0, "the singular value decomposition did not converge at q = {q}");

    let tol = sigma[0] * (100.0 * MACHEP);
    let singular = sigma[k - 1] <= tol;

    // Undo the equilibration, row by row.
    for j in 1..=k {
        let scale = colnor[j];
        for i in 1..=k {
            e[(i - 1) * LD + (j - 1)] /= scale;
        }
    }

    // Solve the least squares problem through the pseudo-inverse.
    let mut dgamma = vec![0.0; k];
    for j in 1..=k {
        dgamma[j - 1] = if tol < sigma[j - 1] {
            ddot(k, &u[(j - 1) * LD..], 1, &eta, 1) / sigma[j - 1]
        } else {
            0.0
        };
    }

    // s(j) reads row j+1 of e, hence the stride of LD.
    let coefficient = |j: usize| -> f64 {
        if j < k {
            ddot(k, &e[j..], LD as isize, &dgamma, 1)
        } else {
            0.0
        }
    };

    VertexFit {
        value: coefficient(0),
        derivative: coefficient(1),
        singular,
    }
}

/// `ehg139`, reduced to what `predict()` needs: fit every vertex in order.
///
/// The `trl` and `setlf` branches are not carried; they build the hat matrix
/// for `one.delta`, `two.delta` and standard errors, none of which
/// `gc_correction.R` reads.
#[must_use]
pub fn fit_vertices(
    vertices: &[f64],
    x: &[f64],
    y: &[f64],
    nf: usize,
    span: f64,
    degree: usize,
) -> Vec<VertexFit> {
    // Initialised once, then inherited from vertex to vertex. See local_fit.
    let mut psi: Vec<usize> = (0..=x.len()).collect();

    vertices
        .iter()
        .map(|&q| local_fit(q, x, y, nf, span, degree, &mut psi))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// On exactly linear data a quadratic local fit must return the line, and
    /// its slope, to rounding.
    #[test]
    fn recovers_a_line() {
        let x: Vec<f64> = (0..40).map(|i| f64::from(i) / 40.0).collect();
        let y: Vec<f64> = x.iter().map(|v| 3.0 + 2.0 * v).collect();
        let mut psi: Vec<usize> = (0..=x.len()).collect();

        let fit = local_fit(0.5, &x, &y, 30, 0.75, 2, &mut psi);
        assert!((fit.value - 4.0).abs() < 1e-10, "value {}", fit.value);
        assert!((fit.derivative - 2.0).abs() < 1e-9, "derivative {}", fit.derivative);
        assert!(!fit.singular);
    }

    /// A quadratic is in the model space, so degree 2 must reproduce it too.
    #[test]
    fn recovers_a_parabola() {
        let x: Vec<f64> = (0..60).map(|i| f64::from(i) / 60.0).collect();
        let y: Vec<f64> = x.iter().map(|v| 1.0 - 4.0 * (v - 0.5) * (v - 0.5)).collect();
        let mut psi: Vec<usize> = (0..=x.len()).collect();

        let fit = local_fit(0.5, &x, &y, 45, 0.75, 2, &mut psi);
        assert!((fit.value - 1.0).abs() < 1e-10, "value {}", fit.value);
        assert!(fit.derivative.abs() < 1e-9, "derivative at the apex {}", fit.derivative);
    }

    /// Constant data gives a constant fit and a zero slope.
    #[test]
    fn a_constant_has_no_slope() {
        let x: Vec<f64> = (0..20).map(f64::from).collect();
        let y = vec![7.0; 20];
        let mut psi: Vec<usize> = (0..=x.len()).collect();

        let fit = local_fit(10.0, &x, &y, 15, 0.75, 2, &mut psi);
        assert!((fit.value - 7.0).abs() < 1e-12);
        assert!(fit.derivative.abs() < 1e-12);
    }

    /// Degree 1 cannot follow curvature, so it must differ from degree 2 on a
    /// parabola. Guards against the degree argument being ignored.
    #[test]
    fn degree_changes_the_fit() {
        let x: Vec<f64> = (0..60).map(|i| f64::from(i) / 60.0).collect();
        let y: Vec<f64> = x.iter().map(|v| 1.0 - 4.0 * (v - 0.5) * (v - 0.5)).collect();

        let mut psi1: Vec<usize> = (0..=x.len()).collect();
        let linear = local_fit(0.5, &x, &y, 45, 0.75, 1, &mut psi1);
        let mut psi2: Vec<usize> = (0..=x.len()).collect();
        let quadratic = local_fit(0.5, &x, &y, 45, 0.75, 2, &mut psi2);

        assert!((linear.value - quadratic.value).abs() > 1e-6);
    }

    /// The permutation must survive as a permutation across vertices, since it
    /// is threaded rather than rebuilt.
    #[test]
    fn the_point_ordering_stays_a_permutation() {
        let x: Vec<f64> = (0..30).map(|i| f64::from(i) / 30.0).collect();
        let y: Vec<f64> = x.iter().map(|v| v.sin()).collect();
        let mut psi: Vec<usize> = (0..=x.len()).collect();

        for q in [0.1, 0.5, 0.9, 0.3] {
            let _ = local_fit(q, &x, &y, 22, 0.75, 2, &mut psi);
            let mut seen: Vec<usize> = psi[1..].to_vec();
            seen.sort_unstable();
            assert_eq!(seen, (1..=30).collect::<Vec<_>>());
        }
    }

    #[test]
    fn fit_vertices_returns_one_fit_per_vertex() {
        let x: Vec<f64> = (0..40).map(|i| f64::from(i) / 40.0).collect();
        let y: Vec<f64> = x.iter().map(|v| 3.0 + 2.0 * v).collect();
        let fits = fit_vertices(&[0.0, 0.25, 0.5, 0.75, 1.0], &x, &y, 30, 0.75, 2);
        assert_eq!(fits.len(), 5);
        for (i, f) in fits.iter().enumerate() {
            let q = f64::from(i as u32) * 0.25;
            assert!((f.value - (3.0 + 2.0 * q)).abs() < 1e-9);
        }
    }
}
