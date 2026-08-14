//! The BLAS calls LINPACK makes, in their reference form.
//!
//! `ehg127` fits each local regression through LINPACK's `dqrdc`, `dqrsl` and
//! `dsvdc`, which in turn call BLAS. Which BLAS is therefore part of the answer.
//!
//! These are transcriptions of the reference implementations R ships in
//! `src/extra/blas` (`blas.f` for the level-1 kernels, `blas2.f90` for `dnrm2`
//! and `drotg`, which moved to Fortran 90 with `LAPACK` 3.10's safe-scaling
//! rewrite). Reference BLAS is what a stock R build uses, and what the cluster
//! behind the 2017 paper would have had.
//!
//! An R linked against an optimised BLAS gives slightly different answers. The
//! R on this machine links `OpenBLAS` 0.3.34, whose NEON dot product accumulates
//! into several partial sums; measured on a 39-element dot product, the size
//! `ehg127`'s neighbourhood actually uses, it lands 8 ULP from sequential
//! summation. That is a property of the installation, not of this port, and it
//! is recorded in `DIVERGENCES.md`.
//!
//! Note on the unrolling: reference `ddot` is unrolled by five, but the
//! expression `dtemp + a + b + c + d + e` associates left, so the accumulation
//! order is exactly sequential and a plain loop is bit-identical to it.

/// `ddot`: the dot product, accumulated strictly left to right.
#[must_use]
pub fn ddot(n: usize, dx: &[f64], incx: isize, dy: &[f64], incy: isize) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let mut dtemp = 0.0;
    let mut ix = start_index(n, incx);
    let mut iy = start_index(n, incy);
    for _ in 0..n {
        dtemp += dx[ix] * dy[iy];
        ix = ix.wrapping_add_signed(incx);
        iy = iy.wrapping_add_signed(incy);
    }
    dtemp
}

/// `dnrm2`: the Euclidean norm, by the three-accumulator safe-scaling scheme.
///
/// Values are binned into small, medium and large by magnitude so that squaring
/// cannot overflow or flush to zero. The binning is not a refinement over a
/// naive `sqrt(sum of squares)`: it changes the result in the last bits for
/// ordinary inputs too, because the medium bin is summed on its own.
#[must_use]
pub fn dnrm2(n: usize, x: &[f64], incx: isize) -> f64 {
    // The thresholds are radix and exponent expressions in the Fortran; for
    // IEEE binary64 they evaluate to these powers of two.
    /// `2^ceiling((minexponent - 1) / 2)`.
    const TSML: f64 = 1.4916681462400413e-154;
    /// `2^floor((maxexponent - digits + 1) / 2)`.
    const TBIG: f64 = 1.9979698819072276e146;
    /// `2^-floor((minexponent - digits) / 2)`.
    const SSML: f64 = 4.4989137945431964e161;
    /// `2^-ceiling((maxexponent + digits - 1) / 2)`.
    const SBIG: f64 = 1.1113793747425387e-162;

    if n == 0 {
        return 0.0;
    }

    let mut notbig = true;
    let mut asml = 0.0;
    let mut amed = 0.0;
    let mut abig = 0.0;

    let mut ix = start_index(n, incx);
    for _ in 0..n {
        let ax = x[ix].abs();
        if ax > TBIG {
            abig += (ax * SBIG) * (ax * SBIG);
            notbig = false;
        } else if ax < TSML {
            if notbig {
                asml += (ax * SSML) * (ax * SSML);
            }
        } else {
            amed += ax * ax;
        }
        ix = ix.wrapping_add_signed(incx);
    }

    let (scl, sumsq) = if abig > 0.0 {
        if amed > 0.0 || amed.is_nan() {
            abig += (amed * SBIG) * SBIG;
        }
        (1.0 / SBIG, abig)
    } else if asml > 0.0 {
        if amed > 0.0 || amed.is_nan() {
            let amed = amed.sqrt();
            let asml = asml.sqrt() / SSML;
            let (ymin, ymax) = if asml > amed {
                (amed, asml)
            } else {
                (asml, amed)
            };
            (1.0, ymax * ymax * (1.0 + (ymin / ymax) * (ymin / ymax)))
        } else {
            (1.0 / SSML, asml)
        }
    } else {
        (1.0, amed)
    };

    scl * sumsq.sqrt()
}

/// `daxpy`: `y := alpha * x + y`.
pub fn daxpy(n: usize, da: f64, dx: &[f64], incx: isize, dy: &mut [f64], incy: isize) {
    if n == 0 || da == 0.0 {
        return;
    }
    let mut ix = start_index(n, incx);
    let mut iy = start_index(n, incy);
    for _ in 0..n {
        dy[iy] += da * dx[ix];
        ix = ix.wrapping_add_signed(incx);
        iy = iy.wrapping_add_signed(incy);
    }
}

/// `dscal`: `x := alpha * x`.
pub fn dscal(n: usize, da: f64, dx: &mut [f64], incx: isize) {
    if n == 0 {
        return;
    }
    let mut ix = start_index(n, incx);
    for _ in 0..n {
        dx[ix] *= da;
        ix = ix.wrapping_add_signed(incx);
    }
}

/// `dcopy`: `y := x`.
pub fn dcopy(n: usize, dx: &[f64], incx: isize, dy: &mut [f64], incy: isize) {
    if n == 0 {
        return;
    }
    let mut ix = start_index(n, incx);
    let mut iy = start_index(n, incy);
    for _ in 0..n {
        dy[iy] = dx[ix];
        ix = ix.wrapping_add_signed(incx);
        iy = iy.wrapping_add_signed(incy);
    }
}

/// `idamax`: the index of the largest absolute value, 0-based here.
///
/// Ties go to the first, as in the Fortran, which matters when `ehg124` picks a
/// split dimension and two dimensions have the same spread.
#[must_use]
pub fn idamax(n: usize, dx: &[f64], incx: isize) -> usize {
    if n == 0 {
        return 0;
    }
    let mut ix = start_index(n, incx);
    let mut dmax = dx[ix].abs();
    let mut best = 0;
    for i in 1..n {
        ix = ix.wrapping_add_signed(incx);
        if dx[ix].abs() > dmax {
            best = i;
            dmax = dx[ix].abs();
        }
    }
    best
}

/// `drotg`: construct a Givens rotation.
///
/// Returns `(r, z, c, s)` for the LINPACK convention, where `a` becomes `r` and
/// `b` becomes `z`. The safe-scaling form from `LAPACK` 3.10, not the older
/// `roe`-based one, which differs in sign conventions on the boundary cases.
#[must_use]
pub fn drotg(a: f64, b: f64) -> (f64, f64, f64, f64) {
    /// `radix^max(minexponent - 1, 1 - maxexponent)`.
    const SAFMIN: f64 = 2.2250738585072014e-308;
    /// `radix^max(1 - minexponent, maxexponent - 1)`.
    const SAFMAX: f64 = 8.98846567431158e307;

    let anorm = a.abs();
    let bnorm = b.abs();

    if bnorm == 0.0 {
        (a, 0.0, 1.0, 0.0)
    } else if anorm == 0.0 {
        (b, 1.0, 0.0, 1.0)
    } else {
        let scl = SAFMAX.min(SAFMIN.max(anorm.max(bnorm)));
        let sigma = if anorm > bnorm {
            sign(1.0, a)
        } else {
            sign(1.0, b)
        };
        let r = sigma * (scl * ((a / scl) * (a / scl) + (b / scl) * (b / scl)).sqrt());
        let c = a / r;
        let s = b / r;
        let z = if anorm > bnorm {
            s
        } else if c != 0.0 {
            1.0 / c
        } else {
            1.0
        };
        (r, z, c, s)
    }
}

/// `drot`: apply a Givens rotation to a pair of vectors.
pub fn drot(n: usize, dx: &mut [f64], incx: isize, dy: &mut [f64], incy: isize, c: f64, s: f64) {
    if n == 0 {
        return;
    }
    let mut ix = start_index(n, incx);
    let mut iy = start_index(n, incy);
    for _ in 0..n {
        let dtemp = c * dx[ix] + s * dy[iy];
        dy[iy] = c * dy[iy] - s * dx[ix];
        dx[ix] = dtemp;
        ix = ix.wrapping_add_signed(incx);
        iy = iy.wrapping_add_signed(incy);
    }
}

/// `dswap`: exchange two vectors.
pub fn dswap(n: usize, dx: &mut [f64], incx: isize, dy: &mut [f64], incy: isize) {
    if n == 0 {
        return;
    }
    let mut ix = start_index(n, incx);
    let mut iy = start_index(n, incy);
    for _ in 0..n {
        std::mem::swap(&mut dx[ix], &mut dy[iy]);
        ix = ix.wrapping_add_signed(incx);
        iy = iy.wrapping_add_signed(incy);
    }
}

/// Fortran's `sign(a, b)`: the magnitude of `a` with the sign of `b`.
///
/// Not `f64::copysign`, which treats -0.0 as negative; Fortran's `sign` gives a
/// positive result for a `b` of -0.0 under the default (non-IEEE) reading, and
/// LINPACK relies on that in `dqrdc`'s pivot sign.
#[must_use]
pub fn sign(a: f64, b: f64) -> f64 {
    if b >= 0.0 { a.abs() } else { -a.abs() }
}

/// Where a strided walk begins, mirroring `if (incx < 0) ix = 1 - (n-1)*incx`.
#[inline]
fn start_index(n: usize, inc: isize) -> usize {
    if inc < 0 {
        (n - 1).wrapping_mul(inc.unsigned_abs())
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ddot_accumulates_left_to_right() {
        // A sum where order is visible: the small terms vanish if added last.
        let a = [1.0, 1e-17, 1e-17];
        let b = [1.0, 1.0, 1.0];
        let mut expected = 0.0;
        for i in 0..3 {
            expected += a[i] * b[i];
        }
        assert_eq!(ddot(3, &a, 1, &b, 1).to_bits(), expected.to_bits());
    }

    #[test]
    fn ddot_handles_strides() {
        let a = [1.0, 9.0, 2.0, 9.0, 3.0];
        let b = [1.0, 1.0, 1.0];
        assert_eq!(ddot(3, &a, 2, &b, 1), 6.0);
    }

    #[test]
    fn dnrm2_matches_the_obvious_answer_when_it_is_safe() {
        assert_eq!(dnrm2(2, &[3.0, 4.0], 1), 5.0);
        assert_eq!(dnrm2(0, &[], 1), 0.0);
        assert_eq!(dnrm2(1, &[-7.0], 1), 7.0);
    }

    /// The point of the three-bin scheme: neither overflow nor flush to zero.
    #[test]
    fn dnrm2_survives_extreme_magnitudes() {
        let big = dnrm2(2, &[1e300, 1e300], 1);
        assert!(big.is_finite() && big > 1e300);
        let small = dnrm2(2, &[1e-300, 1e-300], 1);
        assert!(small > 0.0 && small < 1e-299);
    }

    #[test]
    fn idamax_takes_the_first_of_equal_maxima() {
        assert_eq!(idamax(4, &[1.0, 3.0, -3.0, 2.0], 1), 1);
        assert_eq!(idamax(3, &[-5.0, 1.0, 2.0], 1), 0);
    }

    /// A Givens rotation must zero the second component.
    #[test]
    fn drotg_zeroes_the_second_component() {
        for (a, b) in [(3.0, 4.0), (-3.0, 4.0), (3.0, -4.0), (1e-300, 1e300)] {
            let (r, _z, c, s) = drotg(a, b);
            assert!((c * a + s * b - r).abs() <= 1e-9 * r.abs().max(1.0));
            assert!((c * b - s * a).abs() <= 1e-9 * r.abs().max(1.0));
        }
    }

    #[test]
    fn drotg_handles_the_degenerate_cases() {
        assert_eq!(drotg(2.0, 0.0), (2.0, 0.0, 1.0, 0.0));
        assert_eq!(drotg(0.0, 2.0), (2.0, 1.0, 0.0, 1.0));
    }

    /// Fortran's sign, not copysign: a `b` of -0.0 counts as positive.
    #[test]
    fn sign_treats_negative_zero_as_positive() {
        assert_eq!(sign(3.0, -0.0), 3.0);
        assert_eq!(sign(3.0, 0.0), 3.0);
        assert_eq!(sign(3.0, -1.0), -3.0);
        assert_eq!(sign(-3.0, 1.0), 3.0);
    }

    #[test]
    fn daxpy_and_dscal_are_elementwise() {
        let x = [1.0, 2.0, 3.0];
        let mut y = [10.0, 20.0, 30.0];
        daxpy(3, 2.0, &x, 1, &mut y, 1);
        assert_eq!(y, [12.0, 24.0, 36.0]);
        dscal(3, 0.5, &mut y, 1);
        assert_eq!(y, [6.0, 12.0, 18.0]);
    }

    #[test]
    fn daxpy_with_zero_alpha_is_a_no_op() {
        let x = [f64::NAN, f64::NAN];
        let mut y = [1.0, 2.0];
        daxpy(2, 0.0, &x, 1, &mut y, 1);
        assert_eq!(y, [1.0, 2.0]);
    }

    #[test]
    fn dswap_and_dcopy_move_the_right_elements() {
        let mut a = [1.0, 2.0, 3.0];
        let mut b = [4.0, 5.0, 6.0];
        dswap(3, &mut a, 1, &mut b, 1);
        assert_eq!(a, [4.0, 5.0, 6.0]);
        assert_eq!(b, [1.0, 2.0, 3.0]);

        let mut c = [0.0; 3];
        dcopy(3, &a, 1, &mut c, 1);
        assert_eq!(c, a);
    }

    #[test]
    fn drot_rotates_in_place() {
        let mut x = [1.0, 0.0];
        let mut y = [0.0, 1.0];
        drot(2, &mut x, 1, &mut y, 1, 0.0, 1.0);
        assert_eq!(x, [0.0, 1.0]);
        assert_eq!(y, [-1.0, 0.0]);
    }
}
