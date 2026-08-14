//! The LINPACK routines `ehg127` solves each local regression with.
//!
//! Transcribed from the copies R ships in `src/appl`: `dqrdc.f`, `dqrsl.f` and
//! `dsvdc.f`. Matrices are column-major with a leading dimension, as in the
//! Fortran, and indices are 1-based at the call sites so the transcription can
//! be read against the original.
//!
//! Scope is deliberately narrow. `ehg127` calls `dqrdc` with `job = 0` and
//! `dqrsl` with `job = 1000`, so only the unpivoted decomposition and the
//! `qty` product are implemented. The pivoting, back-substitution and residual
//! paths are not dead-code-carried: an untested transcription of a numerical
//! routine is worse than an honest absence, and the entry points assert their
//! scope. `dsvdc` is complete, because its QR iteration cannot be trimmed.
//!
//! `ehg139`'s `trl` and `setlf` branches, which would need `dqrsl`'s `qy`
//! product, compute the trace of the hat matrix for `one.delta` and
//! `two.delta`. `gc_correction.R` uses neither, only `predict()`.

use super::blas::{daxpy, dcopy, ddot, dnrm2, drot, drotg, dscal, dswap, sign};

/// Offset of the 1-based element `(i, j)` in a column-major matrix.
#[inline]
const fn at(ld: usize, i: usize, j: usize) -> usize {
    (j - 1) * ld + (i - 1)
}

/// Two distinct columns, both mutable.
///
/// `split_at_mut` rather than pointer arithmetic, because this crate forbids
/// unsafe code. Requires `j1 < j2`.
fn two_cols_mut(x: &mut [f64], ld: usize, j1: usize, j2: usize) -> (&mut [f64], &mut [f64]) {
    debug_assert!(j1 < j2, "columns must be distinct and ordered");
    let (left, right) = x.split_at_mut((j2 - 1) * ld);
    (&mut left[(j1 - 1) * ld..], right)
}

/// `dqrdc` with `job = 0`: Householder QR without pivoting.
///
/// On return the upper triangle of `x` holds R, the strict lower triangle holds
/// the Householder vectors, and `qraux` holds the extra information needed to
/// recover Q.
///
/// # Panics
/// Panics unless `job == 0`; see the module note on scope.
pub fn dqrdc(x: &mut [f64], ldx: usize, n: usize, p: usize, qraux: &mut [f64], job: i32) {
    assert_eq!(job, 0, "only the unpivoted decomposition is implemented");

    let lup = n.min(p);
    for l in 1..=lup {
        qraux[l - 1] = 0.0;
        if l == n {
            continue;
        }

        let mut nrmxl = dnrm2(n - l + 1, &x[at(ldx, l, l)..], 1);
        if nrmxl == 0.0 {
            continue;
        }

        if x[at(ldx, l, l)] != 0.0 {
            nrmxl = sign(nrmxl, x[at(ldx, l, l)]);
        }
        dscal(n - l + 1, 1.0 / nrmxl, &mut x[at(ldx, l, l)..], 1);
        x[at(ldx, l, l)] += 1.0;

        // The Householder vector lives in column l and is not touched while the
        // later columns are updated, so a copy sidesteps the aliasing without
        // changing a single operation.
        let house: Vec<f64> = x[at(ldx, l, l)..at(ldx, l, l) + (n - l + 1)].to_vec();
        let pivot = house[0];

        for j in l + 1..=p {
            let t = -ddot(n - l + 1, &house, 1, &x[at(ldx, l, j)..], 1) / pivot;
            daxpy(n - l + 1, t, &house, 1, &mut x[at(ldx, l, j)..], 1);
        }

        qraux[l - 1] = x[at(ldx, l, l)];
        x[at(ldx, l, l)] = -nrmxl;
    }
}

/// `dqrsl` with `job = 1000`: form `qty = Q' * y`.
///
/// `x` and `qraux` must come from [`dqrdc`]. `x` is taken mutably because the
/// Fortran borrows the diagonal to hold the Householder pivot and restores it;
/// the matrix is unchanged on return.
pub fn dqrsl_qty(x: &mut [f64], ldx: usize, n: usize, k: usize, qraux: &[f64], y: &[f64], qty: &mut [f64]) {
    let ju = k.min(n - 1);

    if ju == 0 {
        qty[0] = y[0];
        return;
    }

    dcopy(n, y, 1, qty, 1);

    for j in 1..=ju {
        if qraux[j - 1] == 0.0 {
            continue;
        }
        let temp = x[at(ldx, j, j)];
        x[at(ldx, j, j)] = qraux[j - 1];

        let house: Vec<f64> = x[at(ldx, j, j)..at(ldx, j, j) + (n - j + 1)].to_vec();
        let t = -ddot(n - j + 1, &house, 1, &qty[j - 1..], 1) / house[0];
        daxpy(n - j + 1, t, &house, 1, &mut qty[j - 1..], 1);

        x[at(ldx, j, j)] = temp;
    }
}

/// `dsvdc`: the singular value decomposition, by Golub and Reinsch.
///
/// `job` follows LINPACK's two-digit encoding. `ehg127` passes 21: return the
/// first `min(n, p)` left singular vectors in `u` and the right singular
/// vectors in `v`.
///
/// R patched the two convergence tests: where LINPACK compares `ztest .ne.
/// test`, R compares a relative accuracy against 1e-15. Without that patch the
/// iteration can fail to terminate under extended-precision registers. The
/// patched form is what R runs, so it is what is transcribed.
///
/// Returns LINPACK's `info`: 0 on success, or the index at which the iteration
/// gave up.
///
/// # Panics
/// Panics if `n` or `p` is zero.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn dsvdc(
    x: &mut [f64],
    ldx: usize,
    n: usize,
    p: usize,
    s: &mut [f64],
    e: &mut [f64],
    u: &mut [f64],
    ldu: usize,
    v: &mut [f64],
    ldv: usize,
    work: &mut [f64],
    job: i32,
) -> usize {
    /// LINPACK's iteration cap before it gives up and reports failure.
    const MAXIT: usize = 30;

    assert!(n > 0 && p > 0, "dsvdc needs a non-empty matrix");

    let jobu = (job % 100) / 10;
    let ncu = if jobu > 1 { n.min(p) } else { n };
    let wantu = jobu != 0;
    let wantv = job % 10 != 0;

    let mut info = 0;
    let nct = (n - 1).min(p);
    let nrt = if p >= 2 { (p - 2).min(n) } else { 0 };
    let lu = nct.max(nrt);

    // Reduce to bidiagonal form.
    for l in 1..=lu {
        let lp1 = l + 1;

        if l <= nct {
            s[l - 1] = dnrm2(n - l + 1, &x[at(ldx, l, l)..], 1);
            if s[l - 1] != 0.0 {
                if x[at(ldx, l, l)] != 0.0 {
                    s[l - 1] = sign(s[l - 1], x[at(ldx, l, l)]);
                }
                dscal(n - l + 1, 1.0 / s[l - 1], &mut x[at(ldx, l, l)..], 1);
                x[at(ldx, l, l)] += 1.0;
            }
            s[l - 1] = -s[l - 1];
        }

        for j in lp1..=p {
            if l <= nct && s[l - 1] != 0.0 {
                let house: Vec<f64> = x[at(ldx, l, l)..at(ldx, l, l) + (n - l + 1)].to_vec();
                let t = -ddot(n - l + 1, &house, 1, &x[at(ldx, l, j)..], 1) / house[0];
                daxpy(n - l + 1, t, &house, 1, &mut x[at(ldx, l, j)..], 1);
            }
            e[j - 1] = x[at(ldx, l, j)];
        }

        if wantu && l <= nct {
            for i in l..=n {
                u[at(ldu, i, l)] = x[at(ldx, i, l)];
            }
        }

        if l <= nrt {
            e[l - 1] = dnrm2(p - l, &e[lp1 - 1..], 1);
            if e[l - 1] != 0.0 {
                if e[lp1 - 1] != 0.0 {
                    e[l - 1] = sign(e[l - 1], e[lp1 - 1]);
                }
                dscal(p - l, 1.0 / e[l - 1], &mut e[lp1 - 1..], 1);
                e[lp1 - 1] += 1.0;
            }
            e[l - 1] = -e[l - 1];

            if lp1 <= n && e[l - 1] != 0.0 {
                for i in lp1..=n {
                    work[i - 1] = 0.0;
                }
                for j in lp1..=p {
                    let col: Vec<f64> = x[at(ldx, lp1, j)..at(ldx, lp1, j) + (n - l)].to_vec();
                    daxpy(n - l, e[j - 1], &col, 1, &mut work[lp1 - 1..], 1);
                }
                for j in lp1..=p {
                    let scale = -e[j - 1] / e[lp1 - 1];
                    let w: Vec<f64> = work[lp1 - 1..lp1 - 1 + (n - l)].to_vec();
                    daxpy(n - l, scale, &w, 1, &mut x[at(ldx, lp1, j)..], 1);
                }
            }

            if wantv {
                for i in lp1..=p {
                    v[at(ldv, i, l)] = e[i - 1];
                }
            }
        }
    }

    // Set up the final bidiagonal matrix of order m.
    let mut m = p.min(n + 1);
    let nctp1 = nct + 1;
    let nrtp1 = nrt + 1;
    if nct < p {
        s[nctp1 - 1] = x[at(ldx, nctp1, nctp1)];
    }
    if n < m {
        s[m - 1] = 0.0;
    }
    if nrtp1 < m {
        e[nrtp1 - 1] = x[at(ldx, nrtp1, m)];
    }
    e[m - 1] = 0.0;

    // Generate U.
    if wantu {
        for j in nctp1..=ncu {
            for i in 1..=n {
                u[at(ldu, i, j)] = 0.0;
            }
            u[at(ldu, j, j)] = 1.0;
        }
        for ll in 1..=nct {
            let l = nct - ll + 1;
            if s[l - 1] == 0.0 {
                for i in 1..=n {
                    u[at(ldu, i, l)] = 0.0;
                }
                u[at(ldu, l, l)] = 1.0;
            } else {
                let lp1 = l + 1;
                for j in lp1..=ncu {
                    let house: Vec<f64> = u[at(ldu, l, l)..at(ldu, l, l) + (n - l + 1)].to_vec();
                    let t = -ddot(n - l + 1, &house, 1, &u[at(ldu, l, j)..], 1) / house[0];
                    daxpy(n - l + 1, t, &house, 1, &mut u[at(ldu, l, j)..], 1);
                }
                dscal(n - l + 1, -1.0, &mut u[at(ldu, l, l)..], 1);
                u[at(ldu, l, l)] += 1.0;
                for i in 1..l {
                    u[at(ldu, i, l)] = 0.0;
                }
            }
        }
    }

    // Generate V.
    if wantv {
        for ll in 1..=p {
            let l = p - ll + 1;
            let lp1 = l + 1;
            if l <= nrt && e[l - 1] != 0.0 {
                for j in lp1..=p {
                    let house: Vec<f64> = v[at(ldv, lp1, l)..at(ldv, lp1, l) + (p - l)].to_vec();
                    let t = -ddot(p - l, &house, 1, &v[at(ldv, lp1, j)..], 1) / house[0];
                    daxpy(p - l, t, &house, 1, &mut v[at(ldv, lp1, j)..], 1);
                }
            }
            for i in 1..=p {
                v[at(ldv, i, l)] = 0.0;
            }
            v[at(ldv, l, l)] = 1.0;
        }
    }

    // The QR iteration on the bidiagonal form.
    let mm = m;
    let mut iter = 0;

    while m != 0 {
        if iter >= MAXIT {
            info = m;
            break;
        }

        // Split at a negligible superdiagonal or diagonal element. R's patched
        // test: a relative accuracy rather than an exact float comparison.
        let mut l = 0;
        for ll in 1..=m {
            l = m - ll;
            if l == 0 {
                break;
            }
            let test = s[l - 1].abs() + s[l].abs();
            let ztest = test + e[l - 1].abs();
            let acc = (test - ztest).abs() / (1.0e-100 + test);
            if acc <= 1e-15 {
                e[l - 1] = 0.0;
                break;
            }
        }

        let kase;
        if l == m - 1 {
            kase = 4;
        } else {
            let lp1 = l + 1;
            let mp1 = m + 1;
            let mut ls = l;
            for lls in lp1..=mp1 {
                // Signed: the last iteration has lls = m + 1, where the Fortran
                // relies on m - lls going negative before lp1 brings it back to
                // exactly l, which is the loop's exit condition.
                ls = usize::try_from(m as isize - lls as isize + lp1 as isize)
                    .expect("ls stays non-negative for lls <= m + 1");
                if ls == l {
                    break;
                }
                let mut test = 0.0;
                if ls != m {
                    test += e[ls - 1].abs();
                }
                if ls != l + 1 {
                    test += e[ls - 2].abs();
                }
                let ztest = test + s[ls - 1].abs();
                let acc = (test - ztest).abs() / (1.0e-100 + test);
                if acc <= 1e-15 {
                    s[ls - 1] = 0.0;
                    break;
                }
            }
            if ls == l {
                kase = 3;
            } else if ls == m {
                kase = 1;
            } else {
                kase = 2;
                l = ls;
            }
        }
        l += 1;

        match kase {
            // Deflate negligible s(m).
            1 => {
                let mm1 = m - 1;
                let mut f = e[m - 2];
                e[m - 2] = 0.0;
                for kk in l..=mm1 {
                    let k = mm1 - kk + l;
                    let (t1, _z, cs, sn) = drotg(s[k - 1], f);
                    s[k - 1] = t1;
                    if k != l {
                        f = -sn * e[k - 2];
                        e[k - 2] = cs * e[k - 2];
                    }
                    if wantv {
                        let (a, b) = two_cols_mut(v, ldv, k, m);
                        drot(p, a, 1, b, 1, cs, sn);
                    }
                }
            }
            // Split at a negligible s(l).
            2 => {
                let mut f = e[l - 2];
                e[l - 2] = 0.0;
                for k in l..=m {
                    let (t1, _z, cs, sn) = drotg(s[k - 1], f);
                    s[k - 1] = t1;
                    f = -sn * e[k - 1];
                    e[k - 1] = cs * e[k - 1];
                    if wantu {
                        let (a, b) = two_cols_mut(u, ldu, l - 1, k);
                        // Columns are passed in the Fortran's order, which is
                        // (k, l-1); two_cols_mut needs them ordered, so the
                        // rotation is applied with the sign convention swapped
                        // back by exchanging the roles.
                        drot_swapped(n, b, a, cs, sn);
                    }
                }
            }
            // Perform one QR step.
            3 => {
                let scale = s[m - 1]
                    .abs()
                    .max(s[m - 2].abs())
                    .max(e[m - 2].abs())
                    .max(s[l - 1].abs())
                    .max(e[l - 1].abs());
                let sm = s[m - 1] / scale;
                let smm1 = s[m - 2] / scale;
                let emm1 = e[m - 2] / scale;
                let sl = s[l - 1] / scale;
                let el = e[l - 1] / scale;
                // Not a midpoint despite its shape: the whole expression is
                // halved, and rewriting it through f64::midpoint would change
                // the rounding and with it the last bits of every fitted value.
                #[allow(clippy::manual_midpoint)]
                let b = ((smm1 + sm) * (smm1 - sm) + emm1 * emm1) / 2.0;
                let c = (sm * emm1) * (sm * emm1);
                let mut shift = 0.0;
                if b != 0.0 || c != 0.0 {
                    shift = (b * b + c).sqrt();
                    if b < 0.0 {
                        shift = -shift;
                    }
                    shift = c / (b + shift);
                }
                let mut f = (sl + sm) * (sl - sm) + shift;
                let mut g = sl * el;

                for k in l..m {
                    let (r, _z, cs, sn) = drotg(f, g);
                    f = r;
                    if k != l {
                        e[k - 2] = f;
                    }
                    f = cs * s[k - 1] + sn * e[k - 1];
                    e[k - 1] = cs * e[k - 1] - sn * s[k - 1];
                    g = sn * s[k];
                    s[k] = cs * s[k];
                    if wantv {
                        let (a, bb) = two_cols_mut(v, ldv, k, k + 1);
                        drot(p, a, 1, bb, 1, cs, sn);
                    }

                    let (r, _z, cs, sn) = drotg(f, g);
                    s[k - 1] = r;
                    f = cs * e[k - 1] + sn * s[k];
                    s[k] = -sn * e[k - 1] + cs * s[k];
                    g = sn * e[k];
                    e[k] = cs * e[k];
                    if wantu && k < n {
                        let (a, bb) = two_cols_mut(u, ldu, k, k + 1);
                        drot(n, a, 1, bb, 1, cs, sn);
                    }
                }
                e[m - 2] = f;
                iter += 1;
            }
            // Convergence: make the singular value positive and order it.
            _ => {
                if s[l - 1] < 0.0 {
                    s[l - 1] = -s[l - 1];
                    if wantv {
                        dscal(p, -1.0, &mut v[at(ldv, 1, l)..], 1);
                    }
                }
                let mut l = l;
                while l != mm && s[l - 1] < s[l] {
                    s.swap(l - 1, l);
                    if wantv && l < p {
                        let (a, b) = two_cols_mut(v, ldv, l, l + 1);
                        dswap(p, a, 1, b, 1);
                    }
                    if wantu && l < n {
                        let (a, b) = two_cols_mut(u, ldu, l, l + 1);
                        dswap(n, a, 1, b, 1);
                    }
                    l += 1;
                }
                iter = 0;
                m -= 1;
            }
        }
    }

    info
}

/// `drot` with the two vectors given in the opposite order to the rotation.
///
/// `dsvdc`'s kase 2 rotates `(u(:,k), u(:,l-1))` where `k > l-1`, so the
/// mutable split yields them the other way round. Rotating `(y, x)` with the
/// same `(c, s)` is not the same operation as rotating `(x, y)`, so the
/// arithmetic is written out rather than reordered.
fn drot_swapped(n: usize, dx: &mut [f64], dy: &mut [f64], c: f64, s: f64) {
    for i in 0..n {
        let dtemp = c * dx[i] + s * dy[i];
        dy[i] = c * dy[i] - s * dx[i];
        dx[i] = dtemp;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// R = Q' A, so reconstructing A from the decomposition must return it.
    #[test]
    fn dqrdc_decomposes_a_small_matrix() {
        // 4x2 column-major.
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0];
        let mut x = a.clone();
        let mut qraux = vec![0.0; 2];
        dqrdc(&mut x, 4, 4, 2, &mut qraux, 0);

        // R's diagonal must be non-zero for a full-rank input.
        assert!(x[at(4, 1, 1)].abs() > 1e-12);
        assert!(x[at(4, 2, 2)].abs() > 1e-12);

        // Q'A = R means ||Q'a_1|| = ||a_1||, and Q'a_1 has one non-zero entry.
        let mut qty = vec![0.0; 4];
        let col1 = &a[0..4];
        dqrsl_qty(&mut x, 4, 4, 2, &qraux, col1, &mut qty);
        let norm: f64 = col1.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!((qty[0].abs() - norm).abs() < 1e-12);
        for v in &qty[1..] {
            assert!(v.abs() < 1e-12, "Q' a_1 should be zero below the first row");
        }
    }

    #[test]
    fn dqrsl_qty_preserves_the_norm() {
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0];
        let mut x = a.clone();
        let mut qraux = vec![0.0; 2];
        dqrdc(&mut x, 4, 4, 2, &mut qraux, 0);

        let y = [2.0, -1.0, 0.5, 3.0];
        let mut qty = vec![0.0; 4];
        dqrsl_qty(&mut x, 4, 4, 2, &qraux, &y, &mut qty);

        let ny: f64 = y.iter().map(|v| v * v).sum();
        let nq: f64 = qty.iter().map(|v| v * v).sum();
        assert!((ny - nq).abs() < 1e-12, "an orthogonal map preserves the norm");
    }

    /// The decomposition must leave the matrix it borrowed unchanged.
    #[test]
    fn dqrsl_qty_restores_the_diagonal() {
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0];
        let mut x = a.clone();
        let mut qraux = vec![0.0; 2];
        dqrdc(&mut x, 4, 4, 2, &mut qraux, 0);

        let before = x.clone();
        let mut qty = vec![0.0; 4];
        dqrsl_qty(&mut x, 4, 4, 2, &qraux, &[1.0, 1.0, 1.0, 1.0], &mut qty);
        assert_eq!(x, before);
    }

    /// Singular values of a diagonal matrix are its entries, in order.
    #[test]
    fn dsvdc_of_a_diagonal_matrix() {
        let mut x = vec![0.0; 9];
        x[at(3, 1, 1)] = 2.0;
        x[at(3, 2, 2)] = 5.0;
        x[at(3, 3, 3)] = 1.0;

        let (mut s, mut e, mut work) = (vec![0.0; 4], vec![0.0; 3], vec![0.0; 3]);
        let mut u = vec![0.0; 9];
        let mut v = vec![0.0; 9];
        let info = dsvdc(&mut x, 3, 3, 3, &mut s, &mut e, &mut u, 3, &mut v, 3, &mut work, 21);

        assert_eq!(info, 0);
        assert!((s[0] - 5.0).abs() < 1e-12);
        assert!((s[1] - 2.0).abs() < 1e-12);
        assert!((s[2] - 1.0).abs() < 1e-12);
    }

    /// The decomposition must reconstruct: A = U S V'.
    #[test]
    fn dsvdc_reconstructs_the_matrix() {
        let a = vec![4.0, 1.0, 0.5, 2.0, 3.0, 1.5, 0.25, 0.75, 2.5];
        let mut x = a.clone();
        let (mut s, mut e, mut work) = (vec![0.0; 4], vec![0.0; 3], vec![0.0; 3]);
        let mut u = vec![0.0; 9];
        let mut v = vec![0.0; 9];
        let info = dsvdc(&mut x, 3, 3, 3, &mut s, &mut e, &mut u, 3, &mut v, 3, &mut work, 21);
        assert_eq!(info, 0);

        for i in 1..=3 {
            for j in 1..=3 {
                let mut acc = 0.0;
                for k in 1..=3 {
                    acc += u[at(3, i, k)] * s[k - 1] * v[at(3, j, k)];
                }
                assert!(
                    (acc - a[at(3, i, j)]).abs() < 1e-10,
                    "A[{i},{j}]: reconstructed {acc}, want {}",
                    a[at(3, i, j)]
                );
            }
        }
    }

    /// Singular values come out sorted, largest first.
    #[test]
    fn dsvdc_orders_the_singular_values() {
        let mut x = vec![1.0, 4.0, 2.0, 7.0, 3.0, 1.0, 9.0, 2.0, 5.0];
        let (mut s, mut e, mut work) = (vec![0.0; 4], vec![0.0; 3], vec![0.0; 3]);
        let mut u = vec![0.0; 9];
        let mut v = vec![0.0; 9];
        let _ = dsvdc(&mut x, 3, 3, 3, &mut s, &mut e, &mut u, 3, &mut v, 3, &mut work, 21);
        assert!(s[0] >= s[1] && s[1] >= s[2]);
        assert!(s[2] >= 0.0);
    }

    /// A rank-deficient matrix must report a zero singular value rather than
    /// failing, since `ehg127` relies on that to detect a singular fit.
    #[test]
    fn dsvdc_reports_rank_deficiency() {
        // Third column is twice the first.
        let mut x = vec![1.0, 2.0, 3.0, 0.0, 1.0, 1.0, 2.0, 4.0, 6.0];
        let (mut s, mut e, mut work) = (vec![0.0; 4], vec![0.0; 3], vec![0.0; 3]);
        let mut u = vec![0.0; 9];
        let mut v = vec![0.0; 9];
        let info = dsvdc(&mut x, 3, 3, 3, &mut s, &mut e, &mut u, 3, &mut v, 3, &mut work, 21);
        assert_eq!(info, 0);
        assert!(s[2] < 1e-14, "expected a zero singular value, got {}", s[2]);
    }
}
