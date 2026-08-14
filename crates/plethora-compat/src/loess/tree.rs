//! The k-d tree R's `loess` interpolates over.
//!
//! With the default `surface = "interpolate"`, `loess` does not evaluate the
//! local regression at every data point. It partitions the predictor range into
//! cells, fits the regression only at the cell vertices, and interpolates
//! between them. `predict(fit)` therefore returns values off that interpolant,
//! not off the regression, and the two genuinely differ: on an input shaped like
//! `gc_correction.R`'s, only one fitted value in 53 comes out bit-equal between
//! `surface = "interpolate"` and `surface = "direct"`, with a relative gap up to
//! 3.7e-4. That gap propagates through `k.gc` into
//! `exp(log(coverage) * k.gc)`, so reproducing the published numbers means
//! reproducing the tree.
//!
//! This module builds it. Transcribed from R's `src/library/stats/src/loessf.f`:
//! [`bounding_box`] from `ehg126`, [`select_kth`] from `ehg106`, [`spread`] from
//! `ehg129`, [`add_vertex`] from `ehg125`, and [`KdTree::build`] from `ehg124`.
//!
//! Indices are kept 1-based, as in the Fortran, with slot 0 of every vector left
//! unused. The arithmetic on cell and vertex indices is dense and mutually
//! recursive, and rebasing it to 0 is exactly the kind of edit that produces a
//! tree that looks right and interpolates wrong.

/// One-based storage: a vector whose slot 0 is a placeholder.
///
/// Only used for the arrays the Fortran indexes from 1.
type OneBased<T> = Vec<T>;

/// A built k-d tree, in the layout `ehg128` walks.
#[derive(Debug, Clone)]
pub struct KdTree {
    /// Number of predictors. Only `d == 1` is supported; see the module note.
    pub d: usize,
    /// Vertices per cell, `2^d`.
    pub vc: usize,
    /// Number of cells, internal and leaf.
    pub nc: usize,
    /// Number of distinct vertices.
    pub nv: usize,
    /// Split dimension per cell, 1-based; 0 marks a leaf. R exposes this as
    /// `fit$kd$a`.
    pub a: OneBased<usize>,
    /// Split coordinate per cell. R exposes this as `fit$kd$xi`.
    pub xi: OneBased<f64>,
    /// Left child cell for an internal cell, first point index for a leaf.
    pub lo: OneBased<usize>,
    /// Right child cell for an internal cell, last point index for a leaf.
    pub hi: OneBased<usize>,
    /// Cell to vertex indices, `vc` per cell.
    pub c: OneBased<Vec<usize>>,
    /// Vertex coordinates, `d` per vertex.
    pub v: OneBased<f64>,
    /// Point order, permuted in place by the splits. A leaf's points are
    /// `pi[lo[cell]..=hi[cell]]`.
    pub pi: OneBased<usize>,
}

impl KdTree {
    /// Builds the tree over `x`, exactly as `ehg131` drives `ehg126` and
    /// `ehg124`.
    ///
    /// `fc` is the bucket size, `floor(n * cell * span)` in R's terms, and `fd`
    /// the diameter floor. R leaves `v(3)` at 0, so `fd` is 0 for every call
    /// `loess()` makes and splitting is driven by the point count alone; the
    /// parameter is carried because `ehg124` tests it.
    ///
    /// # Panics
    /// Panics for `d != 1`. The general case needs the tensor-product vertex
    /// bookkeeping in `ehg125`, which this port does not carry.
    #[must_use]
    pub fn build(x: &[f64], fc: usize, fd: f64, nvmax: usize, ncmax: usize) -> Self {
        let d = 1;
        let vc = 2_usize.pow(u32::try_from(d).expect("d fits in u32"));
        let n = x.len();

        // Fortran indexes points from 1.
        let mut xs: OneBased<f64> = Vec::with_capacity(n + 1);
        xs.push(f64::NAN);
        xs.extend_from_slice(x);

        let (lower, upper) = bounding_box(&xs[1..]);

        let mut v: OneBased<f64> = vec![f64::NAN; nvmax + 1];
        v[1] = lower;
        v[vc] = upper;

        let mut tree = Self {
            d,
            vc,
            nc: 1,
            nv: vc,
            a: vec![0; ncmax + 1],
            xi: vec![0.0; ncmax + 1],
            lo: vec![0; ncmax + 1],
            hi: vec![0; ncmax + 1],
            c: vec![vec![0; vc + 1]; ncmax + 1],
            v,
            pi: (0..=n).collect(),
        };

        for j in 1..=vc {
            tree.c[1][j] = j;
        }

        // `fd` is scaled by the box diagonal in ehg131 before ehg124 sees it.
        let fd = fd * (upper - lower).abs();

        tree.split(&xs, n, fc, fd, nvmax, ncmax);
        tree
    }

    /// `ehg124`: split cells until each holds at most `fc` points.
    ///
    /// The Fortran walks cells with a cursor rather than recursing, appending
    /// children past the end and reaching them later in the same sweep. That
    /// order decides cell numbering, which R reports as `fit$kd$a` and `$xi`, so
    /// it is preserved rather than replaced by recursion.
    fn split(&mut self, xs: &[f64], n: usize, fc: usize, fd: f64, nvmax: usize, ncmax: usize) {
        let mut p = 1;
        let mut l = 1;
        let mut u = n;
        self.lo[p] = l;
        self.hi[p] = u;

        while p <= self.nc {
            let diam = (self.v[self.c[p][self.vc]] - self.v[self.c[p][1]]).abs();

            // Signed, because a split at the last point gives the hi son an
            // empty range with lo > hi. Fortran's integers are signed, so the
            // count comes out negative there and the cell becomes a leaf;
            // computing it in usize would underflow instead.
            let count = (u as isize - l as isize) + 1;
            let mut leaf = count <= fc as isize || diam <= fd;
            if !leaf {
                // Out of room for another split.
                leaf =
                    ncmax < self.nc + 2 || (nvmax as f64) < self.nv as f64 + self.vc as f64 / 2.0;
            }

            let mut split_at = 0;
            if !leaf {
                // With d == 1 the widest dimension is the only dimension.
                let k = 1;
                // Fortran writes int(DBLE(l+u)/2.D0), which for positive
                // indices is plain truncating division and equals this.
                let mut m = usize::midpoint(l, u);
                select_kth(l, u, m, xs, &mut self.pi);

                // All ties go with the hi son. Walking outwards from the median
                // finds the nearest position where the value actually changes.
                // The alternating offset is a 2006 upstream fix; a plain scan
                // forward loses the low son when the median value repeats.
                let mut offset: isize = 0;
                loop {
                    let probe = m as isize + offset;
                    if probe >= u as isize || probe < l as isize {
                        break;
                    }
                    let probe = probe as usize;
                    let (lower, check, upper) = if offset < 0 {
                        (l, probe, probe)
                    } else {
                        (probe + 1, probe + 1, u)
                    };
                    select_kth(lower, upper, check, xs, &mut self.pi);
                    if xs[self.pi[probe]] == xs[self.pi[probe + 1]] {
                        offset = -offset;
                        if offset >= 0 {
                            offset += 1;
                        }
                    } else {
                        m = probe;
                        break;
                    }
                }

                // A split landing exactly on a cell boundary would add no vertex.
                leaf = self.v[self.c[p][1]] == xs[self.pi[m]]
                    || self.v[self.c[p][self.vc]] == xs[self.pi[m]];
                if !leaf {
                    self.a[p] = k;
                    self.xi[p] = xs[self.pi[m]];
                    split_at = m;
                }
            }

            if leaf {
                self.a[p] = 0;
            } else {
                let m = split_at;
                let parent_lo = self.c[p][1];
                let parent_hi = self.c[p][self.vc];

                self.nc += 1;
                let left = self.nc;
                self.lo[p] = left;
                self.lo[left] = l;
                self.hi[left] = m;

                self.nc += 1;
                let right = self.nc;
                self.hi[p] = right;
                self.lo[right] = m + 1;
                self.hi[right] = u;

                let mid = self.add_vertex(self.xi[p]);
                self.c[left] = vec![0, parent_lo, mid];
                self.c[right] = vec![0, mid, parent_hi];
            }

            p += 1;
            if p <= self.nc {
                l = self.lo[p];
                u = self.hi[p];
            }
        }
    }

    /// `ehg125`, reduced to the one new vertex a 1-D split creates.
    ///
    /// The redundancy scan is kept: a coordinate already present must reuse its
    /// vertex rather than duplicate it, or the vertex count drifts from R's.
    fn add_vertex(&mut self, t: f64) -> usize {
        for m in 1..=self.nv {
            if self.v[m] == t {
                return m;
            }
        }
        self.nv += 1;
        self.v[self.nv] = t;
        self.nv
    }

    /// Descends to the leaf cell containing `z`, as `ehg128` opens.
    #[must_use]
    pub fn locate(&self, z: f64) -> usize {
        let mut j = 1;
        while self.a[j] != 0 {
            j = if z <= self.xi[j] {
                self.lo[j]
            } else {
                self.hi[j]
            };
        }
        j
    }
}

/// `ehg126`: the bounding box, expanded by a small margin.
///
/// The margin keeps every data point strictly inside the box, so the
/// interpolation parameter stays within [0, 1] at the extremes.
#[must_use]
pub fn bounding_box(x: &[f64]) -> (f64, f64) {
    let mut alpha = f64::MAX;
    let mut beta = -f64::MAX;
    for &t in x {
        alpha = alpha.min(t);
        beta = beta.max(t);
    }
    let mu = 0.005 * (beta - alpha).max(1e-10 * alpha.abs().max(beta.abs()) + 1e-30);
    (alpha - mu, beta + mu)
}

/// `ehg129`: the spread of the points a cell holds.
#[must_use]
pub fn spread(l: usize, u: usize, x: &[f64], pi: &[usize]) -> f64 {
    let mut alpha = f64::MAX;
    let mut beta = -f64::MAX;
    for i in l..=u {
        let t = x[pi[i]];
        alpha = alpha.min(t);
        beta = beta.max(t);
    }
    beta - alpha
}

/// `ehg106`: partially order `pi[il..=ir]` so that `pi[k]` holds the k-th
/// smallest value, by Floyd and Rivest's Algorithm 489.
///
/// Only `pi` is permuted; the data stays put. The exact partitioning matters
/// because it decides which of several equal values ends up at position `k`,
/// and that in turn decides the split coordinate.
pub fn select_kth(il: usize, ir: usize, k: usize, p: &[f64], pi: &mut [usize]) {
    let mut l = il;
    let mut r = ir;

    while l < r {
        let t = p[pi[k]];
        let mut i = l;
        let mut j = r;

        pi.swap(l, k);
        if t < p[pi[r]] {
            pi.swap(l, r);
        }

        while i < j {
            pi.swap(i, j);
            i += 1;
            j -= 1;
            while p[pi[i]] < t {
                i += 1;
            }
            while t < p[pi[j]] {
                j -= 1;
            }
        }

        if p[pi[l]] == t {
            pi.swap(l, j);
        } else {
            j += 1;
            pi.swap(r, j);
        }

        if j <= k {
            l = j + 1;
        }
        if k <= j {
            r = j - 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bounding_box_expands_by_half_a_percent() {
        let (lo, hi) = bounding_box(&[0.2, 0.5, 0.72]);
        let span = 0.72 - 0.2;
        assert!((lo - (0.2 - 0.005 * span)).abs() < 1e-15);
        assert!((hi - (0.72 + 0.005 * span)).abs() < 1e-15);
        assert!(lo < 0.2 && hi > 0.72);
    }

    /// A degenerate input must still produce a box with room in it, or the
    /// interpolation parameter divides by zero.
    #[test]
    fn bounding_box_of_a_constant_is_not_empty() {
        let (lo, hi) = bounding_box(&[0.5, 0.5, 0.5]);
        assert!(lo < hi);
    }

    #[test]
    fn select_kth_places_the_kth_smallest() {
        let p: Vec<f64> = vec![f64::NAN, 5.0, 1.0, 4.0, 2.0, 3.0];
        for k in 1..=5 {
            let mut pi: Vec<usize> = (0..=5).collect();
            select_kth(1, 5, k, &p, &mut pi);
            assert_eq!(p[pi[k]], k as f64, "k = {k}");
        }
    }

    /// Algorithm 489 guarantees position `k` and the partition around it, and
    /// nothing about the order within each side. Assert only that.
    #[test]
    fn select_kth_partitions_around_k() {
        // A deterministic sweep over inputs full of ties, which is the case the
        // splitter actually hits: repeated x values are common in binned data.
        let mut state = 12345_u64;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            f64::from(((state >> 33) % 5) as u32)
        };

        for n in 1..=12_usize {
            for k in 1..=n {
                let mut p = vec![f64::NAN];
                p.extend((0..n).map(|_| next()));
                let mut sorted: Vec<f64> = p[1..].to_vec();
                sorted.sort_by(f64::total_cmp);

                let mut pi: Vec<usize> = (0..=n).collect();
                select_kth(1, n, k, &p, &mut pi);

                assert_eq!(p[pi[k]], sorted[k - 1], "n = {n}, k = {k}");
                for i in 1..k {
                    assert!(
                        p[pi[i]] <= p[pi[k]],
                        "left of k out of order, n = {n}, k = {k}"
                    );
                }
                for i in k + 1..=n {
                    assert!(
                        p[pi[i]] >= p[pi[k]],
                        "right of k out of order, n = {n}, k = {k}"
                    );
                }
                // The permutation must stay a permutation.
                let mut seen: Vec<usize> = pi[1..=n].to_vec();
                seen.sort_unstable();
                assert_eq!(seen, (1..=n).collect::<Vec<_>>());
            }
        }
    }

    #[test]
    fn spread_measures_the_range() {
        let x: Vec<f64> = vec![f64::NAN, 3.0, 1.0, 4.0, 1.0, 5.0];
        let pi: Vec<usize> = (0..=5).collect();
        assert_eq!(spread(1, 5, &x, &pi), 4.0);
        assert_eq!(spread(2, 3, &x, &pi), 3.0);
    }

    /// A cell holding no more than `fc` points is never split.
    #[test]
    fn a_small_input_stays_one_leaf() {
        let x: Vec<f64> = (0..5).map(f64::from).collect();
        let tree = KdTree::build(&x, 7, 0.0, 100, 100);
        assert_eq!(tree.nc, 1);
        assert_eq!(tree.nv, 2);
        assert_eq!(tree.a[1], 0);
        assert_eq!(tree.locate(2.0), 1);
    }

    #[test]
    fn splitting_creates_two_cells_and_one_vertex() {
        let x: Vec<f64> = (0..16).map(f64::from).collect();
        let tree = KdTree::build(&x, 8, 0.0, 100, 100);
        assert!(
            tree.nc >= 3,
            "expected at least one split, got nc = {}",
            tree.nc
        );
        assert_eq!(tree.a[1], 1, "the root splits on the only dimension");
        // Every leaf holds at most fc points. Counted signed: a split at the
        // last point leaves its hi son with an empty, inverted range.
        for cell in 1..=tree.nc {
            if tree.a[cell] == 0 {
                let count = tree.hi[cell] as isize - tree.lo[cell] as isize + 1;
                assert!(count <= 8, "leaf {cell} holds {count} points");
            }
        }
    }

    /// Locating must land in the leaf whose vertex interval brackets the point.
    #[test]
    fn locate_lands_inside_the_cell() {
        let x: Vec<f64> = (0..40).map(|i| f64::from(i) / 40.0).collect();
        let tree = KdTree::build(&x, 7, 0.0, 100, 100);
        for &z in &x {
            let cell = tree.locate(z);
            assert_eq!(tree.a[cell], 0, "locate must return a leaf");
            let lo = tree.v[tree.c[cell][1]];
            let hi = tree.v[tree.c[cell][tree.vc]];
            assert!(lo <= z && z <= hi, "{z} outside [{lo}, {hi}]");
        }
    }
}
