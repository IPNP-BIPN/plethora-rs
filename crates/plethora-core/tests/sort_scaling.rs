//! The parallel sort, checked and measured.
//!
//! `bed::sort` sorts each run across every core. That is only safe because
//! `cmp_k1_k2n` falls back to comparing whole lines, so two lines compare
//! `Equal` only when they are byte-identical and a reordering among equals
//! cannot show in the output. The first test holds that to account.
//!
//! The second is a measurement rather than a test, so it is ignored by default:
//!
//! ```text
//! cargo test --release -p plethora-core --test sort_scaling -- --ignored --nocapture
//! ```
//!
//! On sixteen cores it reports about 6.9x, and the shape is what matters more
//! than the number: the sort is a fraction of a second on the test corpora here
//! and minutes on a whole-genome sample, where it is the stage that dominates.
use std::time::Instant;

use plethora_compat::gnusort::cmp_k1_k2n;
use rayon::slice::ParallelSliceMut as _;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 11
    }
}

fn corpus(n: usize) -> Vec<String> {
    let mut rng = Lcg(0x1234_5678);
    let chroms = ["chr1", "chr2", "chr10", "chrX"];
    (0..n)
        .map(|_| {
            let c = chroms[(rng.next() % 4) as usize];
            let s = rng.next() % 250_000_000;
            format!("{c}\t{s}\t{}\tf{}\t60\t+", s + 300, rng.next() % 1_000_000)
        })
        .collect()
}

/// The two orders agree, which is the whole safety argument.
#[test]
fn the_parallel_order_is_the_sequential_order() {
    let lines = corpus(200_000);
    let mut sequential = lines.clone();
    sequential.sort_unstable_by(|x, y| cmp_k1_k2n(x.as_bytes(), y.as_bytes()));
    let mut parallel = lines;
    parallel.par_sort_unstable_by(|x, y| cmp_k1_k2n(x.as_bytes(), y.as_bytes()));
    assert_eq!(sequential, parallel);
}

#[test]
#[ignore = "a measurement, not a test"]
fn sort_scaling() {
    for n in [1_000_000usize, 10_000_000] {
        let lines = corpus(n);
        let mut a = lines.clone();
        let t = Instant::now();
        a.sort_unstable_by(|x, y| cmp_k1_k2n(x.as_bytes(), y.as_bytes()));
        let seq = t.elapsed();

        let mut b = lines;
        let t = Instant::now();
        b.par_sort_unstable_by(|x, y| cmp_k1_k2n(x.as_bytes(), y.as_bytes()));
        let par = t.elapsed();

        assert_eq!(a, b, "the two orders must agree");
        println!(
            "  {n:>9} lines: sequential {seq:>12.2?}  parallel {par:>12.2?}  {:.1}x",
            seq.as_secs_f64() / par.as_secs_f64()
        );
    }
}
