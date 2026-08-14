//! Differential tests against the real `gsort` and `samtools`.
//!
//! The unit tests inside `gnusort` and `strnum` pin behaviour that was read off
//! the source. These generate bulk data instead and let the actual tools decide,
//! which is what catches a rule that was read correctly but understood wrongly.
//!
//! Both tests skip when their tool is missing, so the suite still runs on a
//! machine without the bioinformatics stack installed. They are not the only
//! coverage of these comparators, so skipping loses breadth, not correctness.

use std::cmp::Ordering;
// Both write traits are in scope: `write!` picks by receiver, so a String
// resolves through fmt and a child's stdin through io.
use std::fmt::Write as _;
use std::io::Write as _;
use std::process::{Command, Stdio};

use plethora_compat::gnusort::cmp_k1_k2n;
use plethora_compat::strnum::cmp_by_qname;

/// A tiny deterministic generator, so the corpus is identical on every run and
/// on every machine without pulling in a dependency.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        // Numerical Recipes' constants for a 64-bit LCG.
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 11
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }

    fn pick<'a>(&mut self, options: &[&'a str]) -> &'a str {
        options[self.below(options.len() as u64) as usize]
    }
}

fn tool_exists(tool: &str) -> bool {
    Command::new(tool)
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

/// A BED-shaped corpus with heavy collisions on both keys, so the last-resort
/// whole-line comparison is what orders most of it.
fn bed_corpus() -> Vec<String> {
    let mut rng = Lcg(0x5eed_1234);
    let chroms = [
        "chr1",
        "chr2",
        "chr10",
        "chr20",
        "chrX",
        "chrY",
        "chrM",
        "chr1_KI270706v1_random",
        "A",
        "a",
    ];
    let names = [
        "read", "read1", "read10", "read2", "a.b", "a-b", "a_b", ".", "zzz",
    ];

    (0..4000)
        .map(|_| {
            // A small coordinate space forces frequent ties on key 2.
            let start = rng.below(50);
            let end = start + 1 + rng.below(100);
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}",
                rng.pick(&chroms),
                start,
                end,
                rng.pick(&names),
                rng.below(256),
                if rng.below(2) == 0 { "+" } else { "-" },
            )
        })
        .collect()
}

#[test]
fn gnusort_matches_gnu_sort() {
    if !tool_exists("gsort") {
        eprintln!("skipping: gsort (GNU coreutils) not installed");
        return;
    }

    let lines = bed_corpus();
    let input = format!("{}\n", lines.join("\n"));

    let mut child = Command::new("gsort")
        .args(["-k1,1", "-k2,2n"])
        .env("LC_ALL", "C")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("spawn gsort");
    child
        .stdin
        .take()
        .expect("gsort stdin")
        .write_all(input.as_bytes())
        .expect("write to gsort");
    let out = child.wait_with_output().expect("gsort output");
    assert!(out.status.success(), "gsort failed");
    let expected: Vec<&str> = std::str::from_utf8(&out.stdout)
        .expect("gsort emitted utf-8")
        .lines()
        .collect();

    let mut got = lines.clone();
    got.sort_by(|a, b| cmp_k1_k2n(a.as_bytes(), b.as_bytes()));

    assert_eq!(got.len(), expected.len(), "line count changed");
    for (i, (g, e)) in got.iter().zip(&expected).enumerate() {
        assert_eq!(
            g, e,
            "line {i} differs from GNU sort\n  ours: {g:?}\n  gsort: {e:?}"
        );
    }
}

/// QNAMEs chosen to exercise every branch of `strnum_cmp`: digit runs of
/// different lengths, leading zeros, digits at the end of a name, numbers wider
/// than `u64`, and pure-ASCII neighbours.
fn qname_corpus() -> Vec<(String, u16)> {
    let mut rng = Lcg(0xc0ffee);
    let stems = ["read", "a", "A", "x", "sim", "HWI-ST745", "", "z9"];
    let flags = [
        0x40_u16,
        0x80,
        0x40 | 0x100,
        0x40 | 0x800,
        0x80 | 0x100,
        0x80 | 0x800,
    ];

    let mut names = vec![
        "9999999999999999999".to_string(),
        "10000000000000000000".to_string(),
        "a1".to_string(),
        "a01".to_string(),
        "a001".to_string(),
        "a0".to_string(),
        "a1b2".to_string(),
        "a1b10".to_string(),
    ];
    for _ in 0..600 {
        let stem = rng.pick(&stems);
        let zeros = "0".repeat(rng.below(3) as usize);
        let n = rng.below(120);
        let tail = if rng.below(3) == 0 { "b" } else { "" };
        names.push(format!("{stem}{zeros}{n}{tail}"));
    }

    names
        .into_iter()
        .flat_map(|n| {
            // Every name appears with two records so the flag key is exercised.
            let f1 = flags[rng.below(flags.len() as u64) as usize];
            let f2 = flags[rng.below(flags.len() as u64) as usize];
            [(n.clone(), f1), (n, f2)]
        })
        .collect()
}

#[test]
fn strnum_matches_samtools_sort_n() {
    if !tool_exists("samtools") {
        eprintln!("skipping: samtools not installed");
        return;
    }

    let records = qname_corpus();

    let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000000\n");
    for (i, (name, flag)) in records.iter().enumerate() {
        // QNAME "*" is reserved, and samtools rejects an empty one.
        let qname = if name.is_empty() { "unnamed" } else { name };
        // POS encodes the input index so a stability failure is legible.
        writeln!(
            sam,
            "{qname}\t{flag}\tchr1\t{}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
            i + 1
        )
        .expect("format into a String cannot fail");
    }

    let mut child = Command::new("samtools")
        .args(["sort", "-n", "-@", "1", "-O", "sam"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("spawn samtools");
    child
        .stdin
        .take()
        .expect("samtools stdin")
        .write_all(sam.as_bytes())
        .expect("write to samtools");
    let out = child.wait_with_output().expect("samtools output");
    assert!(out.status.success(), "samtools sort failed");

    // Recover the input index from POS, giving the permutation samtools chose.
    let expected: Vec<usize> = std::str::from_utf8(&out.stdout)
        .expect("samtools emitted utf-8")
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| {
            let pos: usize = l
                .split('\t')
                .nth(3)
                .expect("POS column")
                .parse()
                .expect("POS is numeric");
            pos - 1
        })
        .collect();

    let mut got: Vec<usize> = (0..records.len()).collect();
    // Stable, because samtools' merge preserves input order on a full tie.
    got.sort_by(|&i, &j| {
        let (ref ni, fi) = records[i];
        let (ref nj, fj) = records[j];
        let ni = if ni.is_empty() { "unnamed" } else { ni };
        let nj = if nj.is_empty() { "unnamed" } else { nj };
        cmp_by_qname(ni.as_bytes(), fi, nj.as_bytes(), fj)
    });

    assert_eq!(got.len(), expected.len(), "record count changed");
    for (rank, (&g, &e)) in got.iter().zip(&expected).enumerate() {
        assert_eq!(
            g, e,
            "rank {rank} differs from samtools\n  ours:     {:?}\n  samtools: {:?}",
            records[g], records[e]
        );
    }
}

/// The comparators must be total orders, or a sort built on them is undefined.
#[test]
fn comparators_are_consistent() {
    let lines = bed_corpus();
    for pair in lines.windows(2) {
        let (a, b) = (pair[0].as_bytes(), pair[1].as_bytes());
        assert_eq!(
            cmp_k1_k2n(a, b),
            cmp_k1_k2n(b, a).reverse(),
            "cmp_k1_k2n is not antisymmetric for {a:?} and {b:?}"
        );
        assert_eq!(cmp_k1_k2n(a, a), Ordering::Equal);
    }
}
