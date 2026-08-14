//! The whole pipeline, from an alignment to copy number.
//!
//! Every other test in this suite checks one stage against the tool it
//! replaces. This one checks that the stages fit together, which is a different
//! failure: a module can be right on its own and still be handed a file the
//! previous one wrote under another name, or compressed when it expected plain.
//! That is exactly the bug this test was written after finding.
//!
//! The assertion at the end is the invariant the whole calibration rests on:
//! the conserved regions are defined as diploid, so their median corrected
//! coverage must come out at 2.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::io::Write as _;
use std::path::Path;
use std::process::{Command, Stdio};

use plethora_core::coverage::{self, Pairing};
use plethora_core::gc;
use plethora_core::io::Compress;

/// A deterministic generator, so the corpus is identical on every run.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 11
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn have(tool: &str) -> bool {
    Command::new(tool)
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

/// Forty conserved domains and ten of interest, evenly spaced.
fn domains() -> String {
    let mut bed = String::new();
    let mut position = 1000;
    for i in 0..40 {
        writeln!(
            bed,
            "chr1\t{position}\t{}\tbaseline_{i:03}\t255\t+",
            position + 1000
        )
        .unwrap();
        position += 2000;
    }
    for i in 0..10 {
        writeln!(
            bed,
            "chr1\t{position}\t{}\tNBPF1_CON1_{i}\t255\t+",
            position + 1000
        )
        .unwrap();
        position += 2000;
    }
    bed
}

/// A GC table for those domains, all in the fitted window so every one reaches
/// the model.
fn gc_table() -> HashMap<String, f64> {
    let mut rng = Lcg(0x6c00_1100);
    let mut table = HashMap::new();
    for i in 0..40 {
        // Spread across the window, in the three-decimal form a real table has.
        let gc = 0.30 + (rng.below(30) as f64) / 100.0;
        table.insert(format!("baseline_{i:03}"), gc);
    }
    for i in 0..10 {
        table.insert(format!("NBPF1_CON1_{i}"), 0.40);
    }
    table
}

/// A SAM of proper pairs landing inside the domains, with the domains of
/// interest given several times the coverage of the conserved ones.
fn alignment_sam() -> String {
    let mut rng = Lcg(0xa119_0000);
    let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:200000\n");
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
    let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";

    let emit = |sam: &mut String, name: &str, start: i64| {
        let mate = start + 200;
        writeln!(
            sam,
            "{name}\t99\tchr1\t{start}\t60\t50M\t=\t{mate}\t250\t{seq}\t{qual}"
        )
        .unwrap();
        writeln!(
            sam,
            "{name}\t147\tchr1\t{mate}\t60\t50M\t=\t{start}\t-250\t{seq}\t{qual}"
        )
        .unwrap();
    };

    let mut fragment = 0;
    // Conserved domains: a steady depth, which sets the haploid unit.
    for i in 0..40 {
        let domain_start = 1000 + i * 2000;
        for _ in 0..30 {
            let start = domain_start + rng.below(700) as i64 + 1;
            emit(&mut sam, &format!("f{fragment:06}"), start);
            fragment += 1;
        }
    }
    // Domains of interest: three times the depth, so six copies.
    for i in 0..10 {
        let domain_start = 81000 + i * 2000;
        for _ in 0..90 {
            let start = domain_start + rng.below(700) as i64 + 1;
            emit(&mut sam, &format!("f{fragment:06}"), start);
            fragment += 1;
        }
    }
    sam
}

/// R's median, which is what the correction divides by. R averages the two
/// middle values as `(a + b) / 2`, not as a midpoint, so this does too.
#[allow(clippy::manual_midpoint)]
fn median(values: &mut [f64]) -> f64 {
    values.sort_by(f64::total_cmp);
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2]
    } else {
        (values[n / 2 - 1] + values[n / 2]) / 2.0
    }
}

#[test]
fn an_alignment_becomes_copy_number() {
    if !have("samtools") {
        eprintln!("skipping: samtools not installed to build the BAM");
        return;
    }

    let dir = tempfile::tempdir().expect("tempdir");
    let d = dir.path();

    let reference = d.join("domains.bed");
    std::fs::write(&reference, domains()).expect("write domains");

    let bam = d.join("sample.bam");
    let mut child = Command::new("samtools")
        .args(["view", "-b", "-o"])
        .arg(&bam)
        .stdin(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn samtools");
    child
        .stdin
        .take()
        .expect("stdin")
        .write_all(alignment_sam().as_bytes())
        .expect("write sam");
    let out = child.wait_with_output().expect("samtools status");
    assert!(
        out.status.success(),
        "samtools view failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // Both ways round, because compression is where the stages meet.
    for compress in [Compress::No, Compress::Yes] {
        let prefix = d.join(format!("out_{compress:?}"));
        let outputs = coverage::make_bed(&bam, &reference, Pairing::Paired, &prefix, d, compress)
            .unwrap_or_else(|e| panic!("make_bed with {compress:?}: {e}"));

        assert!(
            outputs.read_depth.exists(),
            "no read depth at {}",
            outputs.read_depth.display()
        );
        assert_eq!(
            plethora_core::io::is_gzip_name(&outputs.read_depth),
            compress == Compress::Yes,
            "the output name does not match the compression setting"
        );

        let depth = read_pairs(&outputs.read_depth);
        assert_eq!(depth.len(), 50, "every domain should have a row");

        let rows = gc::correction::correct(&depth, &gc_table()).expect("correct");
        assert_eq!(rows.len(), 50);

        // The conserved regions are defined as diploid.
        let mut conserved: Vec<f64> = rows
            .iter()
            .filter(|r| r.domain.starts_with("baseline"))
            .map(|r| r.corrected_coverage)
            .collect();
        let m = median(&mut conserved);
        assert!(
            (m - 2.0).abs() < 1e-9,
            "the conserved regions should sit at two copies, got {m} with {compress:?}"
        );

        // And the domains of interest were given three times the depth.
        let mut interest: Vec<f64> = rows
            .iter()
            .filter(|r| r.domain.starts_with("NBPF1"))
            .map(|r| r.corrected_coverage)
            .collect();
        let copies = median(&mut interest);
        assert!(
            (4.0..8.0).contains(&copies),
            "three times the conserved depth should read as about six copies, got {copies}"
        );
    }
}

/// Reads a two-column table, decompressing if needed.
fn read_pairs(path: &Path) -> Vec<(String, f64)> {
    use std::io::BufRead as _;
    plethora_core::io::open(path)
        .expect("open")
        .lines()
        .map_while(Result::ok)
        .filter_map(|l| {
            let (name, value) = l.split_once('\t')?;
            Some((name.to_string(), value.trim().parse().ok()?))
        })
        .collect()
}
