//! Differential test of the `bamtobed` port against the installed bedtools.
//!
//! The unit tests inside `bamtobed` pin rules read off bedtools' source. This
//! one generates a corpus and lets bedtools itself decide, which is what
//! catches a rule read correctly but understood wrongly.
//!
//! Skips when samtools or bedtools is missing, so the suite still runs without
//! the bioinformatics stack installed.

// Both write traits are in scope: `write!` picks by receiver, so a String
// resolves through fmt and a child's stdin through io.
use std::fmt::Write as _;
use std::io::Write as _;
use std::process::{Command, Stdio};

use plethora_core::bam::bamtobed::{Aln, BedpeIter, bed};
use plethora_core::bam::reader::read_bam;

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

fn tool_exists(tool: &str) -> bool {
    Command::new(tool)
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

/// A SAM covering the cases that decide `PrintBedPE`'s branches: both mates
/// mapped, one unmapped, both unmapped, mates on different chromosomes,
/// read 1 downstream of read 2, differing MAPQ, and CIGARs with soft clips,
/// deletions, insertions, hard clips and skips.
fn corpus_sam() -> String {
    let mut rng = Lcg(0x0bed_0001);
    // Names chosen so string ordering differs from reference-id ordering.
    let chroms = ["chr1", "chr2", "chr10", "chrX"];
    let cigars = [
        "50M",
        "10S40M",
        "40M10S",
        "25M100N25M",
        "20M5D30M",
        "20M5I25M",
        "5H50M",
    ];

    let mut sam = String::from("@HD\tVN:1.6\tSO:queryname\n");
    for c in chroms {
        writeln!(sam, "@SQ\tSN:{c}\tLN:2000000").expect("format into a String cannot fail");
    }

    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
    let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";

    // Names are zero-padded so the file is already grouped when sorted by name.
    for i in 0..400 {
        let name = format!("q{i:04}");
        let kind = rng.below(10);

        let c1 = chroms[rng.below(chroms.len() as u64) as usize];
        let c2 = if kind == 6 {
            chroms[rng.below(chroms.len() as u64) as usize]
        } else {
            c1
        };
        let p1 = 1 + rng.below(500_000) as i64;
        let p2 = 1 + rng.below(500_000) as i64;
        let q1 = rng.below(61) as u8;
        let q2 = rng.below(61) as u8;
        let cig1 = cigars[rng.below(cigars.len() as u64) as usize];
        let cig2 = cigars[rng.below(cigars.len() as u64) as usize];

        match kind {
            // Read 1 unmapped.
            7 => {
                writeln!(sam, "{name}\t77\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}")
                    .expect("format into a String cannot fail");
                writeln!(
                    sam,
                    "{name}\t137\t{c2}\t{p2}\t{q2}\t{cig2}\t*\t0\t0\t{seq}\t{qual}"
                )
                .expect("format into a String cannot fail");
            }
            // Read 2 unmapped.
            8 => {
                writeln!(
                    sam,
                    "{name}\t73\t{c1}\t{p1}\t{q1}\t{cig1}\t*\t0\t0\t{seq}\t{qual}"
                )
                .expect("format into a String cannot fail");
                writeln!(sam, "{name}\t133\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}")
                    .expect("format into a String cannot fail");
            }
            // Both unmapped.
            9 => {
                writeln!(sam, "{name}\t77\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}")
                    .expect("format into a String cannot fail");
                writeln!(sam, "{name}\t141\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}")
                    .expect("format into a String cannot fail");
            }
            // Both mapped, with the strands varying.
            _ => {
                let f1 = if rng.below(2) == 0 { 99 } else { 83 };
                let f2 = if f1 == 99 { 147 } else { 163 };
                writeln!(
                    sam,
                    "{name}\t{f1}\t{c1}\t{p1}\t{q1}\t{cig1}\t=\t{p2}\t0\t{seq}\t{qual}"
                )
                .expect("format into a String cannot fail");
                writeln!(
                    sam,
                    "{name}\t{f2}\t{c2}\t{p2}\t{q2}\t{cig2}\t=\t{p1}\t0\t{seq}\t{qual}"
                )
                .expect("format into a String cannot fail");
            }
        }
    }
    sam
}

/// Writes the corpus through samtools into a BAM, returning its path.
fn build_bam(dir: &std::path::Path) -> std::path::PathBuf {
    let bam = dir.join("corpus.bam");
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
        .expect("samtools stdin")
        .write_all(corpus_sam().as_bytes())
        .expect("write to samtools");
    let out = child.wait_with_output().expect("samtools status");
    assert!(
        out.status.success(),
        "samtools view failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    bam
}

fn bedtools(args: &[&str], bam: &std::path::Path) -> Vec<String> {
    let out = Command::new("bedtools")
        .arg("bamtobed")
        .args(args)
        .arg("-i")
        .arg(bam)
        .stderr(Stdio::null())
        .output()
        .expect("run bedtools");
    assert!(out.status.success(), "bedtools bamtobed failed");
    String::from_utf8(out.stdout)
        .expect("bedtools emitted utf-8")
        .lines()
        .map(String::from)
        .collect()
}

fn skip_unless_tools() -> bool {
    if !tool_exists("samtools") || !tool_exists("bedtools") {
        eprintln!("skipping: samtools or bedtools not installed");
        return true;
    }
    false
}

#[test]
fn bedpe_matches_bedtools() {
    if skip_unless_tools() {
        return;
    }
    let dir = tempfile::tempdir().expect("tempdir");
    let bam = build_bam(dir.path());

    // Upstream passes -split as well; it is ignored for -bedpe, and this
    // asserts that rather than assuming it.
    let expected = bedtools(&["-bedpe"], &bam);
    let expected_split = bedtools(&["-split", "-bedpe"], &bam);
    assert_eq!(expected, expected_split, "-split changed the -bedpe output");

    let records: Vec<Aln> = read_bam(&bam).expect("read bam");
    let got: Vec<String> = BedpeIter::new(records.into_iter())
        .map(|r| r.to_string())
        .collect();

    assert_eq!(
        got.len(),
        expected.len(),
        "line count differs from bedtools"
    );
    for (i, (g, e)) in got.iter().zip(&expected).enumerate() {
        assert_eq!(g, e, "line {i} differs from bedtools");
    }
}

#[test]
fn bed_matches_bedtools() {
    if skip_unless_tools() {
        return;
    }
    let dir = tempfile::tempdir().expect("tempdir");
    let bam = build_bam(dir.path());

    let expected = bedtools(&[], &bam);

    let records: Vec<Aln> = read_bam(&bam).expect("read bam");
    let got: Vec<String> = records
        .iter()
        .filter_map(bed)
        .map(|r| r.to_string())
        .collect();

    assert_eq!(
        got.len(),
        expected.len(),
        "line count differs from bedtools"
    );
    for (i, (g, e)) in got.iter().zip(&expected).enumerate() {
        assert_eq!(g, e, "line {i} differs from bedtools");
    }
}

/// The name sort, against samtools itself.
///
/// The corpus is deliberately shuffled before sorting, so the test exercises
/// the ordering rather than confirming that an already-ordered file stays put.
#[test]
fn namesort_matches_samtools() {
    if skip_unless_tools() {
        return;
    }
    let dir = tempfile::tempdir().expect("tempdir");

    // Shuffle the corpus by writing records in a scrambled order.
    let sam = corpus_sam();
    let (header, body): (Vec<&str>, Vec<&str>) = sam.lines().partition(|l| l.starts_with('@'));
    let mut rng = Lcg(0x5caff01d);
    let mut body: Vec<&str> = body;
    for i in (1..body.len()).rev() {
        body.swap(i, rng.below(i as u64 + 1) as usize);
    }
    let shuffled = format!("{}\n{}\n", header.join("\n"), body.join("\n"));

    let bam = dir.path().join("shuffled.bam");
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
        .write_all(shuffled.as_bytes())
        .expect("write");
    let out = child.wait_with_output().expect("status");
    assert!(
        out.status.success(),
        "samtools view failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // What samtools makes of it.
    let sorted = Command::new("samtools")
        .args(["sort", "-n", "-@", "1", "-O", "sam"])
        .arg(&bam)
        .stderr(Stdio::null())
        .output()
        .expect("run samtools sort");
    assert!(sorted.status.success(), "samtools sort failed");
    let expected: Vec<(String, u16, i64)> = String::from_utf8_lossy(&sorted.stdout)
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (
                f[0].to_string(),
                f[1].parse().expect("flag"),
                f[3].parse().expect("pos"),
            )
        })
        .collect();

    // What we make of it.
    let records = read_bam(&bam).expect("read bam");
    let ours = plethora_core::bam::namesort::sort_by_name(records, 64, dir.path()).expect("sort");
    let got: Vec<(String, u16, i64)> = ours
        .iter()
        .map(|a| {
            (
                String::from_utf8_lossy(&a.name).into_owned(),
                a.flags,
                // BED start is zero-based; SAM POS is one-based, and an
                // unmapped record has POS 0.
                if a.chrom.is_some() { a.start + 1 } else { 0 },
            )
        })
        .collect();

    assert_eq!(got.len(), expected.len(), "record count differs");
    for (i, (g, e)) in got.iter().zip(&expected).enumerate() {
        assert_eq!(g, e, "record {i} out of order relative to samtools");
    }
}
