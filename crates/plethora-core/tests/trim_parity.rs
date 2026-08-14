//! What the move from cutadapt to Trim Galore actually changes.
//!
//! This stage is the one place the port diverges from upstream on purpose, so
//! there is nothing to assert parity against. What there is to do is measure
//! the divergence rather than assert it away.
//!
//! Two comparisons, on upstream's own test reads:
//!
//! - With [`Adapters::None`], which reproduces upstream's dummy `XXX` adapter,
//!   the output should track cutadapt closely. Any gap there is a real
//!   difference in quality trimming or length filtering, not a design choice.
//! - With [`Adapters::Detect`], the difference from cutadapt is the cost of
//!   actually removing adapters, and it is reported rather than hidden.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use plethora_core::trim::{Adapters, trim_pair};

fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

fn have(program: &str) -> bool {
    Command::new(program)
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

/// Reads a gzipped FASTQ into (id, sequence, quality) triples.
fn read_fastq(path: &Path) -> Vec<(String, String, String)> {
    let file = std::fs::File::open(path).expect("open fastq");
    let mut text = String::new();
    {
        use std::io::Read as _;
        let mut decoder = flate2::read::MultiGzDecoder::new(file);
        decoder.read_to_string(&mut text).expect("decode fastq");
    }
    text.lines()
        .collect::<Vec<_>>()
        .chunks(4)
        .filter(|c| c.len() == 4)
        .map(|c| (c[0].to_string(), c[1].to_string(), c[3].to_string()))
        .collect()
}

/// Copies the vendored reads into a scratch directory under upstream's names.
fn stage(dir: &Path) -> (PathBuf, PathBuf) {
    let r1 = dir.join("test_1.fastq.gz");
    let r2 = dir.join("test_2.fastq.gz");
    std::fs::copy(data_dir().join("test_1.fastq.gz"), &r1).expect("copy r1");
    std::fs::copy(data_dir().join("test_2.fastq.gz"), &r2).expect("copy r2");
    (r1, r2)
}

/// Runs upstream's exact cutadapt line.
fn run_cutadapt(dir: &Path, r1: &Path, r2: &Path) -> (PathBuf, PathBuf) {
    let o1 = dir.join("cutadapt_1.fastq.gz");
    let o2 = dir.join("cutadapt_2.fastq.gz");
    let out = Command::new("cutadapt")
        .args([
            "-a",
            "XXX",
            "-A",
            "XXX",
            "-q",
            "10",
            "--minimum-length",
            "80",
            "--trim-n",
            "-o",
        ])
        .arg(&o1)
        .arg("-p")
        .arg(&o2)
        .arg(r1)
        .arg(r2)
        .output()
        .expect("run cutadapt");
    assert!(
        out.status.success(),
        "cutadapt failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    (o1, o2)
}

/// How many reads differ, and by how many bases in total.
fn compare(ours: &[(String, String, String)], theirs: &[(String, String, String)]) -> (usize, i64) {
    assert_eq!(ours.len(), theirs.len(), "read count differs");
    let mut differing = 0;
    let mut base_delta = 0_i64;
    for (a, b) in ours.iter().zip(theirs) {
        if a.1 != b.1 {
            differing += 1;
        }
        base_delta += a.1.len() as i64 - b.1.len() as i64;
    }
    (differing, base_delta)
}

/// The dummy-adapter path against cutadapt: this is where the two should agree.
#[test]
fn without_adapters_it_tracks_cutadapt() {
    if !have("cutadapt") {
        eprintln!("skipping: cutadapt not installed");
        return;
    }

    let dir = tempfile::tempdir().expect("tempdir");
    let (r1, r2) = stage(dir.path());
    let (c1, c2) = run_cutadapt(dir.path(), &r1, &r2);

    let summary = trim_pair(&r1, &r2, &Adapters::None).expect("trim");
    assert_eq!(summary.adapter, "none");

    let theirs_r1 = read_fastq(&c1);
    let theirs_r2 = read_fastq(&c2);
    let ours_r1 = read_fastq(&summary.output_r1);
    let ours_r2 = read_fastq(&summary.output_r2);

    println!(
        "no adapter: cutadapt kept {} pairs, Trim Galore kept {}",
        theirs_r1.len(),
        ours_r1.len()
    );
    assert_eq!(
        ours_r1.len(),
        theirs_r1.len(),
        "the length filter disagrees with cutadapt"
    );

    let (d1, b1) = compare(&ours_r1, &theirs_r1);
    let (d2, b2) = compare(&ours_r2, &theirs_r2);
    println!(
        "no adapter: R1 {d1} reads differ ({b1:+} bases), R2 {d2} reads differ ({b2:+} bases)"
    );

    // Quality trimming and N trimming are the only operations left, and both
    // are well defined, so any disagreement here is a real one worth seeing.
    assert_eq!(d1 + d2, 0, "quality or N trimming disagrees with cutadapt");
}

/// The divergence the project chose: adapters actually removed.
#[test]
fn detection_reports_what_it_changes() {
    if !have("cutadapt") {
        eprintln!("skipping: cutadapt not installed");
        return;
    }

    let dir = tempfile::tempdir().expect("tempdir");
    let (r1, r2) = stage(dir.path());
    let (c1, _c2) = run_cutadapt(dir.path(), &r1, &r2);

    let summary = trim_pair(&r1, &r2, &Adapters::Detect).expect("trim");
    println!("detected adapter: {}", summary.adapter);
    if let Some(message) = &summary.detection {
        println!("detector says: {}", message.lines().next().unwrap_or(""));
    }

    let theirs = read_fastq(&c1);
    let ours = read_fastq(&summary.output_r1);
    println!(
        "detect: cutadapt kept {} pairs, Trim Galore kept {}",
        theirs.len(),
        ours.len()
    );

    // Not an equality: this is the measurement, not a contract. On these
    // simulated reads there is no real adapter to find, so the difference
    // should be small; on a real library it would not be.
    if ours.len() == theirs.len() {
        let (differing, delta) = compare(&ours, &theirs);
        println!(
            "detect: {differing} of {} reads differ ({delta:+} bases)",
            ours.len()
        );
    }

    assert!(summary.pairs_in > 0, "nothing was read");
    assert_eq!(summary.pairs_in, 135, "upstream's test set is 135 pairs");
}

/// The output names the rest of the pipeline globs for.
#[test]
fn the_outputs_land_where_the_pipeline_expects_them() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (r1, r2) = stage(dir.path());
    let summary = trim_pair(&r1, &r2, &Adapters::None).expect("trim");

    assert!(summary.output_r1.ends_with("test_1_filtered.fastq.gz"));
    assert!(summary.output_r2.ends_with("test_2_filtered.fastq.gz"));
    assert!(summary.output_r1.exists());
    assert!(summary.output_r2.exists());
}
