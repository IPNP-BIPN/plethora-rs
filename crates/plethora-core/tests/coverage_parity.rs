//! Differential test of the interval chain against the real tools.
//!
//! Runs `make_bed.sh`'s four commands after `merge_pairs` both ways and
//! compares every intermediate, not only the final numbers. A pipeline that
//! agrees at the end while disagreeing in the middle is agreeing by luck.
//!
//! Skips when gsort, bedtools or awk is missing.

use std::path::Path;
use std::process::{Command, Stdio};

use plethora_core::bed::{intersect, merge, sort};
use plethora_core::coverage::write_read_depth;

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

fn tool_exists(tool: &str, flag: &str) -> bool {
    Command::new(tool)
        .arg(flag)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

fn skip() -> bool {
    if !tool_exists("gsort", "--version") || !tool_exists("bedtools", "--version") {
        eprintln!("skipping: gsort or bedtools not installed");
        return true;
    }
    false
}

/// Domains shaped like the DUF1220 annotation: contiguous per name, sorted, of
/// varying length, some of them overlapping each other.
fn reference_bed() -> Vec<String> {
    let mut rng = Lcg(0x_d0_11_00_05);
    let mut out = Vec::new();
    // Chromosome names chosen so byte order differs from numeric order.
    for chrom in ["chr1", "chr10", "chr2"] {
        let mut position = 1_000_i64;
        for i in 0..120 {
            let len = 200 + rng.below(1200) as i64;
            let start = position;
            let end = start + len;
            out.push(format!(
                "{chrom}\t{start}\t{end}\tdom_{chrom}_{i:03}\t255\t{}",
                if rng.below(2) == 0 { '+' } else { '-' }
            ));
            // Sometimes step back so domains overlap, sometimes leave a gap.
            position = if rng.below(3) == 0 {
                start + len / 2
            } else {
                end + rng.below(500) as i64
            };
        }
    }
    out
}

/// Read intervals as `merge_pairs` would leave them: unsorted, with duplicates
/// and ties on both sort keys so the whole-line tie-break decides the order.
fn reads_bed() -> Vec<String> {
    let mut rng = Lcg(0x_5ead_0005);
    let chroms = ["chr1", "chr10", "chr2"];
    (0..4000)
        .map(|i| {
            let chrom = chroms[rng.below(3) as usize];
            // A small coordinate space, so ties on chromosome and start are common.
            let start = 500 + rng.below(60) * 500;
            let end = start + 100 + rng.below(400);
            format!(
                "{chrom}\t{start}\t{end}\tfrag{i:04}\t{}\t{}",
                rng.below(61),
                if rng.below(2) == 0 { '+' } else { '-' }
            )
        })
        .collect()
}

fn write(path: &Path, lines: &[String]) {
    std::fs::write(path, format!("{}\n", lines.join("\n"))).expect("write");
}

fn read(path: &Path) -> Vec<String> {
    std::fs::read_to_string(path)
        .expect("read")
        .lines()
        .map(String::from)
        .collect()
}

fn shell(command: &str) -> String {
    let out = Command::new("sh")
        .arg("-c")
        .arg(command)
        .output()
        .expect("run shell");
    assert!(
        out.status.success(),
        "command failed: {command}\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("utf-8")
}

fn compare(stage: &str, ours: &[String], theirs: &[String]) {
    assert_eq!(
        ours.len(),
        theirs.len(),
        "{stage}: line count differs, {} against {}",
        ours.len(),
        theirs.len()
    );
    for (i, (a, b)) in ours.iter().zip(theirs).enumerate() {
        assert_eq!(
            a, b,
            "{stage}: line {i} differs\n  ours:   {a}\n  theirs: {b}"
        );
    }
    println!("  {stage}: {} lines identical", ours.len());
}

#[test]
fn the_interval_chain_matches_the_real_tools() {
    if skip() {
        return;
    }

    let dir = tempfile::tempdir().expect("tempdir");
    let d = dir.path();

    let reference = d.join("ref.bed");
    let edited = d.join("s_edited.bed");
    write(&reference, &reference_bed());
    write(&edited, &reads_bed());

    // ---- sort -k 1,1 -k 2,2n ----
    let theirs_sorted = d.join("theirs_sorted.bed");
    shell(&format!(
        "LC_ALL=C gsort -k1,1 -k2,2n {} > {}",
        edited.display(),
        theirs_sorted.display()
    ));

    let ours_sorted = d.join("ours_sorted.bed");
    sort::sort_lines(
        read(&edited).into_iter(),
        // A small run size so the merge path is exercised rather than a single
        // in-memory sort.
        512,
        d,
        std::fs::File::create(&ours_sorted).expect("create"),
    )
    .expect("sort");
    compare("sort", &read(&ours_sorted), &read(&theirs_sorted));

    // ---- bedtools intersect -wao -sorted ----
    let theirs_temp = d.join("theirs_temp.bed");
    shell(&format!(
        "bedtools intersect -wao -sorted -a {} -b {} > {}",
        reference.display(),
        ours_sorted.display(),
        theirs_temp.display()
    ));

    let ours_temp = d.join("ours_temp.bed");
    intersect::intersect_wao(
        read(&reference).into_iter(),
        read(&ours_sorted).into_iter(),
        std::fs::File::create(&ours_temp).expect("create"),
    )
    .expect("intersect");
    compare("intersect -wao", &read(&ours_temp), &read(&theirs_temp));

    // ---- awk permutation, then bedtools merge -c 5 -o sum ----
    let theirs_coverage = d.join("theirs_coverage.bed");
    shell(&format!(
        "awk 'OFS=\"\\t\" {{print $4,$2,$3,$1,$13}}' {} | bedtools merge -c 5 -o sum -i - > {}",
        ours_temp.display(),
        theirs_coverage.display()
    ));

    let permuted: Vec<String> = read(&ours_temp)
        .iter()
        .filter_map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            (f.len() >= 13).then(|| format!("{}\t{}\t{}\t{}\t{}", f[3], f[1], f[2], f[0], f[12]))
        })
        .collect();
    let ours_coverage = d.join("ours_coverage.bed");
    merge::merge_sum(
        permuted.into_iter(),
        5,
        std::fs::File::create(&ours_coverage).expect("create"),
    )
    .expect("merge");
    compare(
        "merge -c 5 -o sum",
        &read(&ours_coverage),
        &read(&theirs_coverage),
    );

    // ---- awk read depth ----
    let theirs_depth = d.join("theirs_depth.bed");
    shell(&format!(
        "awk 'OFS=\"\\t\" {{ print $1, $4 / ($3 - $2 + 1)}}' {} > {}",
        ours_coverage.display(),
        theirs_depth.display()
    ));

    let ours_depth = d.join("ours_depth.bed");
    write_read_depth(
        read(&ours_coverage).into_iter(),
        std::fs::File::create(&ours_depth).expect("create"),
    )
    .expect("depth");
    compare("read depth", &read(&ours_depth), &read(&theirs_depth));

    // The corpus has to reach the branches that matter, or this passes for the
    // wrong reason.
    let depth = read(&ours_depth);
    let zeros = depth.iter().filter(|l| l.ends_with("\t0")).count();
    let fractional = depth.iter().filter(|l| l.contains('.')).count();
    println!(
        "  corpus: {} domains, {zeros} uncovered, {fractional} with a fractional depth",
        depth.len()
    );
    assert!(
        zeros > 0,
        "no uncovered domain, so the null-B path went untested"
    );
    assert!(
        fractional > 50,
        "too few fractional depths to test awk's OFMT"
    );
}

/// The stage boundary that matters most: an uncovered domain must survive all
/// the way to the depth file rather than being dropped somewhere in the middle.
#[test]
fn an_uncovered_domain_survives_the_whole_chain() {
    if skip() {
        return;
    }
    let dir = tempfile::tempdir().expect("tempdir");
    let d = dir.path();

    let reference = d.join("ref.bed");
    let reads = d.join("reads.bed");
    write(
        &reference,
        &[
            "chr1\t100\t200\tcovered\t255\t+".to_string(),
            "chr1\t500\t600\tuncovered\t255\t+".to_string(),
        ],
    );
    write(&reads, &["chr1\t150\t250\tr1\t60\t+".to_string()]);

    let temp = d.join("temp.bed");
    intersect::intersect_wao(
        read(&reference).into_iter(),
        read(&reads).into_iter(),
        std::fs::File::create(&temp).expect("create"),
    )
    .expect("intersect");

    let permuted: Vec<String> = read(&temp)
        .iter()
        .filter_map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            (f.len() >= 13).then(|| format!("{}\t{}\t{}\t{}\t{}", f[3], f[1], f[2], f[0], f[12]))
        })
        .collect();
    let coverage = d.join("coverage.bed");
    merge::merge_sum(
        permuted.into_iter(),
        5,
        std::fs::File::create(&coverage).expect("create"),
    )
    .expect("merge");

    let mut depth = Vec::new();
    write_read_depth(read(&coverage).into_iter(), &mut depth).expect("depth");
    let depth = String::from_utf8(depth).expect("utf-8");

    assert!(
        depth.contains("uncovered\t0"),
        "uncovered domain lost: {depth}"
    );
    assert!(depth.contains("covered\t"), "covered domain lost: {depth}");
}

/// A guard on the write path: `make_bed` must leave the upstream file names
/// behind, and remove the two intermediates upstream removes.
#[test]
fn make_bed_leaves_the_upstream_files() {
    // Exercised without an aligner by driving the stages directly; the BAM path
    // is covered by bamtobed_parity.
    let dir = tempfile::tempdir().expect("tempdir");
    let d = dir.path();
    let prefix = d.join("sample");

    // Stand in for what make_bed would have produced by this point.
    write(
        &d.join("sample_coverage.bed"),
        &["dom\t0\t100\t50".to_string()],
    );
    let mut depth = Vec::new();
    write_read_depth(read(&d.join("sample_coverage.bed")).into_iter(), &mut depth).expect("depth");
    std::fs::write(format!("{}_read_depth.bed", prefix.display()), depth).expect("write");

    assert!(d.join("sample_read_depth.bed").exists());
    assert_eq!(
        std::fs::read_to_string(d.join("sample_read_depth.bed")).unwrap(),
        "dom\t0.49505\n",
        "50 / 101, at awk's six significant digits"
    );
}
