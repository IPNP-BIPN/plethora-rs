//! Differential test of the `merge_pairs` port against the upstream Perl.
//!
//! `tests/oracle/merge_pairs.pl` is `dpastling/plethora`'s script, vendored
//! unchanged and used only as an oracle. It needs `Math::Random`, which the
//! parity work installs into `.oracle/perl5`:
//!
//! ```text
//! env -u PERL5LIB cpanm --notest -l .oracle/perl5 Math::Random
//! ```
//!
//! Skips when either is missing, so the suite still runs without Perl set up.
//!
//! This is the test that exercises the RANDLIB port end to end: every broken
//! pair in the corpus draws an extension seeded from the MD5 of its own line,
//! and the two implementations have to agree on all of them.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use plethora_core::merge_pairs::{emit, measure};

/// A deterministic generator, so the corpus is identical on every run.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 11
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn oracle_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/oracle")
}

/// The project-local Perl library holding `Math::Random`.
fn perl_lib() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../.oracle/perl5/lib/perl5")
        .canonicalize()
        .unwrap_or_else(|_| PathBuf::from("/nonexistent"))
}

fn oracle_available() -> bool {
    if !oracle_dir().join("merge_pairs.pl").exists() {
        eprintln!("skipping: tests/oracle/merge_pairs.pl is missing");
        return false;
    }
    let ok = Command::new("perl")
        .arg(format!("-I{}", perl_lib().display()))
        .args(["-MMath::Random", "-e", "1"])
        .env_remove("PERL5LIB")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success());
    if !ok {
        eprintln!("skipping: perl with Math::Random is not available");
    }
    ok
}

/// A BEDPE corpus shaped like the one `bamtobed` produces, with the branches
/// `merge_pairs.pl` distinguishes represented in force: proper pairs at a range
/// of inner distances, overlapping reads, pairs beyond the ceiling, mates on
/// different chromosomes, same-strand pairs, and half-mapped pairs.
fn corpus() -> Vec<String> {
    let mut rng = Lcg(0x_9e37_79b9);
    let chroms = ["chr1", "chr2", "chr10"];
    let mut lines = Vec::new();

    for i in 0..3000_u64 {
        let name = format!("frag{i:05}");
        let mapq = rng.below(61);
        let kind = rng.below(12);
        let c1 = chroms[rng.below(3) as usize];
        let start1 = 1000 + rng.below(4_000_000) as i64;
        let len1 = 40 + rng.below(20) as i64;
        let end1 = start1 + len1;

        let (c2, start2, end2, s1, s2) = match kind {
            // Half mapped, one way and the other.
            10 => {
                lines.push(format!(
                    ".\t-1\t-1\t{c1}\t{start1}\t{end1}\t{name}\t0\t.\t+"
                ));
                continue;
            }
            11 => {
                lines.push(format!(
                    "{c1}\t{start1}\t{end1}\t.\t-1\t-1\t{name}\t0\t+\t."
                ));
                continue;
            }
            // Mates on different chromosomes.
            9 => {
                let c2 = chroms[((rng.below(2) + 1) % 3) as usize];
                let s = 1000 + rng.below(4_000_000) as i64;
                (c2, s, s + 50, "+", "-")
            }
            // Same strand, so never a proper pair.
            8 => {
                let gap = rng.below(600) as i64;
                (c1, end1 + gap, end1 + gap + 50, "+", "+")
            }
            // Overlapping reads.
            7 => (c1, start1 + 10, start1 + 60, "+", "-"),
            // Beyond the first-pass ceiling but possibly within the second's.
            6 => {
                let gap = 700 + rng.below(600) as i64;
                (c1, end1 + gap, end1 + gap + 50, "+", "-")
            }
            // The bulk: ordinary proper pairs.
            _ => {
                let gap = rng.below(500) as i64;
                (c1, end1 + gap, end1 + gap + 50, "+", "-")
            }
        };

        lines.push(format!(
            "{c1}\t{start1}\t{end1}\t{c2}\t{start2}\t{end2}\t{name}\t{mapq}\t{s1}\t{s2}"
        ));
    }

    lines
}

#[test]
fn merge_pairs_matches_the_perl() {
    if !oracle_available() {
        return;
    }

    let dir = tempfile::tempdir().expect("tempdir");
    let input = dir.path().join("sample.bed");
    let text = format!("{}\n", corpus().join("\n"));
    std::fs::write(&input, &text).expect("write corpus");

    // The Perl writes sample_edited.bed next to its input.
    let status = Command::new("perl")
        .arg(format!("-I{}", perl_lib().display()))
        .arg(oracle_dir().join("merge_pairs.pl"))
        .arg(&input)
        .env_remove("PERL5LIB")
        .stdout(Stdio::null())
        .stderr(Stdio::inherit())
        .status()
        .expect("run merge_pairs.pl");
    assert!(status.success(), "merge_pairs.pl failed");

    let expected = std::fs::read_to_string(dir.path().join("sample_edited.bed")).expect("read perl output");

    let lines = || text.lines().map(String::from);
    let stats = measure(lines());
    let mut got = Vec::new();
    emit(lines(), &stats, &mut got).expect("emit");
    let got = String::from_utf8(got).expect("utf-8");

    // How much of this actually exercised the random draw, so a corpus that
    // quietly stopped producing broken pairs cannot pass unnoticed.
    let broken = text
        .lines()
        .filter(|l| {
            plethora_core::merge_pairs::BedpeLine::parse(l)
                .is_some_and(|p| !p.is_a_proper_pair(stats.max_inner_distance))
        })
        .count();
    println!(
        "corpus: {} lines, {broken} broken pairs each drawing an MD5-seeded extension; \
         distribution n={} mean={} sd={}",
        text.lines().count(),
        stats.n,
        stats.mean,
        stats.sd
    );
    assert!(broken > 500, "the corpus must exercise the random draw, saw {broken}");

    let e: Vec<&str> = expected.lines().collect();
    let g: Vec<&str> = got.lines().collect();

    assert_eq!(g.len(), e.len(), "line count differs from the Perl");
    let mut differing = 0;
    for (i, (gl, el)) in g.iter().zip(&e).enumerate() {
        if gl != el {
            if differing < 5 {
                eprintln!("line {i}:\n  ours: {gl}\n  perl: {el}");
            }
            differing += 1;
        }
    }
    assert_eq!(differing, 0, "{differing} of {} lines differ from the Perl", e.len());
}

/// The measured distribution must match too, not just the output lines: it is
/// what decides the second pass's ceiling.
#[test]
fn the_fragment_distribution_matches_the_perl() {
    if !oracle_available() {
        return;
    }

    let dir = tempfile::tempdir().expect("tempdir");
    let input = dir.path().join("stats.bed");
    let text = format!("{}\n", corpus().join("\n"));
    std::fs::write(&input, &text).expect("write corpus");

    // Ask the Perl for its own mean and standard deviation by running the same
    // first pass in a one-liner, rather than adding a flag to the vendored
    // script, which must stay byte-identical to upstream.
    let script = r#"
        my $max = 800; my @d;
        open(F, $ARGV[0]) or die;
        while (<F>) { chomp; my @e = split(/\t/);
            next if ($e[0] eq "." || $e[3] eq "." || $e[0] ne $e[3] || $e[8] eq $e[9]);
            next if ($e[4] - $e[2] > $max);
            my $ok = (($e[1] >= $e[4] && $e[1] <= $e[5]) || ($e[2] >= $e[4] && $e[2] <= $e[5]));
            $ok = 1 if (!$ok && !($e[1] > $e[4]));
            push(@d, $e[4] - $e[2]) if $ok;
        }
        close(F);
        my $m = 0; $m += $_ for @d; $m = sprintf("%.0f", $m / scalar(@d));
        my $s = 0; for my $x (@d) { my $c = $x - $m; $s += $c * $c; }
        $s = sprintf("%.0f", sqrt($s / scalar(@d)));
        print scalar(@d), "\t$m\t$s\n";
    "#;
    let out = Command::new("perl")
        .args(["-e", script])
        .arg(&input)
        .output()
        .expect("run the first pass in perl");
    assert!(out.status.success(), "perl first pass failed");

    let text_out = String::from_utf8(out.stdout).expect("utf-8");
    let f: Vec<&str> = text_out.trim().split('\t').collect();
    let (n, mean, sd): (usize, i64, i64) = (
        f[0].parse().expect("n"),
        f[1].parse().expect("mean"),
        f[2].parse().expect("sd"),
    );

    let stats = measure(text.lines().map(String::from));
    assert_eq!(stats.n, n, "sample size");
    assert_eq!(stats.mean, mean, "mean inner distance");
    assert_eq!(stats.sd, sd, "standard deviation");
}

/// A guard on the vendored oracle: if it ever drifts from upstream, the parity
/// claim is about the wrong script.
#[test]
fn the_vendored_oracle_is_upstreams_script() {
    let path = oracle_dir().join("merge_pairs.pl");
    if !path.exists() {
        eprintln!("skipping: oracle script is missing");
        return;
    }
    let text = std::fs::read_to_string(&path).expect("read oracle");
    assert!(text.contains("Math::Random"), "the oracle must use Math::Random");
    assert!(
        text.contains("random_set_seed_from_phrase($eed)"),
        "the oracle must seed from the line's MD5"
    );
    assert!(
        text.contains("my $max_inner_distance = 800;"),
        "the oracle must start from the 800 ceiling"
    );
    let mut hasher = std::process::Command::new("shasum");
    hasher.args(["-a", "256"]).arg(&path);
    if let Ok(out) = hasher.output() {
        eprintln!("oracle sha256: {}", String::from_utf8_lossy(&out.stdout).trim());
    }
}
