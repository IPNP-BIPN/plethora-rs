//! Runs both pipelines on the same input and diffs every intermediate.
//!
//! This is what makes the byte-identical claim checkable instead of
//! declarative. It builds a corpus, hands it to upstream's own `make_bed.sh`
//! and `gc_correction.R` with their own tools underneath, hands the identical
//! corpus to `plethora`, and compares what each left on disk.
//!
//! The corpus is synthetic on purpose. The real reference is a 29 MB BED
//! against an hg38 bowtie2 index nobody has lying around, and a comparison that
//! cannot be run is not a comparison. What is synthesised is only the input:
//! every stage downstream of it is the real one, on both sides.

use std::fmt::Write as _;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use anyhow::{Context, Result, bail};

/// The tools upstream's scripts shell out to. Without these there is nothing
/// to compare against.
const REQUIRED: &[&str] = &["samtools", "bedtools", "perl", "Rscript", "sort", "awk"];

/// Files both pipelines leave behind, in the order they are produced.
const COMPARED: &[&str] = &[
    ".bed",
    "_sorted.bed",
    "_coverage.bed",
    "_read_depth.bed",
    "_gc_correct.txt",
];

pub fn run(root: &Path, upstream: Option<&Path>, keep: bool) -> Result<()> {
    let missing: Vec<&str> = REQUIRED.iter().copied().filter(|t| !have(t)).collect();
    if !missing.is_empty() {
        bail!(
            "cannot compare without upstream's own tools; missing: {}",
            missing.join(", ")
        );
    }

    let work = root.join("target/xtask/compare");
    if work.exists() {
        std::fs::remove_dir_all(&work)?;
    }
    std::fs::create_dir_all(work.join("results"))?;

    let upstream = match upstream {
        Some(dir) => dir.to_path_buf(),
        None => clone_upstream(&work)?,
    };
    if !upstream.join("code/make_bed.sh").exists() {
        bail!(
            "{} does not look like a plethora checkout",
            upstream.display()
        );
    }
    // Work on a copy, so nothing here can modify the user's checkout.
    let staged = stage(root, &upstream, &work)?;

    let corpus = Corpus::build(&work)?;
    println!(
        "corpus: {} domains, {} fragments\n",
        corpus.domains, corpus.fragments
    );

    run_upstream(root, &staged, &work, &corpus)?;
    run_plethora(root, &work, &corpus)?;

    let verdict = diff(&work);
    // Kept on a disagreement, which is when somebody wants to look at it.
    if keep || verdict.is_err() {
        println!("\nworking directory kept at {}", work.display());
    } else {
        let _ = std::fs::remove_dir_all(&work);
    }
    verdict
}

/// What was built, and where.
struct Corpus {
    domains_bed: PathBuf,
    gc_table: PathBuf,
    bam: PathBuf,
    domains: usize,
    fragments: usize,
}

/// A fixed generator, so two runs of `xtask compare` compare the same bytes.
struct Lcg(u64);

impl Lcg {
    fn below(&mut self, n: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) % n
    }
}

impl Corpus {
    /// Conserved domains at a steady depth, domains of interest at three times
    /// that, and GC values spread widely enough that the loess has a curve to
    /// fit rather than a handful of bins.
    fn build(work: &Path) -> Result<Self> {
        const CONSERVED: usize = 120;
        const INTEREST: usize = 30;
        const SPACING: i64 = 2000;
        const LENGTH: i64 = 1000;

        let mut bed = String::new();
        let mut gc = String::new();
        let mut rng = Lcg(0x5eed_0f11);

        for i in 0..CONSERVED + INTEREST {
            let start = 1000 + (i as i64) * SPACING;
            let name = if i < CONSERVED {
                format!("baseline_{i:04}")
            } else {
                format!("NBPF1_CON1_{}", i - CONSERVED)
            };
            writeln!(bed, "chr1\t{start}\t{}\t{name}\t255\t+", start + LENGTH)?;
            // Across the model's window, which the R filters at [0.2, 0.73).
            // Three decimals is the form gc_from_fasta.pl writes.
            let percent = 0.22 + (rng.below(480) as f64) / 1000.0;
            writeln!(gc, "{name}\t{percent:.3}")?;
        }

        let domains_bed = work.join("domains.bed");
        let gc_table = work.join("domains_GC.txt");
        std::fs::write(&domains_bed, bed)?;
        std::fs::write(&gc_table, gc)?;

        // Proper pairs, 50 bp mates 200 bp apart, landing inside a domain.
        let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:400000\n");
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
        let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let mut fragments = 0;
        for i in 0..CONSERVED + INTEREST {
            let domain_start = 1000 + (i as i64) * SPACING;
            let depth = if i < CONSERVED { 30 } else { 90 };
            for _ in 0..depth {
                let start = domain_start + rng.below(700) as i64 + 1;
                let mate = start + 200;
                let name = format!("frag{fragments:06}");
                writeln!(
                    sam,
                    "{name}\t99\tchr1\t{start}\t60\t50M\t=\t{mate}\t250\t{seq}\t{qual}"
                )?;
                writeln!(
                    sam,
                    "{name}\t147\tchr1\t{mate}\t60\t50M\t=\t{start}\t-250\t{seq}\t{qual}"
                )?;
                fragments += 1;
            }
        }

        let bam = work.join("sample.bam");
        let mut child = Command::new("samtools")
            .args(["view", "-b", "-o"])
            .arg(&bam)
            .stdin(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .context("spawning samtools to build the BAM")?;
        child
            .stdin
            .take()
            .expect("stdin was piped")
            .write_all(sam.as_bytes())?;
        let out = child.wait_with_output()?;
        if !out.status.success() {
            bail!(
                "samtools view failed: {}",
                String::from_utf8_lossy(&out.stderr)
            );
        }

        Ok(Self {
            domains_bed,
            gc_table,
            bam,
            domains: CONSERVED + INTEREST,
            fragments,
        })
    }
}

/// `code/make_bed.sh` then `code/gc_correction.R`, exactly as upstream calls
/// them. `make_bed.sh` invokes `code/merge_pairs.pl` by a relative path, so it
/// has to run from the checkout root.
fn run_upstream(root: &Path, upstream: &Path, work: &Path, corpus: &Corpus) -> Result<()> {
    println!("running upstream");
    let prefix = work.join("results/upstream");
    let status = Command::new("./code/make_bed.sh")
        .current_dir(upstream)
        .env("PERL5LIB", perl5lib(root))
        .arg("-r")
        .arg(&corpus.domains_bed)
        .args(["-p", "paired"])
        .arg("-b")
        .arg(&corpus.bam)
        .arg("-o")
        .arg(&prefix)
        .status()
        .context("running upstream's make_bed.sh")?;
    if !status.success() {
        bail!("upstream make_bed.sh exited with {status}");
    }

    // gc_correction.R calls as.tbl(), which dplyr made defunct in 1.2.0. The
    // call is a no-op on a data.frame that dplyr verbs accept either way, so
    // dropping it changes nothing but lets the script run on a current dplyr.
    // Anything else in the file is left alone.
    let script = std::fs::read_to_string(upstream.join("code/gc_correction.R"))?;
    let patched = script.replace("X <- as.tbl(X)", "X <- as.data.frame(X)");
    if patched == script {
        println!("  note: as.tbl() not found, running gc_correction.R unmodified");
    }
    let script_path = work.join("gc_correction.R");
    std::fs::write(&script_path, patched)?;

    let read_depth = work.join("results/upstream_read_depth.bed");
    let status = Command::new("Rscript")
        .arg(&script_path)
        .arg(&read_depth)
        .arg(&corpus.gc_table)
        .status()
        .context("running upstream's gc_correction.R")?;
    if !status.success() {
        bail!("upstream gc_correction.R exited with {status}");
    }
    Ok(())
}

fn run_plethora(root: &Path, work: &Path, corpus: &Corpus) -> Result<()> {
    println!("\nrunning plethora");
    let prefix = work.join("results/plethora");

    let cargo = std::env::var("CARGO").unwrap_or_else(|_| "cargo".to_string());
    let mut coverage = Command::new(&cargo);
    coverage
        .current_dir(root)
        .args(["run", "--quiet", "-p", "plethora", "--", "coverage"])
        .arg("-r")
        .arg(&corpus.domains_bed)
        .args(["-p", "paired"])
        .arg("-b")
        .arg(&corpus.bam)
        .arg("-o")
        .arg(&prefix)
        .arg("--temp")
        .arg(work);
    let status = coverage.status().context("running plethora coverage")?;
    if !status.success() {
        bail!("plethora coverage exited with {status}");
    }

    let status = Command::new(&cargo)
        .current_dir(root)
        .args(["run", "--quiet", "-p", "plethora", "--", "gc-correct"])
        .arg(work.join("results/plethora_read_depth.bed"))
        .arg(&corpus.gc_table)
        .status()
        .context("running plethora gc-correct")?;
    if !status.success() {
        bail!("plethora gc-correct exited with {status}");
    }
    Ok(())
}

/// How far two numbers may differ and still count as agreement.
///
/// Not a fudge factor: `DIVERGENCES.md` measures the reference BLAS residual
/// inside `loess` at a worst relative 2.0e-14 over 623,699 real domains, and R
/// compiled on two machines disagrees with itself by more than that. This sits
/// two orders above the observed worst and ten orders below anything that could
/// move a copy-number call.
const TOLERANCE: f64 = 1e-12;

/// The verdict for one file.
enum Verdict {
    Identical(usize),
    /// Same values, last-digit apart, with the worst relative difference and
    /// where it was.
    WithinTolerance {
        rows: usize,
        worst: f64,
        at: String,
    },
    Differs(String),
}

/// Compares byte for byte first, and only when that fails asks whether the two
/// files are the same numbers written differently.
fn diff(work: &Path) -> Result<()> {
    println!("\n{:<20} {:>10}  result", "file", "bytes");
    let mut failures = 0;
    let mut tolerated = 0;

    for suffix in COMPARED {
        let theirs = work.join(format!("results/upstream{suffix}"));
        let ours = work.join(format!("results/plethora{suffix}"));

        let (a, b) = match (std::fs::read(&theirs), std::fs::read(&ours)) {
            (Ok(a), Ok(b)) => (a, b),
            (Err(_), _) => {
                println!("{suffix:<20} {:>10}  upstream did not write it", "-");
                failures += 1;
                continue;
            }
            (_, Err(_)) => {
                println!("{suffix:<20} {:>10}  plethora did not write it", "-");
                failures += 1;
                continue;
            }
        };

        match compare_bytes(&a, &b) {
            Verdict::Identical(n) => println!("{suffix:<20} {n:>10}  identical"),
            Verdict::WithinTolerance { rows, worst, at } => {
                tolerated += 1;
                println!(
                    "{suffix:<20} {:>10}  equal to {worst:.1e} relative ({rows} rows differ in \
                     their last digits, worst at {at})",
                    a.len()
                );
            }
            Verdict::Differs(why) => {
                failures += 1;
                println!("{suffix:<20} {:>10}  DIFFERS", a.len());
                for line in why.lines() {
                    println!("    {line}");
                }
            }
        }
    }

    // _edited.bed and _temp.bed are removed by both, so neither can be read
    // back. _sorted.bed is a total order over _edited.bed's lines, so an
    // identical _sorted.bed says the two agreed on the same multiset of
    // intervals, which is the part that carries into the result.
    println!("\nnot compared: _edited.bed and _temp.bed, which both pipelines delete");

    if failures > 0 {
        bail!("{failures} of {} files disagree", COMPARED.len());
    }
    let identical = COMPARED.len() - tolerated;
    println!(
        "\n{identical} of {} files are byte-identical",
        COMPARED.len()
    );
    if tolerated > 0 {
        println!(
            "{tolerated} agree numerically to better than {TOLERANCE:.0e} relative, which is the\n\
             BLAS residual inside loess that DIVERGENCES.md documents."
        );
    }
    Ok(())
}

fn compare_bytes(a: &[u8], b: &[u8]) -> Verdict {
    if a == b {
        return Verdict::Identical(a.len());
    }
    let (Ok(theirs), Ok(ours)) = (std::str::from_utf8(a), std::str::from_utf8(b)) else {
        return Verdict::Differs("not both UTF-8".to_string());
    };
    compare_tables(theirs, ours)
}

/// Field by field: the row keys must match exactly, numeric fields may differ
/// within tolerance, and anything else is a real disagreement.
fn compare_tables(theirs: &str, ours: &str) -> Verdict {
    let a: Vec<&str> = theirs.lines().collect();
    let b: Vec<&str> = ours.lines().collect();
    if a.len() != b.len() {
        return Verdict::Differs(format!("{} rows against {}", a.len(), b.len()));
    }

    let mut rows = 0;
    let mut worst = 0.0f64;
    let mut at = String::new();

    for (i, (x, y)) in a.iter().zip(&b).enumerate() {
        if x == y {
            continue;
        }
        let fx: Vec<&str> = x.split('\t').collect();
        let fy: Vec<&str> = y.split('\t').collect();
        if fx.len() != fy.len() {
            return Verdict::Differs(format!(
                "line {}: {} fields against {}\n  upstream: {x}\n  plethora: {y}",
                i + 1,
                fx.len(),
                fy.len()
            ));
        }
        rows += 1;

        for (j, (p, q)) in fx.iter().zip(&fy).enumerate() {
            if p == q {
                continue;
            }
            let (Ok(p), Ok(q)) = (p.parse::<f64>(), q.parse::<f64>()) else {
                return Verdict::Differs(format!(
                    "line {}, field {}: {p} against {q}\n  upstream: {x}\n  plethora: {y}",
                    i + 1,
                    j + 1
                ));
            };
            let relative = (p - q).abs() / p.abs().max(f64::MIN_POSITIVE);
            if relative > worst {
                worst = relative;
                // The first field is the domain name on every table here.
                at = format!("{} field {}", fx.first().unwrap_or(&"?"), j + 1);
            }
            if relative > TOLERANCE {
                return Verdict::Differs(format!(
                    "line {}, field {}: relative {relative:.2e} exceeds {TOLERANCE:.0e}\n  \
                     upstream: {x}\n  plethora: {y}",
                    i + 1,
                    j + 1
                ));
            }
        }
    }

    Verdict::WithinTolerance { rows, worst, at }
}

fn clone_upstream(work: &Path) -> Result<PathBuf> {
    let dir = work.join("upstream");
    println!("cloning dpastling/plethora");
    let status = Command::new("git")
        .args(["clone", "--quiet", "--depth", "1"])
        .arg("https://github.com/dpastling/plethora")
        .arg(&dir)
        .status()
        .context("cloning upstream")?;
    if !status.success() {
        bail!("git clone exited with {status}");
    }
    Ok(dir)
}

/// Copies upstream's `code/` into the working directory, and makes sure the
/// Perl in it can actually load `Math::Random`.
///
/// `merge_pairs.pl` carries `#!/usr/bin/perl`, which on macOS is the system
/// Perl, while a `cpanm -l` build of `Math::Random` is compiled against
/// whichever Perl ran cpanm. Loading one into the other fails at the XS
/// handshake. Where that happens the shebang is rewritten to `env perl`, which
/// changes which interpreter runs and nothing about the program; it is reported
/// so a comparison never quietly runs something other than what is on disk.
fn stage(root: &Path, upstream: &Path, work: &Path) -> Result<PathBuf> {
    let staged = work.join("upstream-code");
    std::fs::create_dir_all(staged.join("code"))?;
    for entry in std::fs::read_dir(upstream.join("code"))? {
        let entry = entry?;
        if entry.file_type()?.is_file() {
            let to = staged.join("code").join(entry.file_name());
            std::fs::copy(entry.path(), &to)?;
            // getopts and the pipeline both need these executable.
            #[cfg(unix)]
            {
                use std::os::unix::fs::PermissionsExt as _;
                std::fs::set_permissions(&to, std::fs::Permissions::from_mode(0o755))?;
            }
        }
    }

    let script = staged.join("code/merge_pairs.pl");
    let text = std::fs::read_to_string(&script)?;
    let shebang = text.lines().next().unwrap_or_default().to_string();
    let interpreter = shebang.strip_prefix("#!").unwrap_or("perl").trim();
    let interpreter = interpreter
        .split_whitespace()
        .next_back()
        .unwrap_or("perl")
        .to_string();

    if loads_math_random(root, &interpreter) {
        println!("upstream's {interpreter} has Math::Random");
        return Ok(staged);
    }
    if loads_math_random(root, "perl") {
        let rewritten = text.replacen(&shebang, "#!/usr/bin/env perl", 1);
        std::fs::write(&script, rewritten)?;
        println!(
            "note: {interpreter} cannot load Math::Random, so merge_pairs.pl runs\n\
             \x20     under `env perl` instead. Only the interpreter changes."
        );
        return Ok(staged);
    }
    bail!(
        "merge_pairs.pl needs Math::Random and no Perl here can load it. The\n\
         parity tests install it into .oracle/perl5:\n\n    \
         env -u PERL5LIB cpanm --notest -l .oracle/perl5 Math::Random"
    );
}

/// The project-local Perl library the parity tests install `Math::Random`
/// into, so `xtask compare` and `cargo test` share one oracle.
fn perl5lib(root: &Path) -> PathBuf {
    root.join(".oracle/perl5/lib/perl5")
}

fn loads_math_random(root: &Path, interpreter: &str) -> bool {
    Command::new(interpreter)
        .args(["-MMath::Random", "-e", "1"])
        .env("PERL5LIB", perl5lib(root))
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

fn have(tool: &str) -> bool {
    Command::new(tool)
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}
