//! Differential test of the GC chain against the upstream Perl and R.
//!
//! `tests/oracle/gc_from_fasta.pl` and `tests/oracle/gc_correction.R` are
//! `dpastling/plethora`'s own scripts, vendored unchanged and used only as
//! oracles. The R one needs `dplyr`.
//!
//! This closes the loop on the loess port: `gc_correction.R` fits one, and its
//! fitted values become the correction factor for every domain, so agreeing
//! here means the whole `plethora-compat::loess` chain agrees in the place it
//! is actually used.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use plethora_core::gc::correction::{Row, correct};
use plethora_core::gc::from_fasta::gc_from_fasta;

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

fn have(program: &str, args: &[&str]) -> bool {
    Command::new(program)
        .args(args)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

fn have_r_with_dplyr() -> bool {
    let ok = have("Rscript", &["-e", "library(dplyr)"]);
    if !ok {
        eprintln!("skipping: Rscript with dplyr is not available");
    }
    ok
}

/// A FASTA covering the readings that decide the GC figure: plain uppercase,
/// soft-masked lowercase, uppercase N, lowercase n, and ambiguity codes.
fn corpus_fasta() -> String {
    let mut rng = Lcg(0x_9c_c0_11_00);
    let alphabet = b"ACGTacgtNnRYKMSW";
    let mut out = String::new();
    for i in 0..200 {
        out.push_str(&format!(">dom{i:04}\n"));
        let len = 60 + rng.below(200) as usize;
        let mut seq = String::new();
        for j in 0..len {
            seq.push(char::from(alphabet[rng.below(alphabet.len() as u64) as usize]));
            // Wrap, so the per-line accumulation is exercised.
            if j % 60 == 59 {
                seq.push('\n');
            }
        }
        out.push_str(&seq);
        out.push('\n');
    }
    out
}

#[test]
fn gc_from_fasta_matches_the_perl() {
    let script = oracle_dir().join("gc_from_fasta.pl");
    if !script.exists() || !have("perl", &["-e", "1"]) {
        eprintln!("skipping: perl or the oracle script is missing");
        return;
    }

    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = dir.path().join("domains.fa");
    let text = corpus_fasta();
    std::fs::write(&fasta, &text).expect("write fasta");

    let out = Command::new("perl")
        .arg(&script)
        .arg(&fasta)
        .output()
        .expect("run gc_from_fasta.pl");
    assert!(out.status.success(), "gc_from_fasta.pl failed");

    // The Perl iterates a hash, so its output order is whatever that run's hash
    // randomisation gives. Compare as a mapping, and note the ordering
    // divergence rather than pretending the orders agree.
    let expected: HashMap<String, String> = String::from_utf8_lossy(&out.stdout)
        .lines()
        .filter_map(|l| l.split_once('\t'))
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect();

    let rows = gc_from_fasta(text.as_bytes()).expect("gc from fasta");
    assert_eq!(rows.len(), expected.len(), "sequence count differs");

    let mut soft_masked = 0;
    let mut ambiguous = 0;
    for row in &rows {
        let (_, value) = row.to_line().split_once('\t').map(|(a, b)| (a.to_string(), b.to_string())).unwrap();
        let want = expected.get(&row.name).unwrap_or_else(|| panic!("{} missing from the Perl output", row.name));
        assert_eq!(&value, want, "{}: GC fraction differs", row.name);
        soft_masked += row.counts.soft_masked;
        ambiguous += row.counts.ambiguous;
    }

    println!(
        "corpus: {} sequences, {soft_masked} soft-masked and {ambiguous} ambiguous bases counted as GC",
        rows.len()
    );
    assert!(soft_masked > 100, "the corpus must exercise the lowercase reading");
    assert!(ambiguous > 100, "the corpus must exercise the ambiguity-code reading");
}

/// A read-depth table and a GC table shaped like the real ones: mostly
/// conserved `baseline_*` domains spread across the GC window, plus DUF1220
/// domains at a range of copy numbers.
fn corpus_tables() -> (Vec<(String, f64)>, HashMap<String, f64>) {
    let mut rng = Lcg(0x_9c_c0_22_00);
    let mut depth = Vec::new();
    let mut gc = HashMap::new();

    // Conserved regions: the calibration rests on these.
    for i in 0..600 {
        let name = format!("baseline_{i:04}");
        // GC across the fitted window and a little outside it on both sides.
        let percent = 0.15 + rng.below(650) as f64 / 1000.0;
        // Coverage with a GC-dependent bias, which is what the model corrects.
        let bias = 1.0 - 3.0 * (percent - 0.42) * (percent - 0.42);
        let coverage = 30.0 * bias + rng.below(200) as f64 / 100.0;
        depth.push((name.clone(), coverage));
        gc.insert(name, percent);
    }

    // A few domains whose names contain "uc" without starting with it: they
    // set the haploid unit but not the curve.
    for i in 0..20 {
        let name = format!("NBPF1_uc_{i:02}");
        let percent = 0.30 + rng.below(300) as f64 / 1000.0;
        depth.push((name.clone(), 28.0 + rng.below(400) as f64 / 100.0));
        gc.insert(name, percent);
    }

    // The domains of interest, at varying copy number.
    for i in 0..300 {
        let name = format!("NBPF{}_CON1_{i}", 1 + i % 20);
        let percent = 0.25 + rng.below(400) as f64 / 1000.0;
        let copies = 1 + rng.below(12);
        depth.push((name.clone(), 15.0 * copies as f64 + rng.below(300) as f64 / 100.0));
        gc.insert(name, percent);
    }

    // Some with no coverage at all, which must survive to the output.
    for i in 0..15 {
        let name = format!("NBPF_empty_{i:02}");
        depth.push((name.clone(), 0.0));
        gc.insert(name, 0.40);
    }

    (depth, gc)
}

#[test]
fn gc_correction_matches_the_r() {
    let script = oracle_dir().join("gc_correction.R");
    if !script.exists() || !have_r_with_dplyr() {
        return;
    }

    let dir = tempfile::tempdir().expect("tempdir");
    let (depth, gc) = corpus_tables();

    // The R script derives its output name from the input, so the input has to
    // be spelled the way make_bed.sh spells it.
    let depth_path = dir.path().join("sample_read_depth.bed");
    let gc_path = dir.path().join("domains_GC.txt");
    let mut depth_text = String::new();
    for (d, c) in &depth {
        use std::fmt::Write as _;
        writeln!(depth_text, "{d}\t{}", plethora_compat::awk::print_number(*c)).expect("format");
    }
    std::fs::write(&depth_path, &depth_text).expect("write depth");
    // Written the way gc_from_fasta.pl writes it, at %.15g. Not cosmetic:
    // R's numeric parser is not correctly rounded, and reads a
    // seventeen-digit "0.42500000000000004" back as the neighbouring double
    // 0.425, which then rounds to 0.42 where the exact value rounds to 0.43.
    // Real GC files never carry more than fifteen digits, so writing them any
    // other way would manufacture a disagreement the pipeline cannot have.
    let mut gc_text = String::new();
    for (d, p) in &gc {
        use std::fmt::Write as _;
        writeln!(gc_text, "{d}\t{}", plethora_compat::awk::format_g(*p, 15)).expect("format");
    }
    std::fs::write(&gc_path, &gc_text).expect("write gc");

    // Both sides must start from the same doubles, so read them back from the
    // files rather than reusing the values that generated them. This is not
    // pedantry: 0.15 + 275.0 / 1000.0 and the parse of "0.425" are different
    // doubles, and R's round sends them to 0.43 and 0.42 respectively. Feeding
    // one side the computed value and the other the written text would report a
    // disagreement that is entirely the test's own doing.
    let read_pairs = |path: &Path| -> Vec<(String, f64)> {
        std::fs::read_to_string(path)
            .expect("read table")
            .lines()
            .filter_map(|l| l.split_once('\t'))
            .map(|(d, v)| (d.to_string(), v.parse().expect("a numeric column")))
            .collect()
    };
    let depth: Vec<(String, f64)> = read_pairs(&depth_path);
    let gc: HashMap<String, f64> = read_pairs(&gc_path).into_iter().collect();

    // The vendored script cannot run as it stands: `as.tbl()` was deprecated in
    // dplyr 1.0.0 and is defunct in 1.2.1, so gc_correction.R dies before
    // reading a single row on any dplyr from the last five years. The vendored
    // copy stays byte-identical to upstream and the single substitution is
    // applied here, in the open, so what is being compared against is exactly
    // upstream plus that one line.
    let pristine = std::fs::read_to_string(&script).expect("read the oracle");
    let patched = pristine.replace("X <- as.tbl(X)", "X <- X");
    assert_eq!(
        pristine.matches("as.tbl(X)").count(),
        1,
        "expected exactly one as.tbl call to patch"
    );
    let runnable = dir.path().join("gc_correction_runnable.R");
    std::fs::write(&runnable, &patched).expect("write the patched oracle");

    let out = Command::new("Rscript")
        .arg(&runnable)
        .arg(&depth_path)
        .arg(&gc_path)
        .stderr(Stdio::inherit())
        .output()
        .expect("run gc_correction.R");
    assert!(out.status.success(), "gc_correction.R failed");

    let produced = dir.path().join("sample_gc_correct.txt");
    let expected = std::fs::read_to_string(&produced).expect("read the R output");
    let expected: HashMap<String, Vec<String>> = expected
        .lines()
        .skip(1)
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].to_string(), f[1..].iter().map(|s| (*s).to_string()).collect())
        })
        .collect();

    let rows: Vec<Row> = correct(&depth, &gc).expect("correct");
    assert_eq!(rows.len(), expected.len(), "row count differs from R");

    let mut worst_rel = 0.0_f64;
    let mut worst_where = String::new();
    let mut exact = 0_usize;

    for row in &rows {
        let want = expected
            .get(&row.domain)
            .unwrap_or_else(|| panic!("{} missing from the R output", row.domain));

        // percent.gc and k.gc are compared as text: they come from R's own
        // rounding and formatting, so any difference there is a real one.
        assert_eq!(
            plethora_compat::rmath::format_as_r(row.percent_gc),
            want[1],
            "{}: percent.gc",
            row.domain
        );

        for (label, got, text) in [
            ("k.gc", row.k_gc, &want[2]),
            ("corrected.coverage", row.corrected_coverage, &want[3]),
        ] {
            let want_value: f64 = text.parse().expect("a numeric column");
            if got == want_value {
                exact += 1;
                continue;
            }
            let rel = (got - want_value).abs() / want_value.abs().max(1e-300);
            if rel > worst_rel {
                worst_rel = rel;
                worst_where = format!("{} {label}: got {got:.17e}, want {want_value:.17e}", row.domain);
            }
        }
    }

    println!(
        "gc correction: {} domains, {exact} values bit-exact, worst relative {worst_rel:.3e} ({worst_where})",
        rows.len()
    );
    // The residual is the BLAS difference inside the loess, carried through
    // k.gc into the corrected coverage. See plethora_compat::loess::blas.
    assert!(
        worst_rel < 1e-9,
        "drifted {worst_rel:.3e} from R, far beyond the loess BLAS gap: {worst_where}"
    );
}

/// A guard on the vendored oracles.
#[test]
fn the_vendored_oracles_are_upstreams_scripts() {
    let r = oracle_dir().join("gc_correction.R");
    if r.exists() {
        let text = std::fs::read_to_string(&r).expect("read");
        assert!(text.contains("loess(y ~ x)"), "the oracle must fit a loess");
        assert!(text.contains("min.gc <- 0.2"), "the oracle must use the 0.2 floor");
        assert!(
            text.contains(r#"grepl("^((baseline)|(uc))", domain)"#),
            "the oracle must use the anchored pattern for the model"
        );
        assert!(
            text.contains(r#"grepl("((baseline)|(uc))", domain)"#),
            "the oracle must use the unanchored pattern for the haploid unit"
        );
    }
}
