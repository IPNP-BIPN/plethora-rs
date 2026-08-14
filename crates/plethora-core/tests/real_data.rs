//! Validation against a real upstream run, when one is available.
//!
//! The synthetic corpora elsewhere in this suite are built to reach particular
//! branches. This one does the opposite: it takes files a real cohort produced
//! by running the upstream scripts and asks whether the port reproduces them.
//!
//! The inputs are too large to vendor, so the test reads their paths from the
//! environment and skips without them:
//!
//! ```text
//! PLETHORA_REAL_READ_DEPTH=/path/to/S36742_read_depth.bed \
//! PLETHORA_REAL_GC=/path/to/hg38_duf_full_domains_v2.3_GC.txt \
//! PLETHORA_REAL_EXPECTED=/path/to/S36742_gc_correct.txt \
//!     cargo test -p plethora-core --test real_data -- --nocapture
//! ```
//!
//! Gzipped inputs are accepted, since that is how such files are usually kept.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};
use std::path::PathBuf;

use plethora_core::gc::correction::{Row, correct};

/// Reads a file, transparently decompressing a gzip member.
fn read_lines(path: &PathBuf) -> Vec<String> {
    let mut file = std::fs::File::open(path).expect("open the input");
    let mut magic = [0_u8; 2];
    let gzipped = file.read_exact(&mut magic).is_ok() && magic == [0x1f, 0x8b];

    let file = std::fs::File::open(path).expect("reopen the input");
    let reader: Box<dyn Read> = if gzipped {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    BufReader::new(reader).lines().map_while(Result::ok).collect()
}

fn env_path(key: &str) -> Option<PathBuf> {
    std::env::var_os(key).map(PathBuf::from)
}

#[test]
fn reproduces_a_real_gc_correction() {
    let (Some(depth_path), Some(gc_path), Some(expected_path)) = (
        env_path("PLETHORA_REAL_READ_DEPTH"),
        env_path("PLETHORA_REAL_GC"),
        env_path("PLETHORA_REAL_EXPECTED"),
    ) else {
        eprintln!(
            "skipping: set PLETHORA_REAL_READ_DEPTH, PLETHORA_REAL_GC and \
             PLETHORA_REAL_EXPECTED to a real upstream run"
        );
        return;
    };

    let depth: Vec<(String, f64)> = read_lines(&depth_path)
        .iter()
        .filter_map(|l| l.split_once('\t'))
        .map(|(d, c)| (d.to_string(), c.parse().expect("a coverage figure")))
        .collect();

    let gc: HashMap<String, f64> = read_lines(&gc_path)
        .iter()
        .filter_map(|l| l.split_once('\t'))
        .map(|(d, p)| (d.to_string(), p.parse().expect("a GC fraction")))
        .collect();

    // The upstream output, skipping its header.
    let expected: HashMap<String, Vec<String>> = read_lines(&expected_path)
        .iter()
        .skip(1)
        .filter_map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f.len() >= 5).then(|| {
                (
                    f[0].to_string(),
                    f[1..].iter().map(|s| (*s).to_string()).collect(),
                )
            })
        })
        .collect();

    assert!(!depth.is_empty(), "the read-depth file is empty");
    assert!(!expected.is_empty(), "the expected file is empty");

    let rows: Vec<Row> = correct(&depth, &gc).expect("correct");

    println!(
        "read depth: {} domains, GC table: {} domains, upstream output: {} rows, ours: {} rows",
        depth.len(),
        gc.len(),
        expected.len(),
        rows.len()
    );
    assert_eq!(rows.len(), expected.len(), "row count differs from upstream");

    let mut gc_text_mismatches = 0;
    let mut worst_rel = 0.0_f64;
    let mut worst_where = String::new();
    let mut exact = 0_usize;
    let mut compared = 0_usize;

    for row in &rows {
        let Some(want) = expected.get(&row.domain) else {
            panic!("{} is missing from the upstream output", row.domain);
        };

        // percent.gc is R's own rounding, printed by R. Compared as text.
        if plethora_compat::rmath::format_as_r(row.percent_gc) != want[1] {
            gc_text_mismatches += 1;
        }

        for (label, got, text) in [
            ("k.gc", row.k_gc, &want[2]),
            ("corrected.coverage", row.corrected_coverage, &want[3]),
        ] {
            let Ok(want_value) = text.parse::<f64>() else { continue };
            compared += 1;
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
        "percent.gc: {} of {} disagree as text; numeric: {exact}/{compared} bit-exact, \
         worst relative {worst_rel:.3e} ({worst_where})",
        gc_text_mismatches,
        rows.len()
    );

    assert_eq!(gc_text_mismatches, 0, "R's rounding of percent.gc was not reproduced");
    // The residual is the BLAS difference inside the loess, carried through
    // k.gc into the corrected coverage. See plethora_compat::loess::blas.
    assert!(
        worst_rel < 1e-9,
        "drifted {worst_rel:.3e} from the real upstream output: {worst_where}"
    );
}
