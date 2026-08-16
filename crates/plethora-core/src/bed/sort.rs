//! `sort -k 1,1 -k 2,2n -T ./` over BED lines.
//!
//! ```text
//! sort -k 1,1 -k 2,2n -T ./ ${output}_edited.bed > ${output}_sorted.bed
//! ```
//!
//! The comparison is [`plethora_compat::gnusort::cmp_k1_k2n`], including the
//! whole-line tie-break that orders most of a read-interval file. This module
//! is the sort around it, with the same disk spilling `-T ./` implies: a
//! whole-genome interval file does not fit in memory.

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use rayon::slice::ParallelSliceMut as _;

use plethora_compat::gnusort::cmp_k1_k2n;

/// Lines held in memory before a run is written out.
///
/// About 100 MB of BED text, comfortably below what `sort` allocates by
/// default. Configurable so the tests can force spilling.
pub const DEFAULT_RUN_LINES: usize = 2_000_000;

/// Sorts one run, across every core.
///
/// Unstable, and parallel, both of which are safe here for the same reason:
/// [`cmp_k1_k2n`] falls back to comparing the whole line, so two lines compare
/// `Equal` only when they are byte-identical. A reordering among equal elements
/// is therefore invisible in the output, and the result is the same total order
/// however it was reached.
///
/// Below the threshold the sequential sort wins: rayon's split and join cost
/// more than sorting a few thousand lines.
fn sort_run(buffer: &mut [String]) {
    const PARALLEL_THRESHOLD: usize = 8192;

    if buffer.len() < PARALLEL_THRESHOLD {
        buffer.sort_unstable_by(|a, b| cmp_k1_k2n(a.as_bytes(), b.as_bytes()));
    } else {
        buffer.par_sort_unstable_by(|a, b| cmp_k1_k2n(a.as_bytes(), b.as_bytes()));
    }
}

/// Sorts lines, spilling runs into `tmp_dir` when they grow past `run_lines`.
///
/// GNU `sort` is not stable and breaks ties by comparing whole lines, so the
/// comparison is total and stability is not needed for a faithful result.
///
/// # Errors
/// Returns an error if a run file cannot be written or read back, or if the
/// output cannot be written.
///
/// # Panics
/// Panics if `run_lines` is zero, which would never make progress.
pub fn sort_lines<I, W>(lines: I, run_lines: usize, tmp_dir: &Path, mut out: W) -> io::Result<()>
where
    I: Iterator<Item = String>,
    W: Write,
{
    assert!(run_lines > 0, "a sorting run must hold at least one line");

    // Removed when this goes out of scope rather than after a successful
    // merge, so a failure part way through does not leave them behind.
    let mut runs = Runs::default();
    let mut buffer: Vec<String> = Vec::new();

    let flush_run = |buffer: &mut Vec<String>, runs: &mut Runs| -> io::Result<()> {
        sort_run(buffer);
        let path = tmp_dir.join(format!("bedsort-{}.run", runs.paths.len()));
        let mut w = BufWriter::new(File::create(&path)?);
        for line in buffer.iter() {
            writeln!(w, "{line}")?;
        }
        w.flush()?;
        runs.paths.push(path);
        buffer.clear();
        Ok(())
    };

    for line in lines {
        buffer.push(line);
        if buffer.len() >= run_lines {
            flush_run(&mut buffer, &mut runs)?;
        }
    }

    if runs.paths.is_empty() {
        sort_run(&mut buffer);
        for line in buffer {
            writeln!(out, "{line}")?;
        }
        return out.flush();
    }

    if !buffer.is_empty() {
        flush_run(&mut buffer, &mut runs)?;
    }

    let mut readers: Vec<_> = runs
        .paths
        .iter()
        .map(|p| File::open(p).map(|f| BufReader::new(f).lines()))
        .collect::<io::Result<_>>()?;
    let mut heads: Vec<Option<String>> = readers
        .iter_mut()
        .map(|r| r.next().transpose())
        .collect::<io::Result<_>>()?;

    loop {
        let mut best: Option<usize> = None;
        for (i, head) in heads.iter().enumerate() {
            let Some(candidate) = head else { continue };
            match best {
                None => best = Some(i),
                Some(b) => {
                    let current = heads[b].as_ref().expect("best index holds a line");
                    if cmp_k1_k2n(candidate.as_bytes(), current.as_bytes())
                        == std::cmp::Ordering::Less
                    {
                        best = Some(i);
                    }
                }
            }
        }

        let Some(i) = best else { break };
        writeln!(out, "{}", heads[i].take().expect("chosen run holds a line"))?;
        heads[i] = readers[i].next().transpose()?;
    }

    out.flush()
}

/// The spilled runs, removed when the sort ends however it ends.
#[derive(Default)]
struct Runs {
    paths: Vec<PathBuf>,
}

impl Drop for Runs {
    fn drop(&mut self) {
        for path in &self.paths {
            // A failure to clean up is not a failure to sort.
            let _ = std::fs::remove_file(path);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A merge that fails part way must not leave its spilled runs behind.
    /// They are removed on drop rather than after a successful merge, which is
    /// the only version of this that holds when something goes wrong.
    #[test]
    fn a_failed_merge_leaves_no_runs() {
        struct Fails;
        impl Write for Fails {
            fn write(&mut self, _: &[u8]) -> io::Result<usize> {
                Err(io::Error::other("the sink gave up"))
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }

        let dir = tempfile::tempdir().expect("tempdir");
        let lines = (0..100).map(|i| format!("chr1\t{i}\t{}", i + 10));
        // Two lines a run, so it spills and has a merge to fail during.
        let err = sort_lines(lines, 2, dir.path(), Fails).expect_err("the sink fails");
        assert!(err.to_string().contains("gave up"), "got: {err}");

        let left: Vec<_> = std::fs::read_dir(dir.path())
            .expect("read dir")
            .filter_map(Result::ok)
            .map(|e| e.file_name().to_string_lossy().into_owned())
            .collect();
        assert!(left.is_empty(), "runs left behind: {left:?}");
    }

    fn run(lines: &[&str], run_lines: usize) -> Vec<String> {
        let dir = tempfile::tempdir().unwrap();
        let mut out = Vec::new();
        sort_lines(
            lines.iter().map(|s| (*s).to_string()),
            run_lines,
            dir.path(),
            &mut out,
        )
        .unwrap();
        String::from_utf8(out)
            .unwrap()
            .lines()
            .map(String::from)
            .collect()
    }

    #[test]
    fn sorts_by_chromosome_then_numeric_start() {
        let out = run(&["chr2\t100\ta", "chr1\t200\tb", "chr1\t30\tc"], 100);
        assert_eq!(out, ["chr1\t30\tc", "chr1\t200\tb", "chr2\t100\ta"]);
    }

    /// Chromosomes compare as bytes, so chr10 precedes chr9.
    #[test]
    fn chromosomes_compare_as_bytes() {
        let out = run(&["chr9\t1\ta", "chr10\t1\tb"], 100);
        assert_eq!(out, ["chr10\t1\tb", "chr9\t1\ta"]);
    }

    /// The tie-break that orders most of a real file.
    #[test]
    fn equal_keys_fall_back_to_the_whole_line() {
        let out = run(&["chr1\t100\tzzz", "chr1\t100\taaa"], 100);
        assert_eq!(out, ["chr1\t100\taaa", "chr1\t100\tzzz"]);
    }

    #[test]
    fn spilling_gives_the_same_answer() {
        let lines: Vec<String> = (0..200)
            .map(|i| format!("chr{}\t{}\tread{i}", i % 3, 1000 - i))
            .collect();
        let refs: Vec<&str> = lines.iter().map(String::as_str).collect();
        assert_eq!(run(&refs, 1000), run(&refs, 7));
        assert_eq!(run(&refs, 1000), run(&refs, 1));
    }

    #[test]
    fn an_empty_input_sorts_to_nothing() {
        assert!(run(&[], 10).is_empty());
    }
}
