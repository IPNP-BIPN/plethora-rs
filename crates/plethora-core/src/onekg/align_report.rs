//! `zip.sh`: count the aligned fragments, then compress the BED.
//!
//! ```text
//! bed_file=results/${sample}_sorted.bed
//! n_align=`cut -f 4 $bed_file | sort -T ./ | uniq | wc -l`
//! echo $sample $n_align >> align_report.txt
//! gzip $bed_file
//! ```
//!
//! `cut | sort | uniq | wc -l` counts *distinct* names, not lines: a fragment
//! that overlaps several domains appears once. `trim_qc_report.R` then divides
//! that by the number of reads surviving trimming to get a percent aligned.
//!
//! Two things about the shell version are worth knowing. The count is a full
//! sort of a whole-genome BED's fourth column, which is why upstream points it
//! at the working directory with `-T ./`; here it is a hash set, which needs no
//! spill. And the append is done by every job in a 300-wide array at once,
//! which is a race: see [`append_line`].

use std::collections::HashSet;
use std::io::{BufRead, Write};
use std::path::Path;

/// One line of the report: a sample and its distinct aligned fragments.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Entry {
    pub sample: String,
    pub aligned_fragments: usize,
}

impl Entry {
    /// The line as `echo $sample $n_align` writes it: space-separated, which is
    /// why `trim_qc_report.R` reads the file with `sep = " "`.
    #[must_use]
    pub fn to_line(&self) -> String {
        format!("{} {}", self.sample, self.aligned_fragments)
    }

    /// Parses a line back.
    ///
    /// `zip.sh` writes the file name rather than the sample, so a report
    /// written by upstream reads `HG00250_sorted.bed 55000000`. The R strips
    /// that suffix on the way in (`gsub("_sorted.bed", "", SAMPLE_NAME)`) and
    /// so does this, which is what lets one report hold lines from both.
    #[must_use]
    pub fn parse(line: &str) -> Option<Self> {
        let (sample, count) = line.split_once(' ')?;
        let sample = sample.strip_suffix("_sorted.bed").unwrap_or(sample);
        Some(Self {
            sample: sample.to_string(),
            aligned_fragments: count.trim().parse().ok()?,
        })
    }
}

/// Counts distinct names in the fourth column.
///
/// A `HashSet` rather than a sort: the answer is a cardinality, and the sort
/// upstream performs exists only because `uniq` needs its input grouped.
///
/// # Errors
/// Returns an error if the input cannot be read.
pub fn count_distinct_names<R: BufRead>(input: R) -> std::io::Result<usize> {
    let mut seen: HashSet<String> = HashSet::new();
    for line in input.lines() {
        let line = line?;
        // `cut -f 4` yields an empty field for a short line rather than
        // skipping it, and `uniq` then counts that empty string once.
        let name = line.split('\t').nth(3).unwrap_or("");
        if !seen.contains(name) {
            seen.insert(name.to_string());
        }
    }
    Ok(seen.len())
}

/// Counts a sample's aligned fragments from its sorted BED.
///
/// # Errors
/// Returns an error if the file cannot be read.
pub fn count_for_sample(sorted_bed: &Path) -> std::io::Result<usize> {
    count_distinct_names(crate::io::open(sorted_bed)?)
}

/// Appends one entry to the report, atomically enough for a job array.
///
/// Upstream appends with `>>` from every job in a 300-wide LSF array. Two jobs
/// finishing together can interleave inside one line, and the result is a
/// corrupt row that `trim_qc_report.R` reads as a missing sample. Opening in
/// append mode and issuing the line as a single `write_all` keeps it whole:
/// POSIX guarantees an atomic append for a write below `PIPE_BUF`, and a line
/// here is a sample name and an integer.
///
/// # Errors
/// Returns an error if the file cannot be opened or written.
pub fn append_line(report: &Path, entry: &Entry) -> std::io::Result<()> {
    let mut file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(report)?;
    let line = format!("{}\n", entry.to_line());
    file.write_all(line.as_bytes())?;
    file.flush()
}

/// Reads a report back, skipping rows that did not survive being written.
///
/// # Errors
/// Returns an error if the file cannot be read.
pub fn read_report(path: &Path) -> std::io::Result<Vec<Entry>> {
    let mut out = Vec::new();
    for line in crate::io::open(path)?.lines() {
        if let Some(entry) = Entry::parse(&line?) {
            out.push(entry);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distinct_names_are_counted_not_lines() {
        // Three lines, two names: a fragment overlapping two domains.
        let bed = "chr1\t1\t2\tfragA\t60\t+\nchr1\t3\t4\tfragA\t60\t+\nchr1\t5\t6\tfragB\t60\t-\n";
        assert_eq!(count_distinct_names(bed.as_bytes()).unwrap(), 2);
    }

    #[test]
    fn an_empty_file_counts_nothing() {
        assert_eq!(count_distinct_names("".as_bytes()).unwrap(), 0);
    }

    /// `cut -f 4` on a short line gives an empty field, and `uniq` counts it,
    /// so a malformed BED inflates the count by exactly one.
    #[test]
    fn short_lines_contribute_one_empty_name() {
        let bed = "chr1\t1\t2\tfragA\t60\t+\nchr1\t3\t4\nchr1\t5\t6\n";
        assert_eq!(count_distinct_names(bed.as_bytes()).unwrap(), 2);
    }

    #[test]
    fn the_line_is_space_separated() {
        let e = Entry {
            sample: "HG00250".into(),
            aligned_fragments: 123_456,
        };
        assert_eq!(e.to_line(), "HG00250 123456");
        assert_eq!(Entry::parse("HG00250 123456"), Some(e));
    }

    #[test]
    fn a_malformed_row_is_skipped_rather_than_guessed() {
        assert_eq!(Entry::parse("HG00250"), None);
        assert_eq!(Entry::parse("HG00250 notanumber"), None);
        assert_eq!(Entry::parse(""), None);
    }

    #[test]
    fn entries_append_and_read_back() {
        let dir = tempfile::tempdir().unwrap();
        let report = dir.path().join("align_report.txt");
        for (sample, n) in [("HG00250", 10), ("HG00251", 20)] {
            append_line(
                &report,
                &Entry {
                    sample: sample.into(),
                    aligned_fragments: n,
                },
            )
            .unwrap();
        }
        let entries = read_report(&report).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[1].sample, "HG00251");
        assert_eq!(entries[1].aligned_fragments, 20);
    }

    /// The case the shell version loses: many writers appending at once. Each
    /// line must arrive whole, whatever order they land in.
    #[test]
    fn concurrent_appends_do_not_interleave() {
        let dir = tempfile::tempdir().unwrap();
        let report = dir.path().join("align_report.txt");

        std::thread::scope(|scope| {
            for i in 0..32 {
                let report = report.clone();
                scope.spawn(move || {
                    for j in 0..32 {
                        append_line(
                            &report,
                            &Entry {
                                sample: format!("S{i:02}_{j:02}"),
                                aligned_fragments: i * 1000 + j,
                            },
                        )
                        .unwrap();
                    }
                });
            }
        });

        let entries = read_report(&report).unwrap();
        assert_eq!(entries.len(), 32 * 32, "a line was lost or split");
        // Every line parsed, so none was cut in half.
        let text = std::fs::read_to_string(&report).unwrap();
        assert_eq!(text.lines().count(), 32 * 32);
        for line in text.lines() {
            assert!(Entry::parse(line).is_some(), "corrupt line: {line:?}");
        }
    }

    /// A gzipped BED is read transparently, which is what `zip.sh` leaves
    /// behind after it runs.
    #[test]
    fn a_gzipped_bed_counts_the_same() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("s_sorted.bed.gz");
        let mut w = crate::io::create(&path).unwrap();
        w.write_all(b"chr1\t1\t2\tfragA\t60\t+\nchr1\t3\t4\tfragB\t60\t+\n")
            .unwrap();
        drop(w);
        assert_eq!(count_for_sample(&path).unwrap(), 2);
    }
}
