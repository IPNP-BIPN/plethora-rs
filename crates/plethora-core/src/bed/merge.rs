//! `bedtools merge -c 5 -o sum`, on a file whose first column is not a
//! chromosome.
//!
//! ```text
//! awk 'OFS="\t" {print $4,$2,$3,$1,$13}' ${output}_temp.bed | bedtools merge -c 5 -o sum -i -
//! ```
//!
//! The awk in front is doing something worth spelling out. It permutes the
//! thirteen intersect columns into `domain name, domain start, domain end, real
//! chromosome, overlap`, so the *domain name* takes the place bedtools reads as
//! a chromosome. Merging then groups by domain rather than by locus, and the sum
//! it produces is the total number of overlapping bases per domain.
//!
//! That works only because every row for one domain is contiguous in the file,
//! which it is: `intersect` emits A's intervals in order and all the hits for
//! one interval together. bedtools checks, and a domain name that reappears
//! after another one has intervened is a hard error, not a silent regroup. Since
//! the DUF1220 annotation carries a name per domain instance, that error would
//! mean the reference file has a duplicated name.

use std::io::Write;

/// A row of the permuted file: an interval keyed by name, carrying a value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Row {
    /// Column one, which bedtools treats as the chromosome. Here it is the
    /// domain name.
    pub key: String,
    pub start: i64,
    pub end: i64,
    /// The column named by `-c`, summed by `-o sum`.
    pub value: i64,
}

impl Row {
    /// Parses a permuted row: key, start, end, then the summed column last.
    #[must_use]
    pub fn parse(line: &str, value_column: usize) -> Option<Self> {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < value_column {
            return None;
        }
        Some(Self {
            key: f[0].to_string(),
            start: f[1].parse().ok()?,
            end: f[2].parse().ok()?,
            value: f[value_column - 1].parse().ok()?,
        })
    }
}

/// Input whose keys are not grouped, which bedtools refuses.
#[derive(Debug, thiserror::Error)]
#[error("sorted input specified, but the file has the following out of order record\n{record}")]
pub struct OutOfOrder {
    pub record: String,
}

/// Merges overlapping or book-ended rows sharing a key, summing the value.
///
/// Book-ended rows do merge: bedtools' default distance is zero, which means
/// touching counts. In this pipeline every row for a domain carries that
/// domain's own coordinates, so they are identical rather than merely touching
/// and collapse to one.
///
/// # Errors
/// Returns an error if a key reappears after another key intervened, or if
/// writing fails.
pub fn merge_sum<I, W>(
    rows: I,
    value_column: usize,
    mut out: W,
) -> Result<(), Box<dyn std::error::Error>>
where
    I: Iterator<Item = String>,
    W: Write,
{
    let mut seen: Vec<String> = Vec::new();
    let mut current: Option<Row> = None;

    for line in rows {
        let Some(row) = Row::parse(&line, value_column) else {
            continue;
        };

        match &mut current {
            Some(open) if open.key == row.key => {
                if row.start < open.start {
                    return Err(Box::new(OutOfOrder { record: line }));
                }
                if row.start <= open.end {
                    // Overlapping or book-ended: extend and accumulate.
                    open.end = open.end.max(row.end);
                    open.value += row.value;
                } else {
                    writeln!(
                        out,
                        "{}\t{}\t{}\t{}",
                        open.key, open.start, open.end, open.value
                    )?;
                    *open = row;
                }
            }
            _ => {
                if let Some(open) = current.take() {
                    writeln!(
                        out,
                        "{}\t{}\t{}\t{}",
                        open.key, open.start, open.end, open.value
                    )?;
                }
                if seen.contains(&row.key) {
                    return Err(Box::new(OutOfOrder { record: line }));
                }
                seen.push(row.key.clone());
                current = Some(row);
            }
        }
    }

    if let Some(open) = current {
        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            open.key, open.start, open.end, open.value
        )?;
    }

    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(rows: &[&str]) -> Result<Vec<String>, String> {
        let mut out = Vec::new();
        merge_sum(rows.iter().map(|s| (*s).to_string()), 5, &mut out).map_err(|e| e.to_string())?;
        Ok(String::from_utf8(out)
            .unwrap()
            .lines()
            .map(String::from)
            .collect())
    }

    /// The shape this stage actually meets: identical intervals, one per read
    /// that hit the domain.
    #[test]
    fn identical_intervals_collapse_and_sum() {
        let out = run(&[
            "domA\t100\t200\tchr1\t50",
            "domA\t100\t200\tchr1\t10",
            "domB\t300\t400\tchr1\t20",
        ])
        .unwrap();
        assert_eq!(out, ["domA\t100\t200\t60", "domB\t300\t400\t20"]);
    }

    /// A domain with no coverage still produces a row, summing to zero.
    #[test]
    fn an_uncovered_domain_sums_to_zero() {
        let out = run(&["domC\t500\t600\tchr1\t0"]).unwrap();
        assert_eq!(out, ["domC\t500\t600\t0"]);
    }

    /// Same key, disjoint intervals: bedtools keeps them apart rather than
    /// merging everything under one name.
    #[test]
    fn disjoint_intervals_under_one_key_stay_apart() {
        let out = run(&["domA\t100\t200\tchr1\t50", "domA\t700\t800\tchr1\t7"]).unwrap();
        assert_eq!(out, ["domA\t100\t200\t50", "domA\t700\t800\t7"]);
    }

    /// Touching intervals merge, since the default distance is zero.
    #[test]
    fn book_ended_intervals_merge() {
        let out = run(&["domA\t100\t200\tchr1\t5", "domA\t200\t300\tchr1\t7"]).unwrap();
        assert_eq!(out, ["domA\t100\t300\t12"]);
    }

    /// The error bedtools raises, reproduced: a name that comes back after
    /// another one means the reference has a duplicated domain name.
    #[test]
    fn a_key_that_reappears_is_an_error() {
        let err = run(&[
            "domA\t100\t200\tchr1\t50",
            "domB\t300\t400\tchr1\t20",
            "domA\t100\t200\tchr1\t10",
        ])
        .unwrap_err();
        assert!(err.contains("out of order"), "unexpected message: {err}");
    }

    #[test]
    fn a_backwards_start_within_one_key_is_an_error() {
        let err = run(&["domA\t300\t400\tchr1\t5", "domA\t100\t200\tchr1\t5"]).unwrap_err();
        assert!(err.contains("out of order"), "unexpected message: {err}");
    }

    #[test]
    fn an_empty_input_produces_nothing() {
        assert!(run(&[]).unwrap().is_empty());
    }
}
