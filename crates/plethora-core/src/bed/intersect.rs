//! `bedtools intersect -wao -sorted`.
//!
//! ```text
//! bedtools intersect -wao -sorted -a $reference_bed -b ${output}_sorted.bed > ${output}_temp.bed
//! ```
//!
//! A is the reference: one interval per domain of interest. B is the sample's
//! reads. Every A interval produces at least one output line, so a domain with
//! no coverage is still represented, with a null B and an overlap of zero. That
//! is what `-wao` means and it is why the downstream sum has a row for every
//! domain rather than only the covered ones.
//!
//! Output is A's six columns, B's six columns, and the overlap: thirteen fields.
//! A domain with no overlapping read gets `. -1 -1 . -1 .` and `0`.

use std::io::Write;

/// A BED6 interval, kept as parsed coordinates plus the original trailing text.
///
/// The tail is carried verbatim rather than reformatted, so name, score and
/// strand reach the output exactly as they arrived.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    /// Columns four onwards, tab-separated, exactly as read.
    pub rest: String,
}

impl Interval {
    /// Parses a BED line. Returns `None` for a line with fewer than three
    /// fields, which is what bedtools rejects as malformed.
    #[must_use]
    pub fn parse(line: &str) -> Option<Self> {
        let mut it = line.splitn(4, '\t');
        let chrom = it.next()?;
        let start = it.next()?.parse().ok()?;
        let end = it.next()?.parse().ok()?;
        let rest = it.next().unwrap_or("");
        Some(Self {
            chrom: chrom.to_string(),
            start,
            end,
            rest: rest.to_string(),
        })
    }

    /// The line as bedtools would echo it.
    #[must_use]
    pub fn to_line(&self) -> String {
        if self.rest.is_empty() {
            format!("{}\t{}\t{}", self.chrom, self.start, self.end)
        } else {
            format!("{}\t{}\t{}\t{}", self.chrom, self.start, self.end, self.rest)
        }
    }

    /// Overlap in base pairs, zero when the intervals only touch.
    ///
    /// bedtools requires at least one base by default, so book-ended intervals
    /// do not count as intersecting.
    #[must_use]
    pub fn overlap(&self, other: &Self) -> i64 {
        if self.chrom != other.chrom {
            return 0;
        }
        (self.end.min(other.end) - self.start.max(other.start)).max(0)
    }
}

/// The null B record bedtools writes for an A interval with no overlap.
///
/// Six fields shaped to BED6: the string columns become `.` and the numeric
/// ones `-1`. bedtools derives this from B's column count; the pipeline always
/// feeds it BED6, so the shape is fixed here.
pub const NULL_B: &str = ".\t-1\t-1\t.\t-1\t.";

/// An input whose chromosomes are not in a consistent order.
#[derive(Debug, thiserror::Error)]
#[error("sorted input specified, but {file} has the following out of order record\n{record}")]
pub struct OutOfOrder {
    pub file: &'static str,
    pub record: String,
}

/// Streams `-wao` output for A against B, both sorted by chromosome then start.
///
/// The sweep keeps only the B intervals that can still overlap something: once
/// a B interval ends at or before the current A interval's start, no later A
/// can reach it, because A's starts are non-decreasing.
///
/// # Errors
/// Returns an error if either input is out of order, or if writing fails.
///
/// # Panics
/// Panics if the peeked B interval vanishes before it is consumed, which the
/// iterator contract forbids.
pub fn intersect_wao<A, B, W>(a: A, b: B, mut out: W) -> Result<(), Box<dyn std::error::Error>>
where
    A: Iterator<Item = String>,
    B: Iterator<Item = String>,
    W: Write,
{
    let mut b = b.filter_map(|l| Interval::parse(&l)).peekable();
    let mut cache: Vec<Interval> = Vec::new();

    let mut last_a: Option<(String, i64)> = None;
    let mut b_chroms_seen: Vec<String> = Vec::new();

    for line in a {
        let Some(feature) = Interval::parse(&line) else {
            continue;
        };

        if let Some((chrom, start)) = &last_a {
            let ordered = match chrom.as_str().cmp(feature.chrom.as_str()) {
                std::cmp::Ordering::Less => true,
                std::cmp::Ordering::Equal => *start <= feature.start,
                std::cmp::Ordering::Greater => false,
            };
            if !ordered {
                return Err(Box::new(OutOfOrder {
                    file: "the -a file",
                    record: feature.to_line(),
                }));
            }
        }
        last_a = Some((feature.chrom.clone(), feature.start));

        // Drop what can no longer overlap: a different chromosome, or an end at
        // or before this interval's start.
        cache.retain(|c| c.chrom == feature.chrom && c.end > feature.start);

        // Pull in every B interval that starts before this one ends.
        while let Some(next) = b.peek() {
            match next.chrom.as_str().cmp(feature.chrom.as_str()) {
                // B is still on an earlier chromosome: consume and discard.
                std::cmp::Ordering::Less => {
                    let consumed = b.next().expect("peeked");
                    if !b_chroms_seen.contains(&consumed.chrom) {
                        b_chroms_seen.push(consumed.chrom);
                    }
                }
                std::cmp::Ordering::Equal => {
                    if next.start < feature.end {
                        let consumed = b.next().expect("peeked");
                        if !b_chroms_seen.contains(&consumed.chrom) {
                            b_chroms_seen.push(consumed.chrom.clone());
                        }
                        if consumed.end > feature.start {
                            cache.push(consumed);
                        }
                    } else {
                        break;
                    }
                }
                // B has moved past this chromosome.
                std::cmp::Ordering::Greater => break,
            }
        }

        let mut hits = 0;
        for candidate in &cache {
            let overlap = feature.overlap(candidate);
            if overlap > 0 {
                writeln!(out, "{}\t{}\t{overlap}", feature.to_line(), candidate.to_line())?;
                hits += 1;
            }
        }
        if hits == 0 {
            writeln!(out, "{}\t{NULL_B}\t0", feature.to_line())?;
        }
    }

    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(a: &[&str], b: &[&str]) -> Vec<String> {
        let mut out = Vec::new();
        intersect_wao(
            a.iter().map(|s| (*s).to_string()),
            b.iter().map(|s| (*s).to_string()),
            &mut out,
        )
        .expect("intersect");
        String::from_utf8(out).unwrap().lines().map(String::from).collect()
    }

    #[test]
    fn reports_every_overlap_separately() {
        let out = run(
            &["chr1\t100\t200\tdomA\t255\t+"],
            &["chr1\t150\t250\tr1\t60\t+", "chr1\t160\t170\tr2\t30\t-"],
        );
        assert_eq!(
            out,
            [
                "chr1\t100\t200\tdomA\t255\t+\tchr1\t150\t250\tr1\t60\t+\t50",
                "chr1\t100\t200\tdomA\t255\t+\tchr1\t160\t170\tr2\t30\t-\t10",
            ]
        );
    }

    /// The point of -wao: an uncovered domain still produces a row, so the sum
    /// downstream has an entry for it.
    #[test]
    fn an_uncovered_interval_still_produces_a_row() {
        let out = run(&["chr1\t500\t600\tdomC\t255\t-"], &["chr1\t100\t200\tr1\t60\t+"]);
        assert_eq!(out, ["chr1\t500\t600\tdomC\t255\t-\t.\t-1\t-1\t.\t-1\t.\t0"]);
    }

    /// Book-ended intervals share no base, so bedtools does not count them.
    #[test]
    fn touching_intervals_do_not_intersect() {
        let out = run(&["chr1\t100\t200\tdomA\t255\t+"], &["chr1\t200\t300\tr1\t60\t+"]);
        assert!(out[0].ends_with("\t.\t-1\t-1\t.\t-1\t.\t0"));
    }

    #[test]
    fn different_chromosomes_never_overlap() {
        let out = run(&["chr1\t100\t200\tdomA\t255\t+"], &["chr2\t100\t200\tr1\t60\t+"]);
        assert!(out[0].ends_with("\t0"));
        assert_eq!(out.len(), 1);
    }

    /// Overlapping reference intervals must both see a read that spans them,
    /// which is what the retained cache is for.
    #[test]
    fn a_read_can_serve_several_overlapping_domains() {
        let out = run(
            &["chr1\t100\t300\tdomA\t255\t+", "chr1\t200\t400\tdomB\t255\t+"],
            &["chr1\t150\t350\tr1\t60\t+"],
        );
        assert_eq!(out.len(), 2);
        assert!(out[0].ends_with("\t150"), "domA sees 150 bp: {}", out[0]);
        assert!(out[1].ends_with("\t150"), "domB sees 150 bp: {}", out[1]);
    }

    #[test]
    fn several_chromosomes_stream_in_order() {
        let out = run(
            &[
                "chr1\t100\t200\tdomA\t255\t+",
                "chr2\t100\t200\tdomD\t255\t+",
            ],
            &["chr1\t150\t250\tr1\t60\t+", "chr2\t150\t160\tr4\t10\t+"],
        );
        assert_eq!(out.len(), 2);
        assert!(out[0].ends_with("\t50"));
        assert!(out[1].ends_with("\t10"));
    }

    #[test]
    fn an_out_of_order_reference_is_rejected() {
        let mut out = Vec::new();
        let result = intersect_wao(
            ["chr1\t500\t600\tx\t0\t+", "chr1\t100\t200\ty\t0\t+"]
                .into_iter()
                .map(String::from),
            std::iter::empty(),
            &mut out,
        );
        assert!(result.is_err(), "an unsorted -a file must be refused");
    }

    #[test]
    fn overlap_is_the_shared_span() {
        let a = Interval::parse("chr1\t100\t200\tx\t0\t+").unwrap();
        let b = Interval::parse("chr1\t150\t250\ty\t0\t+").unwrap();
        assert_eq!(a.overlap(&b), 50);
        assert_eq!(b.overlap(&a), 50);
        let c = Interval::parse("chr1\t250\t300\tz\t0\t+").unwrap();
        assert_eq!(a.overlap(&c), 0);
    }
}
