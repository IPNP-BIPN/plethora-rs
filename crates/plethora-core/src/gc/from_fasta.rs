//! `gc_from_fasta.pl`: percent GC per sequence.
//!
//! ```text
//! code/gc_from_fasta.pl $result.fa > ${result}_GC.txt
//! ```
//!
//! The script is four lines of Perl and two of them are surprising:
//!
//! ```text
//! $genome_size{$chr} += length($line);
//! $line =~ s/[ATN]//g;
//! $gc{$chr} += length($line);
//! ```
//!
//! What counts as GC is "everything that is not an uppercase A, T or N". So a
//! soft-masked repeat, which `bedtools getfasta` writes in lowercase, counts as
//! GC in full: `acgt` contributes four GC bases, not two. So do the IUPAC
//! ambiguity codes, R, Y, K, M and the rest. The denominator is the untouched
//! sequence length.
//!
//! This is almost certainly not what was meant, and it is what produced the
//! published GC model, so it is what this reproduces. The `GcCounts` type
//! reports the masked and ambiguous bases separately so a caller can see how
//! much of a given figure comes from that reading.

use std::io::{self, BufRead};

use plethora_compat::awk::format_g;

/// Perl's precision when it stringifies a double.
///
/// `print 1/3` gives `0.333333333333333` and `print 1e20` gives `1e+20`, which
/// is C's `%g` at precision 15. Note this is not R's rule: Perl writes
/// `0.0001` and `100000` where R's `as.character` writes `1e-04` and `1e+05`.
const PERL_PRECISION: usize = 15;

/// GC content of one sequence, with the parts that make it debatable split out.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct GcCounts {
    /// Total bases, the denominator.
    pub length: usize,
    /// Bases counted as GC, which is everything but uppercase A, T and N.
    pub gc: usize,
    /// Of those, how many are lowercase. Soft-masked repeat sequence.
    pub soft_masked: usize,
    /// Of those, how many are neither A, C, G, T nor N in either case.
    pub ambiguous: usize,
}

impl GcCounts {
    /// Adds one line of sequence.
    pub fn push_line(&mut self, line: &[u8]) {
        self.length += line.len();
        for &b in line {
            if matches!(b, b'A' | b'T' | b'N') {
                continue;
            }
            self.gc += 1;
            if b.is_ascii_lowercase() {
                self.soft_masked += 1;
            }
            if !matches!(b, b'C' | b'G' | b'c' | b'g' | b'a' | b't' | b'n') {
                self.ambiguous += 1;
            }
        }
    }

    /// The fraction the script prints.
    ///
    /// Zero-length sequences divide by zero in Perl and produce an "Illegal
    /// division by zero" death; here they yield NaN, which the caller can
    /// reject with a better message.
    #[must_use]
    pub fn fraction(&self) -> f64 {
        self.gc as f64 / self.length as f64
    }
}

/// One output row: a sequence name and its GC fraction.
#[derive(Debug, Clone, PartialEq)]
pub struct GcRow {
    pub name: String,
    pub counts: GcCounts,
}

impl GcRow {
    /// The line as the Perl prints it, with Perl's own number formatting.
    #[must_use]
    pub fn to_line(&self) -> String {
        format!("{}\t{}", self.name, format_g(self.counts.fraction(), PERL_PRECISION))
    }
}

/// Reads a FASTA and returns GC content per sequence.
///
/// Order follows the file. The Perl iterates a hash, so its output order is
/// whatever Perl's hash randomisation gives that run: a divergence in ordering
/// only, and one that makes the upstream file non-reproducible run to run.
///
/// A repeated sequence name accumulates into one entry, as the Perl's hash
/// does, rather than producing two rows.
///
/// # Errors
/// Returns an error if the input cannot be read.
pub fn gc_from_fasta<R: BufRead>(input: R) -> io::Result<Vec<GcRow>> {
    let mut rows: Vec<GcRow> = Vec::new();
    let mut index: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut current: Option<usize> = None;

    for line in input.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            // The whole header line is the key, `>` removed. bedtools writes
            // just the name under `-name` in 2.17, but `name::chrom:start-end`
            // from 2.27 on; see `build_model`.
            let name = header.to_string();
            let slot = *index.entry(name.clone()).or_insert_with(|| {
                rows.push(GcRow { name, counts: GcCounts::default() });
                rows.len() - 1
            });
            current = Some(slot);
            continue;
        }

        // Sequence before any header is dropped, as the Perl drops it: it would
        // accumulate under an undefined key.
        if let Some(slot) = current {
            rows[slot].counts.push_line(line.as_bytes());
        }
    }

    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts(sequence: &str) -> GcCounts {
        let mut c = GcCounts::default();
        c.push_line(sequence.as_bytes());
        c
    }

    #[test]
    fn plain_uppercase_sequence_counts_normally() {
        let c = counts("ACGT");
        assert_eq!(c.length, 4);
        assert_eq!(c.gc, 2);
        assert_eq!(c.fraction(), 0.5);
    }

    /// The behaviour worth knowing about: lowercase counts as GC in full.
    #[test]
    fn soft_masked_bases_all_count_as_gc() {
        let c = counts("acgt");
        assert_eq!(c.gc, 4, "every lowercase base counts, not just c and g");
        assert_eq!(c.soft_masked, 4);
        assert_eq!(c.fraction(), 1.0);
    }

    /// So do the ambiguity codes.
    #[test]
    fn ambiguity_codes_count_as_gc() {
        let c = counts("RYKM");
        assert_eq!(c.gc, 4);
        assert_eq!(c.ambiguous, 4);
    }

    /// N is excluded from the numerator but stays in the denominator, so a run
    /// of Ns lowers the reported GC rather than being skipped.
    #[test]
    fn uppercase_n_lowers_the_fraction() {
        let c = counts("GCNN");
        assert_eq!(c.length, 4);
        assert_eq!(c.gc, 2);
        assert_eq!(c.fraction(), 0.5);
    }

    /// Lowercase n, unlike uppercase N, counts as GC. The substitution in the
    /// Perl is case-sensitive.
    #[test]
    fn lowercase_n_counts_but_uppercase_does_not() {
        assert_eq!(counts("nnnn").gc, 4);
        assert_eq!(counts("NNNN").gc, 0);
    }

    /// The exact figures a mixed sequence gives, checked against the Perl's
    /// arithmetic by hand.
    #[test]
    fn a_mixed_sequence() {
        let c = counts("ACGTacgtNNNN");
        assert_eq!(c.length, 12);
        // C, G kept; a, c, g, t kept; A, T, and the four N dropped.
        assert_eq!(c.gc, 6);
        assert_eq!(c.fraction(), 0.5);
    }

    #[test]
    fn sequences_are_read_in_file_order() {
        let fasta = ">domB\nGGGG\n>domA\nAAAA\n";
        let rows = gc_from_fasta(fasta.as_bytes()).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].name, "domB");
        assert_eq!(rows[1].name, "domA");
    }

    /// Multi-line sequences accumulate, as the Perl accumulates per line.
    #[test]
    fn wrapped_sequences_accumulate() {
        let rows = gc_from_fasta(">d\nGG\nCC\nAA\n".as_bytes()).unwrap();
        assert_eq!(rows[0].counts.length, 6);
        assert_eq!(rows[0].counts.gc, 4);
    }

    /// A repeated name folds into one entry, as a Perl hash key would.
    #[test]
    fn a_repeated_name_accumulates_into_one_row() {
        let rows = gc_from_fasta(">d\nGG\n>d\nAA\n".as_bytes()).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].counts.length, 4);
        assert_eq!(rows[0].counts.gc, 2);
    }

    /// Perl's number formatting, which is %.15g and not R's rule.
    #[test]
    fn numbers_are_spelled_the_way_perl_spells_them() {
        let row = GcRow { name: "d".into(), counts: counts("GCA") };
        // 2/3 at fifteen significant digits.
        assert_eq!(row.to_line(), "d\t0.666666666666667");

        let row = GcRow { name: "d".into(), counts: counts("GCAT") };
        assert_eq!(row.to_line(), "d\t0.5");
    }

    /// A 1000 bp domain gives a three-decimal figure, which is why the real GC
    /// file is full of values that land exactly on a rounding boundary.
    #[test]
    fn a_thousand_base_domain_gives_three_decimals() {
        let mut c = GcCounts::default();
        c.push_line(&b"G".repeat(348));
        c.push_line(&b"A".repeat(652));
        let row = GcRow { name: "baseline_1_1".into(), counts: c };
        assert_eq!(row.to_line(), "baseline_1_1\t0.348");
    }
}
