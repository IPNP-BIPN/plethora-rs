//! `build_gc_model.sh`: measure GC per domain straight from the genome.
//!
//! ```text
//! bedtools getfasta -name -fi $genome -bed $bed -fo $result.fa
//! code/gc_from_fasta.pl $result.fa > ${result}_GC.txt
//! rm $result.fa
//! ```
//!
//! The FASTA in the middle exists only to be read back and deleted, so this
//! goes straight from the BED and the genome to the GC table. A
//! [`write_fasta`] is kept for interoperability with anything that still wants
//! the intermediate.
//!
//! One thing to know before regenerating a GC table: the naming that
//! `getfasta -name` produces changed. Up to bedtools 2.26 the header was the
//! BED name alone; from 2.27 it is `name::chrom:start-end`, which no longer
//! matches the domain names in the read-depth file, so the join in
//! `gc_correction.R` silently produces nothing. This module writes the bare
//! name, which is what the rest of the pipeline expects.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

use noodles_fasta::fai;

use super::from_fasta::{GcCounts, GcRow};

/// One interval to extract: a name and its coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    pub name: String,
    pub chrom: String,
    /// Zero-based, half-open, as in BED.
    pub start: u64,
    pub end: u64,
}

impl Region {
    /// Parses a BED line, taking column four as the name.
    ///
    /// A line without a name column is skipped: `getfasta -name` has nothing to
    /// call the sequence, and an unnamed row could not join with anything
    /// downstream.
    #[must_use]
    pub fn parse(line: &str) -> Option<Self> {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 4 {
            return None;
        }
        Some(Self {
            name: f[3].to_string(),
            chrom: f[0].to_string(),
            start: f[1].parse().ok()?,
            end: f[2].parse().ok()?,
        })
    }

    /// The header `getfasta` writes for this region.
    ///
    /// `Naming::Bare` is bedtools 2.26 and earlier, and what the pipeline
    /// joins on. `Naming::WithCoordinates` is 2.27 and later.
    #[must_use]
    pub fn header(&self, naming: Naming) -> String {
        match naming {
            Naming::Bare => self.name.clone(),
            Naming::WithCoordinates => {
                format!("{}::{}:{}-{}", self.name, self.chrom, self.start, self.end)
            }
        }
    }
}

/// Which of `getfasta -name`'s two header conventions to write.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Naming {
    /// Just the BED name, as bedtools 2.26 and earlier wrote it. The default,
    /// because it is the only one the rest of the pipeline can join on.
    #[default]
    Bare,
    /// `name::chrom:start-end`, as bedtools 2.27 and later write it.
    WithCoordinates,
}

/// A genome that regions can be pulled out of, by its FASTA index.
pub struct Genome {
    reader: noodles_fasta::io::IndexedReader<BufReader<File>>,
}

impl Genome {
    /// Opens a FASTA and its `.fai` index.
    ///
    /// The index must exist. bedtools builds one on demand; requiring it here
    /// keeps this from silently spending minutes indexing a human genome.
    ///
    /// # Errors
    /// Returns an error if either file cannot be opened or the index cannot be
    /// parsed.
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let path = path.as_ref();
        let index_path = path.with_extension(format!(
            "{}.fai",
            path.extension().map_or(String::new(), |e| e.to_string_lossy().into_owned())
        ));
        let index = fai::fs::read(&index_path).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("cannot read the FASTA index {}: {e}", index_path.display()),
            )
        })?;
        let inner = BufReader::new(File::open(path)?);
        Ok(Self {
            reader: noodles_fasta::io::IndexedReader::new(inner, index),
        })
    }

    /// Pulls one region's sequence out of the genome.
    ///
    /// Not reverse-complemented for a minus-strand feature, matching
    /// `getfasta` without `-s`, which is how `build_gc_model.sh` calls it. GC
    /// content is the same on both strands, so this cannot change the result;
    /// it is stated because a reader will wonder.
    ///
    /// # Errors
    /// Returns an error if the region is not in the index or cannot be read.
    pub fn sequence(&mut self, region: &Region) -> io::Result<Vec<u8>> {
        // noodles regions are one-based and inclusive; BED is zero-based and
        // half-open.
        let spec = format!("{}:{}-{}", region.chrom, region.start + 1, region.end);
        let parsed: noodles_core::Region = spec
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, format!("{spec}: {e}")))?;
        let record = self.reader.query(&parsed)?;
        Ok(record.sequence().as_ref().to_vec())
    }
}

/// Measures GC for every region in a BED file.
///
/// # Errors
/// Returns an error if the BED cannot be read or a region is not in the genome.
pub fn build<R: BufRead>(bed: R, genome: &mut Genome) -> io::Result<Vec<GcRow>> {
    let mut rows = Vec::new();
    for line in bed.lines() {
        let line = line?;
        let Some(region) = Region::parse(&line) else {
            continue;
        };
        let sequence = genome.sequence(&region)?;
        let mut counts = GcCounts::default();
        counts.push_line(&sequence);
        rows.push(GcRow {
            name: region.name,
            counts,
        });
    }
    Ok(rows)
}

/// Writes the intermediate FASTA `getfasta` would have written.
///
/// Not used by [`build`], which skips the round trip. Kept for anything that
/// still wants the file, and because it is where the two naming conventions
/// become visible.
///
/// # Errors
/// Returns an error if the BED cannot be read, a region is missing, or writing
/// fails.
pub fn write_fasta<R: BufRead, W: Write>(
    bed: R,
    genome: &mut Genome,
    naming: Naming,
    mut out: W,
) -> io::Result<()> {
    for line in bed.lines() {
        let line = line?;
        let Some(region) = Region::parse(&line) else {
            continue;
        };
        let sequence = genome.sequence(&region)?;
        writeln!(out, ">{}", region.header(naming))?;
        // bedtools writes the sequence on one line, however long it is.
        out.write_all(&sequence)?;
        writeln!(out)?;
    }
    out.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Builds a tiny genome and its index, the way `samtools faidx` would.
    fn tiny_genome(dir: &Path) -> std::path::PathBuf {
        let fasta = dir.join("g.fa");
        let mut f = File::create(&fasta).unwrap();
        // Two records, fixed line length so the index is easy to state.
        writeln!(f, ">chr1").unwrap();
        writeln!(f, "ACGTACGTAC").unwrap();
        writeln!(f, "GGGGCCCCAA").unwrap();
        writeln!(f, ">chr2").unwrap();
        writeln!(f, "acgtNNNNRY").unwrap();
        drop(f);

        // name, length, offset of the first base, bases per line, bytes per line
        let mut fai = File::create(dir.join("g.fa.fai")).unwrap();
        writeln!(fai, "chr1\t20\t6\t10\t11").unwrap();
        writeln!(fai, "chr2\t10\t34\t10\t11").unwrap();
        fasta
    }

    #[test]
    fn a_bed_line_becomes_a_region() {
        let r = Region::parse("chr1\t100\t200\tdomA\t255\t+").unwrap();
        assert_eq!(
            r,
            Region {
                name: "domA".into(),
                chrom: "chr1".into(),
                start: 100,
                end: 200
            }
        );
    }

    #[test]
    fn a_bed_line_without_a_name_is_skipped() {
        assert!(Region::parse("chr1\t100\t200").is_none());
    }

    /// The two conventions, side by side. The bare one is the only one the
    /// pipeline can join on.
    #[test]
    fn the_two_naming_conventions() {
        let r = Region::parse("chr1\t100\t200\tdomA\t255\t+").unwrap();
        assert_eq!(r.header(Naming::Bare), "domA");
        assert_eq!(r.header(Naming::WithCoordinates), "domA::chr1:100-200");
        assert_eq!(Naming::default(), Naming::Bare);
    }

    #[test]
    fn extracts_the_right_bases() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = tiny_genome(dir.path());
        let mut genome = Genome::open(&fasta).unwrap();

        let region = Region::parse("chr1\t0\t4\td\t0\t+").unwrap();
        assert_eq!(genome.sequence(&region).unwrap(), b"ACGT");

        // Spanning the line break in the FASTA.
        let region = Region::parse("chr1\t8\t12\td\t0\t+").unwrap();
        assert_eq!(genome.sequence(&region).unwrap(), b"ACGG");
    }

    /// The whole point: GC measured straight from the genome, with the same
    /// reading `gc_from_fasta.pl` uses.
    #[test]
    fn measures_gc_for_each_domain() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = tiny_genome(dir.path());
        let mut genome = Genome::open(&fasta).unwrap();

        let bed = "chr1\t10\t20\thigh\t255\t+\nchr2\t0\t10\tmixed\t255\t-\n";
        let rows = build(bed.as_bytes(), &mut genome).unwrap();

        assert_eq!(rows.len(), 2);
        // GGGGCCCCAA: eight of ten are not A, T or N.
        assert_eq!(rows[0].name, "high");
        assert_eq!(rows[0].counts.gc, 8);
        assert_eq!(rows[0].counts.length, 10);

        // acgtNNNNRY: the four lowercase and the two ambiguity codes count,
        // the four uppercase N do not.
        assert_eq!(rows[1].name, "mixed");
        assert_eq!(rows[1].counts.gc, 6);
        assert_eq!(rows[1].counts.soft_masked, 4);
        assert_eq!(rows[1].counts.ambiguous, 2);
    }

    /// A minus-strand feature is not reverse-complemented, matching `getfasta`
    /// without `-s`. GC is unchanged either way, which is the point.
    #[test]
    fn minus_strand_features_are_not_reverse_complemented() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = tiny_genome(dir.path());
        let mut genome = Genome::open(&fasta).unwrap();

        let region = Region::parse("chr1\t0\t4\td\t0\t-").unwrap();
        assert_eq!(genome.sequence(&region).unwrap(), b"ACGT", "not TGCA");
    }

    #[test]
    fn the_intermediate_fasta_uses_the_bare_name_by_default() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = tiny_genome(dir.path());
        let mut genome = Genome::open(&fasta).unwrap();

        let mut out = Vec::new();
        write_fasta(
            "chr1\t0\t4\tdomA\t255\t+\n".as_bytes(),
            &mut genome,
            Naming::Bare,
            &mut out,
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), ">domA\nACGT\n");
    }

    #[test]
    fn a_missing_index_is_reported_clearly() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("no_index.fa");
        File::create(&fasta).unwrap();
        let Err(err) = Genome::open(&fasta) else {
            panic!("opening a FASTA without an index should fail");
        };
        assert!(
            err.to_string().contains("FASTA index"),
            "unhelpful message: {err}"
        );
    }
}
