//! `samtools sort -n`, including the part that reaches the output.
//!
//! ```text
//! samtools sort -n -@ 12 -m 2G -o ${output}_sorted.bam $bam
//! ```
//!
//! The obvious purpose is to put mates next to each other so `bamtobed -bedpe`
//! can pair them. Both bowtie2 and BWA-MEM already emit mates adjacently, so
//! for that purpose the sort is nearly a no-op.
//!
//! The purpose that actually matters is the order itself. `merge_pairs.pl`
//! stops accumulating fragment lengths after 50 million records, and a
//! whole-genome sample has more pairs than that, so the mean and standard
//! deviation it derives, and therefore the length of every extended single-end
//! read in the output, depend on which 50 million came first. Getting the order
//! wrong produces a file that looks right and carries different numbers.
//!
//! The comparison is [`plethora_compat::strnum::cmp_by_qname`]; this module is
//! the sort around it. samtools merges stably, so records tying on both name
//! and pair rank keep their input order, and that is reproduced.

use std::cmp::Ordering;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;

use plethora_compat::strnum::cmp_by_qname;

use super::bamtobed::Aln;

/// How many records a sorting run holds before it spills to disk.
///
/// At roughly 80 bytes a record this is about 400 MB, in the range of the
/// `-m 2G` per thread upstream gives samtools. Configurable because the tests
/// need to force spilling on inputs small enough to eyeball.
pub const DEFAULT_RUN_RECORDS: usize = 5_000_000;

/// Sorts records by name, spilling runs to disk when they exceed `run_records`.
///
/// Returns the records in `samtools sort -n` order. Stability is part of the
/// contract, not an accident: samtools' merge preserves input order on a full
/// tie, and `merge_pairs.pl` reads the result positionally.
///
/// # Errors
/// Returns an error if a temporary run file cannot be written or read back.
///
/// # Panics
/// Panics if `run_records` is zero, which would never make progress.
pub fn sort_by_name<I>(records: I, run_records: usize, tmp_dir: &Path) -> io::Result<Vec<Aln>>
where
    I: IntoIterator<Item = Aln>,
{
    assert!(run_records > 0, "a sorting run must hold at least one record");

    let mut runs: Vec<std::path::PathBuf> = Vec::new();
    let mut buffer: Vec<Aln> = Vec::new();
    let mut spilled = 0_usize;

    for record in records {
        buffer.push(record);
        if buffer.len() >= run_records {
            sort_run(&mut buffer);
            let path = tmp_dir.join(format!("namesort-{spilled}.run"));
            write_run(&buffer, &path)?;
            runs.push(path);
            spilled += 1;
            buffer.clear();
        }
    }

    if runs.is_empty() {
        // Everything fit, so no merge is needed and no temporary file was made.
        sort_run(&mut buffer);
        return Ok(buffer);
    }

    if !buffer.is_empty() {
        sort_run(&mut buffer);
        let path = tmp_dir.join(format!("namesort-{spilled}.run"));
        write_run(&buffer, &path)?;
        runs.push(path);
    }

    let merged = merge_runs(&runs)?;
    for path in &runs {
        // A failure to clean up is not a failure to sort.
        let _ = std::fs::remove_file(path);
    }
    Ok(merged)
}

/// Sorts one run in place, stably.
fn sort_run(run: &mut [Aln]) {
    run.sort_by(|a, b| cmp_by_qname(&a.name, a.flags, &b.name, b.flags));
}

/// Merges sorted runs, preserving input order across runs on a full tie.
///
/// A k-way merge that always takes the earliest run among equal keys, which is
/// what keeps the whole sort stable: run `i` holds records that came before run
/// `j > i` in the input.
fn merge_runs(paths: &[std::path::PathBuf]) -> io::Result<Vec<Aln>> {
    let mut readers: Vec<RunReader> = paths.iter().map(RunReader::open).collect::<io::Result<_>>()?;
    let mut heads: Vec<Option<Aln>> = readers.iter_mut().map(RunReader::next).collect::<io::Result<_>>()?;

    let mut out = Vec::new();
    loop {
        let mut best: Option<usize> = None;
        for (i, head) in heads.iter().enumerate() {
            let Some(candidate) = head else { continue };
            match best {
                None => best = Some(i),
                Some(b) => {
                    let current = heads[b].as_ref().expect("best index holds a record");
                    // Strictly less, so an equal key leaves the earlier run in
                    // front and the merge stays stable.
                    if cmp_by_qname(&candidate.name, candidate.flags, &current.name, current.flags)
                        == Ordering::Less
                    {
                        best = Some(i);
                    }
                }
            }
        }

        let Some(i) = best else { break };
        out.push(heads[i].take().expect("chosen run holds a record"));
        heads[i] = readers[i].next()?;
    }

    Ok(out)
}

/// Writes a sorted run in a compact length-prefixed form.
fn write_run(run: &[Aln], path: &Path) -> io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    for a in run {
        let name_len = u16::try_from(a.name.len())
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "QNAME longer than 65535"))?;
        w.write_all(&name_len.to_le_bytes())?;
        w.write_all(&a.name)?;
        w.write_all(&a.flags.to_le_bytes())?;

        match &a.chrom {
            Some(c) => {
                let len = u16::try_from(c.len()).map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidData, "reference name longer than 65535")
                })?;
                w.write_all(&len.to_le_bytes())?;
                w.write_all(c.as_bytes())?;
            }
            None => w.write_all(&u16::MAX.to_le_bytes())?,
        }

        w.write_all(&a.start.to_le_bytes())?;
        w.write_all(&a.end.to_le_bytes())?;
        w.write_all(&[a.mapq])?;
    }
    w.flush()
}

/// Reads back a run written by [`write_run`].
struct RunReader {
    inner: BufReader<File>,
}

impl RunReader {
    fn open(path: &std::path::PathBuf) -> io::Result<Self> {
        Ok(Self {
            inner: BufReader::new(File::open(path)?),
        })
    }

    fn next(&mut self) -> io::Result<Option<Aln>> {
        let mut len = [0_u8; 2];
        match self.inner.read_exact(&mut len) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
            Err(e) => return Err(e),
        }

        let mut name = vec![0_u8; usize::from(u16::from_le_bytes(len))];
        self.inner.read_exact(&mut name)?;

        let mut flags = [0_u8; 2];
        self.inner.read_exact(&mut flags)?;

        let mut chrom_len = [0_u8; 2];
        self.inner.read_exact(&mut chrom_len)?;
        let chrom_len = u16::from_le_bytes(chrom_len);
        let chrom = if chrom_len == u16::MAX {
            None
        } else {
            let mut buf = vec![0_u8; usize::from(chrom_len)];
            self.inner.read_exact(&mut buf)?;
            Some(String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?)
        };

        let mut start = [0_u8; 8];
        self.inner.read_exact(&mut start)?;
        let mut end = [0_u8; 8];
        self.inner.read_exact(&mut end)?;
        let mut mapq = [0_u8; 1];
        self.inner.read_exact(&mut mapq)?;

        Ok(Some(Aln {
            name,
            flags: u16::from_le_bytes(flags),
            chrom,
            start: i64::from_le_bytes(start),
            end: i64::from_le_bytes(end),
            mapq: mapq[0],
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::bamtobed::{FIRST_MATE, PAIRED, SECOND_MATE};

    fn aln(name: &str, flags: u16, start: i64) -> Aln {
        Aln {
            name: name.as_bytes().to_vec(),
            flags,
            chrom: Some("chr1".to_string()),
            start,
            end: start + 50,
            mapq: 60,
        }
    }

    fn names(sorted: &[Aln]) -> Vec<String> {
        sorted
            .iter()
            .map(|a| String::from_utf8_lossy(&a.name).into_owned())
            .collect()
    }

    #[test]
    fn sorts_naturally_not_lexically() {
        let dir = tempfile::tempdir().unwrap();
        let input = vec![
            aln("r10", PAIRED | FIRST_MATE, 0),
            aln("r2", PAIRED | FIRST_MATE, 0),
            aln("r1", PAIRED | FIRST_MATE, 0),
        ];
        let out = sort_by_name(input, 100, dir.path()).unwrap();
        assert_eq!(names(&out), ["r1", "r2", "r10"]);
    }

    #[test]
    fn read1_precedes_read2_under_one_name() {
        let dir = tempfile::tempdir().unwrap();
        let input = vec![
            aln("r1", PAIRED | SECOND_MATE, 300),
            aln("r1", PAIRED | FIRST_MATE, 100),
        ];
        let out = sort_by_name(input, 100, dir.path()).unwrap();
        assert_eq!(out[0].flags & FIRST_MATE, FIRST_MATE);
        assert_eq!(out[1].flags & SECOND_MATE, SECOND_MATE);
    }

    /// Full ties keep input order, which is what samtools' merge does and what
    /// `merge_pairs.pl` relies on.
    #[test]
    fn equal_keys_keep_input_order() {
        let dir = tempfile::tempdir().unwrap();
        let input = vec![
            aln("same", PAIRED | FIRST_MATE, 1),
            aln("same", PAIRED | FIRST_MATE, 2),
            aln("same", PAIRED | FIRST_MATE, 3),
        ];
        let out = sort_by_name(input, 100, dir.path()).unwrap();
        assert_eq!(out.iter().map(|a| a.start).collect::<Vec<_>>(), [1, 2, 3]);
    }

    /// Spilling must not change the answer, and must stay stable across runs.
    #[test]
    fn spilling_to_disk_gives_the_same_order() {
        let dir = tempfile::tempdir().unwrap();
        let mut input = Vec::new();
        for i in 0..50 {
            input.push(aln(&format!("q{}", 50 - i), PAIRED | FIRST_MATE, i));
            input.push(aln(&format!("q{}", 50 - i), PAIRED | SECOND_MATE, i));
        }

        let in_memory = sort_by_name(input.clone(), 1000, dir.path()).unwrap();
        let spilled = sort_by_name(input.clone(), 7, dir.path()).unwrap();
        assert_eq!(in_memory, spilled, "spilling changed the order");

        // And a run size of one, so every record is its own run.
        let each = sort_by_name(input, 1, dir.path()).unwrap();
        assert_eq!(in_memory, each);
    }

    /// Stability must survive the merge, not only the in-run sort.
    #[test]
    fn stability_survives_spilling() {
        let dir = tempfile::tempdir().unwrap();
        let input: Vec<Aln> = (0..20).map(|i| aln("same", PAIRED | FIRST_MATE, i)).collect();
        let spilled = sort_by_name(input, 3, dir.path()).unwrap();
        assert_eq!(
            spilled.iter().map(|a| a.start).collect::<Vec<_>>(),
            (0..20).collect::<Vec<_>>()
        );
    }

    #[test]
    fn unmapped_records_round_trip_through_a_run() {
        let dir = tempfile::tempdir().unwrap();
        let mut unmapped = aln("u", PAIRED | FIRST_MATE, 0);
        unmapped.chrom = None;
        let input = vec![unmapped.clone(), aln("a", PAIRED | FIRST_MATE, 5)];
        let out = sort_by_name(input, 1, dir.path()).unwrap();
        assert_eq!(names(&out), ["a", "u"]);
        assert_eq!(out[1], unmapped);
    }

    #[test]
    fn an_empty_input_sorts_to_nothing() {
        let dir = tempfile::tempdir().unwrap();
        assert!(sort_by_name(Vec::new(), 10, dir.path()).unwrap().is_empty());
    }
}
