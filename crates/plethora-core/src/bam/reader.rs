//! Reading a BAM into the narrow view [`Aln`] that the interval stages need.
//!
//! The whole pipeline touches six fields of an alignment. Projecting onto them
//! at the boundary keeps the bedtools and merge-pairs rules testable without
//! constructing BAM files, and keeps `noodles` out of everything downstream.

use std::io;
use std::path::Path;

use noodles_sam::alignment::Record as _;

use super::bamtobed::Aln;

/// Reads a BAM file, yielding one [`Aln`] per record in file order.
///
/// # Errors
/// Returns an error if the file cannot be opened, the header cannot be parsed,
/// or a record is malformed.
pub fn read_bam<P: AsRef<Path>>(path: P) -> io::Result<Vec<Aln>> {
    read_bam_streaming(path)?.collect()
}

/// [`read_bam`] without holding the file in memory.
///
/// A 30x whole genome is upwards of a billion records, and at 181 bytes each a
/// `Vec` of them is 137 GB. Every stage downstream reads them once, in order,
/// so nothing needed the `Vec` in the first place.
///
/// # Errors
/// Returns an error if the file cannot be opened or the header cannot be
/// parsed. A malformed record surfaces from the iterator.
pub fn read_bam_streaming<P: AsRef<Path>>(
    path: P,
) -> io::Result<impl Iterator<Item = io::Result<Aln>>> {
    let mut reader = noodles_bam::io::reader::Builder.build_from_path(path)?;
    let header = reader.read_header()?;

    // Reference names, indexed by reference id, resolved once.
    let names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
        .collect();

    // `records()` borrows the reader, so the reader moves into the iterator
    // and yields projected records as it goes.
    Ok(std::iter::from_fn(move || {
        let mut record = noodles_bam::Record::default();
        match reader.read_record(&mut record) {
            Ok(0) => None,
            Ok(_) => Some(project(&record, &names)),
            Err(e) => Some(Err(e)),
        }
    }))
}

/// Projects one record onto the fields the pipeline reads.
///
/// The end position comes from `alignment_end`, which sums only the
/// reference-consuming CIGAR operations (M, D, N, =, X). That is the same set
/// `BamTools`' `GetEndPosition(false)` uses, so a soft-clipped read spans exactly
/// as much here as it does under `bedtools bamtobed`.
fn project(record: &noodles_bam::Record, names: &[String]) -> io::Result<Aln> {
    let flags = u16::from(record.flags());

    let name = record.name().map(|n| n.to_vec()).unwrap_or_default();
    let mapq = record.mapping_quality().map_or(0, u8::from);

    let mapped = flags & super::bamtobed::UNMAPPED == 0;
    let (chrom, start, end) = if mapped {
        let id = record.reference_sequence_id().transpose()?.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "mapped record has no reference")
        })?;
        let chrom = names.get(id).cloned().ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "reference id out of range")
        })?;

        let start = record.alignment_start().transpose()?.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "mapped record has no position")
        })?;
        let end = record.alignment_end().transpose()?.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "mapped record has no end")
        })?;

        // BED is zero-based half-open; noodles positions are one-based inclusive.
        (
            Some(chrom),
            (usize::from(start) - 1) as i64,
            usize::from(end) as i64,
        )
    } else {
        (None, 0, 0)
    };

    Ok(Aln {
        name,
        flags,
        chrom,
        start,
        end,
        mapq,
    })
}
