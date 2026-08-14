//! Alignment, in process, through [bwa-mem4](https://crates.io/crates/bwa-mem4).
//!
//! This is the one stage the paper's pipeline shells out for and this one does
//! not. `bwa-mem4` exposes its command implementations as a library, so the
//! alignment runs on the same code path its binary takes and produces the same
//! bytes, without a subprocess to spawn or a stdout to parse.
//!
//! **This is not bowtie2, and the copy numbers will differ.** `DIVERGENCES.md`
//! says so at length: BWA-MEM aligns locally with soft-clipping and emits
//! supplementary records, and on a locus with more than three hundred paralogous
//! copies the multimappers land differently. Use [`crate::coverage`] on a
//! bowtie2 BAM for the paper's numbers; use this when a self-consistent
//! measurement across a cohort is what is wanted.
//!
//! # Why the SAM detour
//!
//! `bwa-mem4` can write BAM itself, but only through its `multi-format`
//! feature, which pulls htslib. Asking it for SAM and transcoding with noodles
//! keeps the C out of the picture, and costs one pass over a file that is being
//! written and read back on the same machine anyway.

use std::io::BufRead;
use std::path::{Path, PathBuf};

// The BAM writer takes a SAM record through this trait rather than through
// `write_record`, which wants a `bam::Record`: the record is encoded straight
// into the BAM block, so its fields are never materialised in between.
use noodles_sam::alignment::io::Write as _;

use bwa_mem4::cmd_index::IndexArgs;
use bwa_mem4::cmd_mem::MemArgs;
use noodles_bam as bam;
use noodles_sam as sam;

/// What went wrong, in the caller's terms rather than the aligner's.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("the index {0} is missing; build it with `plethora index <genome.fa>`")]
    NoIndex(PathBuf),
    #[error("{0} does not exist")]
    NoReads(PathBuf),
    #[error("the aligner failed: {0}")]
    Aligner(#[from] anyhow::Error),
    #[error("converting the alignment to BAM: {0}")]
    Convert(#[from] std::io::Error),
}

/// How to run one alignment.
#[derive(Debug, Clone)]
pub struct Options {
    /// The prefix given to `plethora index`, which is the FASTA path itself
    /// unless `-p` said otherwise.
    pub index_prefix: PathBuf,
    pub read1: PathBuf,
    /// `None` aligns single-end.
    pub read2: Option<PathBuf>,
    pub threads: usize,
    /// `-K`, the batch size in input bases.
    ///
    /// Worth setting for anything that has to be reproducible. Paired-end
    /// insert-size statistics are estimated once per batch, so moving a batch
    /// boundary moves MAPQ and pairing decisions for the reads near it; pinning
    /// `-K` makes a run repeatable across thread counts. Left `None` the
    /// aligner uses `10M * threads`, which is thread-dependent by construction.
    pub batch_bases: Option<i64>,
    /// `@RG`, which a cohort usually wants and a single sample rarely does.
    pub read_group: Option<String>,
}

impl Options {
    /// The common case: one paired sample against one index.
    #[must_use]
    pub fn paired(index_prefix: &Path, read1: &Path, read2: &Path, threads: usize) -> Self {
        Self {
            index_prefix: index_prefix.to_path_buf(),
            read1: read1.to_path_buf(),
            read2: Some(read2.to_path_buf()),
            threads,
            batch_bases: None,
            read_group: None,
        }
    }
}

/// The five files `bwa-mem4 index` writes beside the prefix.
///
/// Checked before the aligner runs so a missing index is one clear sentence
/// rather than whatever the FM-index loader says about a file it cannot open.
const INDEX_SUFFIXES: &[&str] = &[".pac", ".ann", ".amb", ".bwt.2bit.64", ".0123"];

/// Builds the index, which is `bwa-mem4 index`.
///
/// # Errors
/// Returns an error if the FASTA cannot be read or the index cannot be written.
pub fn build_index(fasta: &Path, prefix: Option<&Path>) -> Result<(), Error> {
    if !fasta.exists() {
        return Err(Error::NoReads(fasta.to_path_buf()));
    }
    bwa_mem4::cmd_index::run(IndexArgs {
        fasta: fasta.to_path_buf(),
        prefix: prefix.map(Path::to_path_buf),
        // The minimiser index the long-read path wants. Not built: plethora
        // aligns short reads, and it roughly doubles index time.
        mmi: Vec::new(),
    })?;
    Ok(())
}

/// Aligns one sample and writes a BAM.
///
/// The SAM the aligner produces goes to a temporary file beside the output and
/// is transcoded; nothing is left behind on either path.
///
/// # Errors
/// Returns an error if the index or the reads are missing, if the aligner
/// fails, or if the BAM cannot be written.
pub fn align(options: &Options, output: &Path) -> Result<(), Error> {
    for suffix in INDEX_SUFFIXES {
        let mut path = options.index_prefix.clone().into_os_string();
        path.push(suffix);
        let path = PathBuf::from(path);
        if !path.exists() {
            return Err(Error::NoIndex(path));
        }
    }
    for reads in [Some(&options.read1), options.read2.as_ref()]
        .into_iter()
        .flatten()
    {
        if !reads.exists() {
            return Err(Error::NoReads(reads.clone()));
        }
    }
    if let Some(parent) = output.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)?;
    }

    // Beside the output rather than in the system temporary directory: a
    // whole-genome SAM is tens of gigabytes, and /tmp is usually the smallest
    // filesystem on a cluster node.
    let sam_path = output.with_extension("sam.tmp");
    let result = run_aligner(options, &sam_path)
        .and_then(|()| sam_to_bam(&sam_path, output).map_err(Error::Convert));
    let _ = std::fs::remove_file(&sam_path);
    result
}

fn run_aligner(options: &Options, sam_path: &Path) -> Result<(), Error> {
    let args = MemArgs {
        index_prefix: options.index_prefix.clone(),
        reads: options.read1.clone(),
        reads2: options.read2.clone(),
        threads: i32::try_from(options.threads.max(1)).unwrap_or(i32::MAX),
        k_batch: options.batch_bases,
        read_group: options.read_group.clone(),
        output: Some(sam_path.to_path_buf()),
        ..Default::default()
    };
    // This lands in the SAM `@PG CL:` record and nowhere else: it is provenance
    // for the file, not a second source of options.
    let argv = vec![
        "plethora".to_string(),
        "align".to_string(),
        options.index_prefix.display().to_string(),
    ];
    bwa_mem4::cmd_mem::run(args, &argv)?;
    Ok(())
}

/// Transcodes SAM to BAM, record for record.
fn sam_to_bam(sam_path: &Path, output: &Path) -> std::io::Result<()> {
    let mut reader = sam::io::Reader::new(crate::io::open(sam_path)?);
    let header = read_header(&mut reader)?;

    let mut writer = bam::io::Writer::new(std::fs::File::create(output)?);
    writer.write_alignment_header(&header)?;

    let mut record = sam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        writer.write_alignment_record(&header, &record)?;
    }
    Ok(())
}

/// The header, with the one thing a SAM from an aligner can lack.
///
/// `noodles` will not write a BAM without reference sequences, and it takes
/// them from the header rather than from the records. bwa-mem4 always writes
/// `@SQ` lines, so this only ever fails on a truncated file, but failing with
/// the reason beats failing with an index out of range further down.
fn read_header<R: BufRead>(reader: &mut sam::io::Reader<R>) -> std::io::Result<sam::Header> {
    let header = reader.read_header()?;
    if header.reference_sequences().is_empty() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "the alignment has no @SQ lines, so it names no reference sequences",
        ));
    }
    Ok(header)
}
