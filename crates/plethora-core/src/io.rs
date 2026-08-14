//! Reading and writing the pipeline's files, compressed or not.
//!
//! Upstream writes every intermediate uncompressed and leaves compression to
//! whoever archives the run: `zip.sh` gzips the sorted BED afterwards, and
//! cohorts keep `_read_depth.bed.gz` and `_gc_correct.txt.gz`. A 623,699-domain
//! `_gc_correct.txt` is 40 MB and about 10 MB gzipped, and there is one per
//! sample.
//!
//! So compression is decided by the file name here, in one place: a path ending
//! in `.gz` is written through gzip, and any input starting with the gzip magic
//! is read through it whatever it is called. Nothing else in the pipeline needs
//! to know.

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;

/// The gzip magic number, which is what decides how an input is read.
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Compression level for the outputs this writes.
///
/// Level 1 rather than the default 6: these files are written once and read
/// once, and the difference in size is a few percent against several times the
/// CPU. `trim-galore` makes the same choice for the same reason.
pub const COMPRESSION_LEVEL: u32 = 1;

/// Opens a file for reading, decompressing if it is gzipped.
///
/// Detection is on the content, not the name, so a `.bed` that turns out to be
/// gzipped still reads. Multi-member archives are handled, which matters
/// because a file concatenated from several gzip streams is otherwise
/// truncated silently at the first member's end.
///
/// # Errors
/// Returns an error if the file cannot be opened or its first bytes cannot be
/// read.
pub fn open(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let mut probe = File::open(path)?;
    let mut magic = [0_u8; 2];
    let gzipped = match probe.read_exact(&mut magic) {
        Ok(()) => magic == GZIP_MAGIC,
        // A file shorter than two bytes holds no gzip member.
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => false,
        Err(e) => return Err(e),
    };

    let file = File::open(path)?;
    Ok(if gzipped {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    })
}

/// Creates a file for writing, compressing when the name ends in `.gz`.
///
/// # Errors
/// Returns an error if the file cannot be created.
pub fn create(path: &Path) -> io::Result<Box<dyn Write>> {
    let file = File::create(path)?;
    Ok(if is_gzip_name(path) {
        Box::new(GzEncoder::new(
            BufWriter::new(file),
            Compression::new(COMPRESSION_LEVEL),
        ))
    } else {
        Box::new(BufWriter::new(file))
    })
}

/// True when the name asks for compression.
#[must_use]
pub fn is_gzip_name(path: &Path) -> bool {
    path.extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"))
}

/// Adds `.gz` to a path, or leaves it alone if it already ends that way.
#[must_use]
pub fn with_gzip(path: &Path) -> PathBuf {
    if is_gzip_name(path) {
        return path.to_path_buf();
    }
    let mut name = path.as_os_str().to_os_string();
    name.push(".gz");
    PathBuf::from(name)
}

/// Removes a trailing `.gz`, or leaves the path alone.
#[must_use]
pub fn without_gzip(path: &Path) -> PathBuf {
    if !is_gzip_name(path) {
        return path.to_path_buf();
    }
    path.with_extension("")
}

/// Whether an output should be written compressed.
///
/// Carried through the pipeline as one value rather than a `.gz` suffix decided
/// per stage, so a run's outputs are all one way or all the other.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Compress {
    /// Write plain files, as upstream does.
    #[default]
    No,
    /// Append `.gz` to every output and write through gzip.
    Yes,
}

impl Compress {
    /// The name an output should take under this setting.
    #[must_use]
    pub fn apply(self, path: &Path) -> PathBuf {
        match self {
            Self::No => without_gzip(path),
            Self::Yes => with_gzip(path),
        }
    }
}

impl std::str::FromStr for Compress {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "yes" | "gz" | "gzip" | "true" => Ok(Self::Yes),
            "no" | "none" | "plain" | "false" => Ok(Self::No),
            other => Err(format!("expected yes or no, got {other}")),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_plain_file_round_trips() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("plain.bed");
        let mut w = create(&path).unwrap();
        writeln!(w, "chr1\t1\t2").unwrap();
        drop(w);

        let mut text = String::new();
        open(&path).unwrap().read_to_string(&mut text).unwrap();
        assert_eq!(text, "chr1\t1\t2\n");
    }

    #[test]
    fn a_gz_name_is_written_compressed_and_read_back() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("compressed.bed.gz");
        let mut w = create(&path).unwrap();
        for i in 0..1000 {
            writeln!(w, "chr1\t{i}\t{}", i + 100).unwrap();
        }
        drop(w);

        // Really compressed, not merely named that way.
        let raw = std::fs::read(&path).unwrap();
        assert_eq!(&raw[..2], &GZIP_MAGIC, "no gzip header");
        assert!(
            raw.len() < 5000,
            "not actually compressed: {} bytes",
            raw.len()
        );

        let mut text = String::new();
        open(&path).unwrap().read_to_string(&mut text).unwrap();
        assert_eq!(text.lines().count(), 1000);
        assert!(text.starts_with("chr1\t0\t100\n"));
    }

    /// Detection is on the content, so a gzipped file with a misleading name
    /// still reads. Cohort archives are full of these.
    #[test]
    fn a_gzipped_file_reads_whatever_it_is_called() {
        let dir = tempfile::tempdir().unwrap();
        let real = dir.path().join("x.gz");
        let mut w = create(&real).unwrap();
        writeln!(w, "hello").unwrap();
        drop(w);

        let misleading = dir.path().join("x.txt");
        std::fs::rename(&real, &misleading).unwrap();

        let mut text = String::new();
        open(&misleading)
            .unwrap()
            .read_to_string(&mut text)
            .unwrap();
        assert_eq!(text, "hello\n");
    }

    /// Several gzip members concatenated, which is what `cat a.gz b.gz` gives
    /// and what a plain decoder truncates at the first.
    #[test]
    fn a_multi_member_archive_reads_in_full() {
        let dir = tempfile::tempdir().unwrap();
        let a = dir.path().join("a.gz");
        let b = dir.path().join("b.gz");
        for (path, line) in [(&a, "first"), (&b, "second")] {
            let mut w = create(path).unwrap();
            writeln!(w, "{line}").unwrap();
        }

        let joined = dir.path().join("joined.gz");
        let mut bytes = std::fs::read(&a).unwrap();
        bytes.extend(std::fs::read(&b).unwrap());
        std::fs::write(&joined, bytes).unwrap();

        let mut text = String::new();
        open(&joined).unwrap().read_to_string(&mut text).unwrap();
        assert_eq!(text, "first\nsecond\n");
    }

    #[test]
    fn an_empty_file_is_not_mistaken_for_gzip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.bed");
        std::fs::write(&path, b"").unwrap();
        let mut text = String::new();
        open(&path).unwrap().read_to_string(&mut text).unwrap();
        assert!(text.is_empty());
    }

    #[test]
    fn names_gain_and_lose_the_suffix_idempotently() {
        let plain = Path::new("results/S1_read_depth.bed");
        let gz = Path::new("results/S1_read_depth.bed.gz");

        assert_eq!(with_gzip(plain), gz);
        assert_eq!(with_gzip(gz), gz, "already compressed, left alone");
        assert_eq!(without_gzip(gz), plain);
        assert_eq!(without_gzip(plain), plain);
    }

    #[test]
    fn the_setting_decides_the_name() {
        let plain = Path::new("results/S1_gc_correct.txt");
        assert_eq!(Compress::No.apply(plain), plain);
        assert_eq!(
            Compress::Yes.apply(plain),
            Path::new("results/S1_gc_correct.txt.gz")
        );
        // And it is idempotent from either starting point.
        let gz = Path::new("results/S1_gc_correct.txt.gz");
        assert_eq!(Compress::Yes.apply(gz), gz);
        assert_eq!(Compress::No.apply(gz), plain);
    }

    #[test]
    fn the_setting_parses_the_spellings_a_caller_might_use() {
        for yes in ["yes", "gz", "gzip", "true"] {
            assert_eq!(yes.parse::<Compress>().unwrap(), Compress::Yes);
        }
        for no in ["no", "none", "plain", "false"] {
            assert_eq!(no.parse::<Compress>().unwrap(), Compress::No);
        }
        assert!("maybe".parse::<Compress>().is_err());
        assert_eq!(Compress::default(), Compress::No);
    }
}
