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

use flate2::read::MultiGzDecoder;

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
    let file = match header(path)? {
        // BGZF is a chain of independently deflated blocks, so it decodes on
        // every core. Measured on a 58 MB intermediate: 35 ms sequentially,
        // 6.5 ms in parallel.
        Header::Bgzf => {
            let decoder = rapidgzip_core::Decoder::builder()
                .build()
                .map_err(io::Error::other)?;
            return Ok(Box::new(BufReader::new(
                decoder.open(path).map_err(io::Error::other)?,
            )));
        }
        // A single deflate stream has no block boundaries to split on, so the
        // parallel decoder has to guess where they are and comes out slower:
        // 46 ms against flate2's 36 ms on the same 58 MB. Foreign `.gz` inputs
        // are usually this, so they keep the sequential path.
        Header::Gzip => {
            return Ok(Box::new(BufReader::new(MultiGzDecoder::new(File::open(
                path,
            )?))));
        }
        Header::Plain => File::open(path)?,
    };
    Ok(Box::new(BufReader::new(file)))
}

/// What the first bytes say the file is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Header {
    /// gzip carrying the BGZF `BC` extra subfield, so block-framed.
    Bgzf,
    /// gzip, one deflate stream as far as the header says.
    Gzip,
    Plain,
}

/// Sniffs the header, since the name is not evidence and BGZF is spelled `.gz`.
///
/// BGZF is gzip with `FEXTRA` set and a `BC` subfield giving the compressed
/// block size, which is what makes the blocks findable without decoding them.
/// The layout is fixed by the SAM specification: eighteen bytes of header, with
/// the subfield identifier at offsets twelve and thirteen.
fn header(path: &Path) -> io::Result<Header> {
    // FLG bit 2 is FEXTRA. Without it there is no subfield to look for.
    const FEXTRA: u8 = 0b0000_0100;

    let mut probe = File::open(path)?;
    let mut head = [0_u8; 16];
    let read = read_up_to(&mut probe, &mut head)?;

    if read < 2 || head[..2] != GZIP_MAGIC {
        return Ok(Header::Plain);
    }
    if read >= 14 && head[3] & FEXTRA != 0 && head[12] == b'B' && head[13] == b'C' {
        return Ok(Header::Bgzf);
    }
    Ok(Header::Gzip)
}

/// Reads as much as the buffer holds, returning how much arrived.
fn read_up_to(file: &mut File, buffer: &mut [u8]) -> io::Result<usize> {
    let mut filled = 0;
    while filled < buffer.len() {
        match file.read(&mut buffer[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(filled)
}

/// Creates a file for writing, compressing when the name ends in `.gz`.
///
/// # Errors
/// Returns an error if the file cannot be created.
pub fn create(path: &Path) -> io::Result<Box<dyn Write>> {
    let file = File::create(path)?;
    if !is_gzip_name(path) {
        return Ok(Box::new(BufWriter::new(file)));
    }
    // BGZF rather than a single deflate stream, which every gzip reader still
    // takes. It compresses on every core and, because its blocks are framed,
    // decodes on every core too. Measured on a 58 MB intermediate: 210 ms to
    // write as gzip against 51 ms as BGZF, 13.8 MB against 11.6 MB, and 46 ms
    // to read back in parallel against 6.5 ms.
    let workers = std::thread::available_parallelism().unwrap_or(std::num::NonZero::<usize>::MIN);
    let level = noodles_bgzf::io::writer::CompressionLevel::new(COMPRESSION_LEVEL as u8)
        .unwrap_or(noodles_bgzf::io::writer::CompressionLevel::FAST);
    Ok(Box::new(
        noodles_bgzf::io::multithreaded_writer::Builder::default()
            .set_worker_count(workers)
            .set_compression_level(level)
            .build_from_writer(file),
    ))
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

        // Really compressed, not merely named that way, and BGZF rather than a
        // single deflate stream. On a file this small the block framing costs
        // more than it saves, which is why the bound is loose; on the files
        // that matter it is 11.6 MB where gzip is 13.8.
        let raw = std::fs::read(&path).unwrap();
        assert_eq!(&raw[..2], &GZIP_MAGIC, "no gzip header");
        assert_eq!(
            header(&path).unwrap(),
            Header::Bgzf,
            "the writer should frame its blocks"
        );
        assert!(
            raw.len() < 12_990,
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
