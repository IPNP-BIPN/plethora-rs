//! `download_sample.pl`: fetch a sample's FASTQ files and check them.
//!
//! ```text
//! code/download_sample.pl HG00250 data/1000Genomes_samples.txt
//! ```
//!
//! The script's real content is not the fetching, it is what it does when a
//! fetch goes wrong: verify the MD5 from the sequence index, retry up to ten
//! times across the sample as a whole, and clean up wget's numbered leftovers
//! before trying again. That logic is what is ported here, and the fetching
//! itself sits behind [`Fetcher`] so it can be exercised without a network.
//!
//! One thing has to change. The index gives `ftp://` URLs, and the 1000 Genomes
//! FTP endpoints have not been reliable for years; both EBI and NCBI serve the
//! same paths over HTTPS. [`https_url`] rewrites them, and says so.

use std::io::Read;
use std::path::{Path, PathBuf};

use md5::{Digest, Md5};

use super::sample_index::Record;

/// How many failures across one sample before giving up.
///
/// Upstream counts across the whole sample rather than per file, and resets the
/// counter after a success, so a sample with one bad file and nine good ones
/// still completes.
pub const MAX_TRIES: usize = 10;

/// Something that can fetch a URL to a path.
///
/// Behind a trait because the interesting behaviour is the retrying and the
/// verification, and neither should need a network to test.
pub trait Fetcher {
    /// Fetches `url` into `destination`, replacing whatever is there.
    ///
    /// # Errors
    /// Returns an error if the transfer fails. A partial file may be left
    /// behind; the caller removes it before retrying, as upstream does.
    fn fetch(&self, url: &str, destination: &Path) -> anyhow::Result<()>;
}

/// What became of one file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Outcome {
    /// Already present with the right checksum, so not fetched.
    AlreadyValid,
    /// Fetched and verified.
    Fetched {
        /// How many attempts it took, one for a clean first try.
        attempts: usize,
    },
}

/// A file that could not be obtained.
#[derive(Debug, thiserror::Error)]
pub enum DownloadError {
    #[error("exceeded the maximum number of tries ({MAX_TRIES}) for sample {sample}")]
    TooManyTries { sample: String },

    #[error("invalid checksum for {file}: should be {expected}, but is {found}")]
    BadChecksum {
        file: String,
        expected: String,
        found: String,
    },

    #[error(transparent)]
    Io(#[from] std::io::Error),

    #[error(transparent)]
    Other(#[from] anyhow::Error),
}

/// Rewrites an `ftp://` URL from the sequence index to `https://`.
///
/// The index was written when both EBI and NCBI served these paths over FTP.
/// Both now serve the same paths over HTTPS, and the FTP endpoints time out
/// often enough that following the file as written is the main reason a
/// download fails. Anything already `http`, and anything else, is left alone.
#[must_use]
pub fn https_url(url: &str) -> String {
    url.strip_prefix("ftp://")
        .map_or_else(|| url.to_string(), |rest| format!("https://{rest}"))
}

/// The MD5 of a file, lowercase hex, as `md5sum` prints it.
///
/// # Errors
/// Returns an error if the file cannot be read.
pub fn md5_of(path: &Path) -> std::io::Result<String> {
    let mut file = std::fs::File::open(path)?;
    let mut hasher = Md5::new();
    let mut buffer = vec![0_u8; 1 << 20];
    loop {
        let read = file.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

/// Where a sample's files go: `fastq/<sample>/`.
#[must_use]
pub fn sample_dir(root: &Path, sample: &str) -> PathBuf {
    root.join("fastq").join(sample)
}

/// Removes the numbered leftovers a retried transfer can leave.
///
/// wget appends `.1`, `.2` and so on rather than overwriting, and upstream
/// deletes `$file*` before retrying. Deleting by glob would also take the file
/// itself, which is what upstream intends; this does the same, but only for
/// names that are the file or the file plus a numeric suffix, so an unrelated
/// name sharing the prefix survives.
///
/// # Errors
/// Returns an error if the directory cannot be read.
pub fn remove_partials(dir: &Path, file_name: &str) -> std::io::Result<usize> {
    let mut removed = 0;
    if !dir.exists() {
        return Ok(removed);
    }
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let Some(name) = entry.file_name().to_str().map(String::from) else {
            continue;
        };
        let matches = name == file_name
            || name
                .strip_prefix(file_name)
                .and_then(|s| s.strip_prefix('.'))
                .is_some_and(|s| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit()));
        if matches && std::fs::remove_file(entry.path()).is_ok() {
            removed += 1;
        }
    }
    Ok(removed)
}

/// Fetches one sample's files, verifying each against the index.
///
/// A file already present with the right checksum is left alone, which is what
/// makes the script safe to re-run over a partly downloaded cohort. The retry
/// budget is shared across the sample and reset by each success, as upstream
/// counts it.
///
/// # Errors
/// Returns an error if the budget runs out, or if the filesystem refuses.
pub fn fetch_sample<F: Fetcher>(
    fetcher: &F,
    root: &Path,
    sample: &str,
    records: &[&Record],
) -> Result<Vec<(String, Outcome)>, DownloadError> {
    let dir = sample_dir(root, sample);
    std::fs::create_dir_all(&dir)?;

    let mut outcomes = Vec::new();
    let mut tries = 0_usize;

    for record in records {
        let name = record.file_name().to_string();
        let path = dir.join(&name);

        if path.exists() && md5_of(&path)? == record.md5 {
            outcomes.push((name, Outcome::AlreadyValid));
            continue;
        }

        let url = https_url(&record.fastq_file);
        let mut attempts = 0;
        loop {
            remove_partials(&dir, &name)?;
            attempts += 1;

            let transferred = fetcher.fetch(&url, &path).is_ok();
            let verified = transferred && path.exists() && md5_of(&path)? == record.md5;
            if verified {
                // A success clears the budget, as upstream's reset does.
                tries = 0;
                outcomes.push((name.clone(), Outcome::Fetched { attempts }));
                break;
            }

            tries += 1;
            if tries >= MAX_TRIES {
                return Err(DownloadError::TooManyTries {
                    sample: sample.to_string(),
                });
            }
        }
    }

    Ok(outcomes)
}

/// A [`Fetcher`] that really downloads.
#[derive(Debug, Clone, Default)]
pub struct HttpFetcher;

impl Fetcher for HttpFetcher {
    fn fetch(&self, url: &str, destination: &Path) -> anyhow::Result<()> {
        let mut response = ureq::get(url).call()?;
        let mut body = response.body_mut().as_reader();
        let mut file = std::fs::File::create(destination)?;
        std::io::copy(&mut body, &mut file)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;

    /// A fetcher that writes what it is told to, in the order it is told.
    struct Scripted {
        /// Bodies to write, one per call; `None` fails the call.
        bodies: RefCell<Vec<Option<Vec<u8>>>>,
        calls: RefCell<usize>,
    }

    impl Scripted {
        fn new(bodies: Vec<Option<&str>>) -> Self {
            Self {
                bodies: RefCell::new(
                    bodies
                        .into_iter()
                        .map(|b| b.map(|s| s.as_bytes().to_vec()))
                        .collect(),
                ),
                calls: RefCell::new(0),
            }
        }
    }

    impl Fetcher for Scripted {
        fn fetch(&self, _url: &str, destination: &Path) -> anyhow::Result<()> {
            let mut calls = self.calls.borrow_mut();
            let body = self.bodies.borrow().get(*calls).cloned().flatten();
            *calls += 1;
            match body {
                Some(bytes) => {
                    std::fs::write(destination, bytes)?;
                    Ok(())
                }
                None => anyhow::bail!("transfer failed"),
            }
        }
    }

    fn record(name: &str, md5: &str) -> Record {
        let mut f: Vec<String> = vec![String::new(); 26];
        f[0] = format!("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR000/{name}");
        f[1] = md5.to_string();
        f[9] = "HG00250".to_string();
        Record::parse(&f.join("\t")).expect("a well-formed row")
    }

    /// The MD5 of "hello\n", which `md5sum` agrees on.
    const HELLO_MD5: &str = "b1946ac92492d2347c6235b4d2611184";

    #[test]
    fn ftp_urls_are_rewritten_to_https() {
        assert_eq!(
            https_url("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/x_1.fastq.gz"),
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/x_1.fastq.gz"
        );
        // Anything else is left alone.
        let already = "https://example.org/x.gz";
        assert_eq!(https_url(already), already);
        assert_eq!(https_url("file:///tmp/x.gz"), "file:///tmp/x.gz");
    }

    #[test]
    fn the_checksum_matches_md5sum() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("f");
        std::fs::write(&path, b"hello\n").unwrap();
        assert_eq!(md5_of(&path).unwrap(), HELLO_MD5);
    }

    /// A file already present and correct is not fetched again, which is what
    /// makes a re-run over a partly downloaded cohort cheap.
    #[test]
    fn an_existing_valid_file_is_left_alone() {
        let dir = tempfile::tempdir().unwrap();
        let sample = sample_dir(dir.path(), "HG00250");
        std::fs::create_dir_all(&sample).unwrap();
        std::fs::write(sample.join("a_1.fastq.gz"), b"hello\n").unwrap();

        let fetcher = Scripted::new(vec![]);
        let r = record("a_1.fastq.gz", HELLO_MD5);
        let out = fetch_sample(&fetcher, dir.path(), "HG00250", &[&r]).unwrap();

        assert_eq!(out, [("a_1.fastq.gz".to_string(), Outcome::AlreadyValid)]);
        assert_eq!(*fetcher.calls.borrow(), 0, "it should not have fetched");
    }

    #[test]
    fn a_clean_fetch_takes_one_attempt() {
        let dir = tempfile::tempdir().unwrap();
        let fetcher = Scripted::new(vec![Some("hello\n")]);
        let r = record("a_1.fastq.gz", HELLO_MD5);
        let out = fetch_sample(&fetcher, dir.path(), "HG00250", &[&r]).unwrap();
        assert_eq!(out[0].1, Outcome::Fetched { attempts: 1 });
    }

    /// A transfer that succeeds but produces the wrong bytes is caught by the
    /// checksum and retried, which is the case a plain retry-on-error misses.
    #[test]
    fn a_corrupt_transfer_is_retried() {
        let dir = tempfile::tempdir().unwrap();
        let fetcher = Scripted::new(vec![Some("wrong\n"), Some("hello\n")]);
        let r = record("a_1.fastq.gz", HELLO_MD5);
        let out = fetch_sample(&fetcher, dir.path(), "HG00250", &[&r]).unwrap();
        assert_eq!(out[0].1, Outcome::Fetched { attempts: 2 });
    }

    #[test]
    fn a_failing_transfer_is_retried() {
        let dir = tempfile::tempdir().unwrap();
        let fetcher = Scripted::new(vec![None, None, Some("hello\n")]);
        let r = record("a_1.fastq.gz", HELLO_MD5);
        let out = fetch_sample(&fetcher, dir.path(), "HG00250", &[&r]).unwrap();
        assert_eq!(out[0].1, Outcome::Fetched { attempts: 3 });
    }

    #[test]
    fn the_budget_runs_out_eventually() {
        let dir = tempfile::tempdir().unwrap();
        let fetcher = Scripted::new(vec![None; MAX_TRIES + 5]);
        let r = record("a_1.fastq.gz", HELLO_MD5);
        let err = fetch_sample(&fetcher, dir.path(), "HG00250", &[&r]).unwrap_err();
        assert!(matches!(err, DownloadError::TooManyTries { .. }));
    }

    /// The budget is shared across the sample and reset by a success, so a
    /// sample with a few flaky files still completes.
    #[test]
    fn a_success_clears_the_budget() {
        let dir = tempfile::tempdir().unwrap();
        // Nine failures spread over three files, none of them consecutive
        // enough to exhaust a reset budget.
        let fetcher = Scripted::new(vec![
            None,
            None,
            None,
            Some("hello\n"),
            None,
            None,
            None,
            Some("hello\n"),
            None,
            None,
            None,
            Some("hello\n"),
        ]);
        let records: Vec<Record> = (1..=3)
            .map(|i| record(&format!("a_{i}.fastq.gz"), HELLO_MD5))
            .collect();
        let refs: Vec<&Record> = records.iter().collect();
        let out = fetch_sample(&fetcher, dir.path(), "HG00250", &refs).unwrap();
        assert_eq!(out.len(), 3);
        assert!(
            out.iter()
                .all(|(_, o)| matches!(o, Outcome::Fetched { .. }))
        );
    }

    /// wget's numbered leftovers go before a retry, and nothing else does.
    #[test]
    fn partials_are_removed_and_neighbours_are_not() {
        let dir = tempfile::tempdir().unwrap();
        for name in [
            "a_1.fastq.gz",
            "a_1.fastq.gz.1",
            "a_1.fastq.gz.12",
            "a_1.fastq.gz.part",
            "a_1.fastq.gz.bak",
            "b_1.fastq.gz",
        ] {
            std::fs::write(dir.path().join(name), b"x").unwrap();
        }

        let removed = remove_partials(dir.path(), "a_1.fastq.gz").unwrap();
        assert_eq!(removed, 3, "the file and its two numbered leftovers");
        assert!(dir.path().join("a_1.fastq.gz.part").exists());
        assert!(dir.path().join("a_1.fastq.gz.bak").exists());
        assert!(dir.path().join("b_1.fastq.gz").exists());
    }

    #[test]
    fn removing_partials_from_a_missing_directory_is_not_an_error() {
        let dir = tempfile::tempdir().unwrap();
        assert_eq!(
            remove_partials(&dir.path().join("nothing-here"), "a").unwrap(),
            0
        );
    }

    #[test]
    fn files_land_under_the_sample_directory() {
        assert_eq!(
            sample_dir(Path::new("/work"), "HG00250"),
            Path::new("/work/fastq/HG00250")
        );
    }
}
