//! `plethora.toml`, which replaces `config.sh`.
//!
//! Upstream keeps a project's settings in a shell script that every batch step
//! sources:
//!
//! ```text
//! sample_index=data/1000Genomes_samples.txt
//! genome=$HOME/genomes/bowtie2.2.9_indicies/hg38/hg38
//! master_ref=data/hg38_duf_full_domains_v2.3.bed
//! bowtie_params=""
//! alignment_dir=alignments
//! bed_dir=results
//! SAMPLES=( NA19914 HG00623 ... )
//! ```
//!
//! Sourcing it is what makes the sample list reachable from a job array, and
//! also what makes it unreadable to anything but bash: `trim_qc_report.R` has to
//! slice the array out of the file with `readLines` and a pair of `grep`s. A
//! declarative file is readable by every stage without pretending to be a shell.
//!
//! `config.sh` also creates its output directories as a side effect of being
//! sourced. That happens here too, but only when asked, through
//! [`Config::create_directories`].

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::io::Compress;

/// A project's settings.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    /// Where the reference and the aligner index live.
    pub reference: Reference,
    /// Where the pipeline reads and writes.
    #[serde(default)]
    pub paths: Paths,
    /// How reads are processed.
    #[serde(default)]
    pub options: Options,
    /// The samples to run, in order. A job array indexes into this.
    #[serde(default)]
    pub samples: Vec<String>,
}

/// The reference data a run needs.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Reference {
    /// The BED of domains to measure. `master_ref` upstream.
    pub domains: PathBuf,
    /// The per-domain GC table the correction joins against.
    ///
    /// Upstream derives this from `master_ref` by substituting `_GC.txt` for
    /// `.bed`, which is why the two must be named in step. Named here instead,
    /// so a table built elsewhere can be used without renaming the BED.
    pub gc_table: PathBuf,
    /// The aligner index prefix. `genome` upstream.
    #[serde(default)]
    pub index: Option<PathBuf>,
    /// The genome FASTA, needed only to rebuild the GC table.
    #[serde(default)]
    pub fasta: Option<PathBuf>,
    /// The 1000 Genomes sequence index, for the download and selection steps.
    #[serde(default)]
    pub sample_index: Option<PathBuf>,
}

/// Where things go. The defaults are upstream's.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Paths {
    pub fastq: PathBuf,
    /// `alignment_dir` upstream.
    pub alignments: PathBuf,
    /// `bed_dir` upstream.
    pub results: PathBuf,
    pub logs: PathBuf,
    /// Where the sorts spill. Upstream passes `-T ./`, so the working
    /// directory, which is rarely what anyone wants on a cluster.
    pub temp: PathBuf,
}

impl Default for Paths {
    fn default() -> Self {
        Self {
            fastq: PathBuf::from("fastq"),
            alignments: PathBuf::from("alignments"),
            results: PathBuf::from("results"),
            logs: PathBuf::from("logs"),
            temp: PathBuf::from("."),
        }
    }
}

/// How reads are processed.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Options {
    /// Reads are paired, which is what the paper's data is.
    #[serde(default = "default_true")]
    pub paired: bool,
    /// Quality cutoff for trimming. Upstream's `-q 10`.
    #[serde(default = "default_quality")]
    pub quality_cutoff: u8,
    /// Minimum length after trimming. Upstream's `--minimum-length 80`.
    #[serde(default = "default_min_length")]
    pub min_length: usize,
    /// Whether outputs are gzipped.
    #[serde(default)]
    pub compress: bool,
    /// Threads per sample, passed to the aligner.
    #[serde(default = "default_threads")]
    pub threads: usize,
    /// Samples processed at once by the local runner.
    #[serde(default = "default_jobs")]
    pub jobs: usize,
    /// Extra arguments for the aligner. `bowtie_params` upstream.
    #[serde(default)]
    pub aligner_args: Vec<String>,
}

const fn default_true() -> bool {
    true
}
const fn default_quality() -> u8 {
    10
}
const fn default_min_length() -> usize {
    80
}
const fn default_threads() -> usize {
    12
}
const fn default_jobs() -> usize {
    1
}

impl Default for Options {
    fn default() -> Self {
        Self {
            paired: true,
            quality_cutoff: default_quality(),
            min_length: default_min_length(),
            compress: false,
            threads: default_threads(),
            jobs: default_jobs(),
            aligner_args: Vec::new(),
        }
    }
}

/// A configuration that cannot be run.
#[derive(Debug, thiserror::Error)]
pub enum ConfigError {
    #[error("cannot read {path}: {source}")]
    Read {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("cannot parse {path}: {source}")]
    Parse {
        path: PathBuf,
        #[source]
        source: toml::de::Error,
    },

    #[error("the sample list has {count} duplicate(s), starting with {first}")]
    DuplicateSamples { count: usize, first: String },

    #[error("{what} is not readable: {path}")]
    Missing { what: &'static str, path: PathBuf },
}

impl Config {
    /// Reads a `plethora.toml`.
    ///
    /// # Errors
    /// Returns an error if the file cannot be read or parsed, or if the sample
    /// list repeats a name.
    pub fn load(path: &Path) -> Result<Self, ConfigError> {
        let text = std::fs::read_to_string(path).map_err(|source| ConfigError::Read {
            path: path.to_path_buf(),
            source,
        })?;
        let config: Self = toml::from_str(&text).map_err(|source| ConfigError::Parse {
            path: path.to_path_buf(),
            source,
        })?;
        config.check_samples()?;
        Ok(config)
    }

    /// Rejects a repeated sample name.
    ///
    /// Upstream's list carries duplicates: `HG00261` and `HG00270` each appear
    /// twice in `config.sh`, which means two job array entries process the same
    /// sample into the same output files at the same time. Worth refusing
    /// rather than racing.
    ///
    /// # Errors
    /// Returns an error naming the first repeat.
    pub fn check_samples(&self) -> Result<(), ConfigError> {
        let mut seen: HashSet<&str> = HashSet::new();
        let mut duplicates: Vec<&str> = Vec::new();
        for sample in &self.samples {
            if !seen.insert(sample.as_str()) {
                duplicates.push(sample.as_str());
            }
        }
        if duplicates.is_empty() {
            return Ok(());
        }
        Err(ConfigError::DuplicateSamples {
            count: duplicates.len(),
            first: duplicates[0].to_string(),
        })
    }

    /// Checks that the reference files exist.
    ///
    /// # Errors
    /// Returns an error naming the first file that is not readable.
    pub fn check_reference(&self) -> Result<(), ConfigError> {
        for (what, path) in [
            ("the domain BED", &self.reference.domains),
            ("the GC table", &self.reference.gc_table),
        ] {
            if !path.exists() {
                return Err(ConfigError::Missing {
                    what,
                    path: path.clone(),
                });
            }
        }
        Ok(())
    }

    /// Creates the output directories, as sourcing `config.sh` does.
    ///
    /// # Errors
    /// Returns an error if a directory cannot be created.
    pub fn create_directories(&self) -> std::io::Result<()> {
        for dir in [
            &self.paths.alignments,
            &self.paths.results,
            &self.paths.logs,
        ] {
            std::fs::create_dir_all(dir)?;
        }
        Ok(())
    }

    /// Whether outputs are compressed.
    #[must_use]
    pub const fn compress(&self) -> Compress {
        if self.options.compress {
            Compress::Yes
        } else {
            Compress::No
        }
    }

    /// The output prefix for a sample: `results/<sample>`.
    #[must_use]
    pub fn result_prefix(&self, sample: &str) -> PathBuf {
        self.paths.results.join(sample)
    }

    /// The alignment for a sample: `alignments/<sample>.bam`.
    #[must_use]
    pub fn alignment(&self, sample: &str) -> PathBuf {
        self.paths.alignments.join(format!("{sample}.bam"))
    }

    /// Where a sample's reads live: `fastq/<sample>/`.
    #[must_use]
    pub fn fastq_dir(&self, sample: &str) -> PathBuf {
        self.paths.fastq.join(sample)
    }

    /// Serialises back to TOML, for `plethora config init`.
    ///
    /// # Errors
    /// Returns an error if the value cannot be serialised.
    pub fn to_toml(&self) -> Result<String, toml::ser::Error> {
        toml::to_string_pretty(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MINIMAL: &str = r#"
[reference]
domains = "data/hg38_duf_full_domains_v2.3.bed"
gc_table = "data/hg38_duf_full_domains_v2.3_GC.txt"
"#;

    #[test]
    fn a_minimal_config_takes_upstreams_defaults() {
        let c: Config = toml::from_str(MINIMAL).unwrap();
        assert_eq!(c.paths.alignments, Path::new("alignments"));
        assert_eq!(c.paths.results, Path::new("results"));
        assert_eq!(c.paths.fastq, Path::new("fastq"));
        assert!(c.options.paired);
        assert_eq!(c.options.quality_cutoff, 10, "upstream's -q 10");
        assert_eq!(c.options.min_length, 80, "upstream's --minimum-length 80");
        assert_eq!(c.options.threads, 12, "upstream's bowtie2 -p 12");
        assert!(c.samples.is_empty());
    }

    /// A misspelled key is refused rather than silently ignored, which is the
    /// failure mode a shell config cannot avoid.
    #[test]
    fn an_unknown_key_is_refused() {
        let text = format!("{MINIMAL}\n[options]\nqualty_cutoff = 20\n");
        let err = toml::from_str::<Config>(&text).unwrap_err();
        assert!(
            err.to_string().contains("qualty_cutoff"),
            "the message should name the key: {err}"
        );
    }

    /// Upstream's own list repeats HG00261 and HG00270, which would have two
    /// job array entries writing the same outputs at once.
    #[test]
    fn a_repeated_sample_is_refused() {
        // Top-level keys go before the tables, or TOML reads them as part of
        // the last table opened.
        let text = format!("samples = [\"HG00261\", \"HG00270\", \"HG00261\"]\n{MINIMAL}");
        let c: Config = toml::from_str(&text).unwrap();
        let err = c.check_samples().unwrap_err();
        assert!(matches!(
            err,
            ConfigError::DuplicateSamples { count: 1, ref first } if first == "HG00261"
        ));
    }

    #[test]
    fn a_list_without_repeats_passes() {
        let text = format!("samples = [\"A\", \"B\", \"C\"]\n{MINIMAL}");
        let c: Config = toml::from_str(&text).unwrap();
        assert!(c.check_samples().is_ok());
        assert_eq!(c.samples.len(), 3);
    }

    #[test]
    fn paths_are_derived_the_way_upstream_names_them() {
        let c: Config = toml::from_str(MINIMAL).unwrap();
        assert_eq!(c.result_prefix("HG00250"), Path::new("results/HG00250"));
        assert_eq!(c.alignment("HG00250"), Path::new("alignments/HG00250.bam"));
        assert_eq!(c.fastq_dir("HG00250"), Path::new("fastq/HG00250"));
    }

    #[test]
    fn compression_follows_the_option() {
        let mut c: Config = toml::from_str(MINIMAL).unwrap();
        assert_eq!(c.compress(), Compress::No);
        c.options.compress = true;
        assert_eq!(c.compress(), Compress::Yes);
    }

    #[test]
    fn a_missing_reference_is_named() {
        let c: Config = toml::from_str(MINIMAL).unwrap();
        let err = c.check_reference().unwrap_err();
        assert!(
            err.to_string().contains("domain BED"),
            "the message should say which file: {err}"
        );
    }

    #[test]
    fn directories_are_created_only_when_asked() {
        let dir = tempfile::tempdir().unwrap();
        let text = format!(
            "{MINIMAL}\n[paths]\nfastq = \"{d}/fastq\"\nalignments = \"{d}/aln\"\nresults = \"{d}/res\"\nlogs = \"{d}/logs\"\ntemp = \"{d}\"\n",
            d = dir.path().display()
        );
        let c: Config = toml::from_str(&text).unwrap();
        assert!(!dir.path().join("aln").exists(), "not created by parsing");

        c.create_directories().unwrap();
        assert!(dir.path().join("aln").is_dir());
        assert!(dir.path().join("res").is_dir());
        assert!(dir.path().join("logs").is_dir());
    }

    #[test]
    fn a_config_round_trips_through_toml() {
        let c: Config = toml::from_str(MINIMAL).unwrap();
        let text = c.to_toml().unwrap();
        let back: Config = toml::from_str(&text).unwrap();
        assert_eq!(back.reference.domains, c.reference.domains);
        assert_eq!(back.options.min_length, c.options.min_length);
    }

    #[test]
    fn loading_names_the_file_it_could_not_read() {
        let err = Config::load(Path::new("/nonexistent/plethora.toml")).unwrap_err();
        assert!(err.to_string().contains("plethora.toml"), "{err}");
    }
}
