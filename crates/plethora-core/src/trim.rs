//! Read trimming, on Trim Galore rather than cutadapt.
//!
//! Upstream runs:
//!
//! ```text
//! cutadapt -a XXX -A XXX -q 10 --minimum-length 80 --trim-n \
//!  -o fastq/test_1_filtered.fastq.gz -p fastq/test_2_filtered.fastq.gz \
//!  fastq/test_1.fastq.gz fastq/test_2.fastq.gz
//! ```
//!
//! The adapter is the literal string `XXX`. X is not a nucleotide, so it
//! cannot match anything: measured on upstream's own test data, 0 of 135 read
//! pairs had an adapter found. What cutadapt actually does in this pipeline is
//! quality trimming at Q10, `--trim-n`, and a paired length filter at 80 bp.
//! Adapter contamination is left in the reads.
//!
//! This module keeps the parameters that matter to the pipeline and fixes the
//! part that does not work. Adapters are detected from the data by Trim Galore
//! and actually removed; the quality cutoff stays at upstream's 10, and so do
//! the 80 bp floor and `--trim-n`.
//!
//! That is a deliberate divergence, and it moves reads: a pair whose adapter is
//! now removed aligns differently, so copy numbers shift slightly against the
//! published figures. Recorded in `DIVERGENCES.md`.
//!
//! ## Licence
//!
//! `trim-galore` is GPL-3.0-only, so Plethora-rs is GPL-3.0-only. The upstream
//! plethora scripts are MIT, which is compatible in this direction: MIT code
//! may be incorporated into a GPL work.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use trim_galore::adapter::{self, AdapterPreset, DetectionResult};
use trim_galore::filters::UnpairedLengths;
use trim_galore::trimmer::{self, TrimConfig};

/// Quality cutoff, from upstream's `-q 10`.
pub const QUALITY_CUTOFF: u8 = 10;
/// Minimum length for a surviving pair, from upstream's `--minimum-length 80`.
///
/// Not a detail: the domains being measured are 1000 bp repeats that differ by
/// a handful of bases, so a short read has nowhere unique to go and only adds
/// ambiguity.
pub const MIN_LENGTH: usize = 80;
/// Phred+33, which is every Illumina run since 1.8.
pub const PHRED_OFFSET: u8 = 33;

/// How adapters are chosen.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum Adapters {
    /// Detect from the first reads of the file, as `trim_galore` does by
    /// default. What this pipeline uses.
    #[default]
    Detect,
    /// A named preset, when the library preparation is already known.
    Preset(&'static str),
    /// An explicit sequence.
    Sequence(Vec<u8>),
    /// None at all, reproducing upstream's dummy `XXX` adapter.
    ///
    /// Kept so the published behaviour can be reproduced exactly when
    /// comparing against an existing run.
    None,
}

/// What one trimming run did.
#[derive(Debug, Clone)]
pub struct TrimSummary {
    /// Pairs read.
    pub pairs_in: usize,
    /// Pairs surviving both filters.
    pub pairs_out: usize,
    /// Which adapter was used, and how it was chosen.
    pub adapter: String,
    /// The detector's own report, when detection ran.
    pub detection: Option<String>,
    pub output_r1: PathBuf,
    pub output_r2: PathBuf,
}

/// A named adapter sequence, in the shape `TrimConfig` wants.
pub type NamedAdapter = (String, Vec<u8>);

/// Resolves the adapter to trim with, and the detector's report when it ran.
///
/// # Errors
/// Returns an error if detection cannot read the input.
pub fn resolve_adapters(
    choice: &Adapters,
    r1: &Path,
) -> Result<(Vec<NamedAdapter>, Option<DetectionResult>)> {
    match choice {
        Adapters::None => Ok((Vec::new(), None)),
        Adapters::Sequence(seq) => Ok((vec![("user".to_string(), seq.clone())], None)),
        Adapters::Preset(name) => {
            let preset = preset_by_name(name)
                .with_context(|| format!("unknown adapter preset: {name}"))?;
            Ok((vec![(preset.name.to_string(), preset.seq.as_bytes().to_vec())], None))
        }
        Adapters::Detect => {
            let result = adapter::autodetect_adapter(r1, None)
                .with_context(|| format!("detecting the adapter in {}", r1.display()))?;
            let named = vec![(
                result.adapter.name.to_string(),
                result.adapter.seq.as_bytes().to_vec(),
            )];
            Ok((named, Some(result)))
        }
    }
}

/// The presets Trim Galore ships, by name.
#[must_use]
pub fn preset_by_name(name: &str) -> Option<AdapterPreset> {
    match name.to_ascii_lowercase().as_str() {
        "illumina" => Some(adapter::ILLUMINA),
        "nextera" => Some(adapter::NEXTERA),
        "small_rna" | "smallrna" => Some(adapter::SMALL_RNA),
        "stranded_illumina" => Some(adapter::STRANDED_ILLUMINA),
        "bgiseq" => Some(adapter::BGISEQ),
        _ => None,
    }
}

/// The trimming configuration this pipeline uses.
///
/// `TrimConfig` has no `Default`, so every field is set here. The four at the
/// top come from upstream's cutadapt line; the rest are Trim Galore's own CLI
/// defaults, written out so a reader can see that nothing was quietly changed.
#[must_use]
pub fn config(adapters: Vec<NamedAdapter>) -> TrimConfig {
    TrimConfig {
        // From upstream.
        adapters,
        quality_cutoff: QUALITY_CUTOFF,
        length_cutoff: MIN_LENGTH,
        trim_n: true,

        // Fixed by the pipeline.
        phred_offset: PHRED_OFFSET,
        is_paired: true,

        // Trim Galore's defaults, unchanged.
        adapters_r2: Vec::new(),
        times: 1,
        error_rate: 0.1,
        min_overlap: 1,
        max_length: None,
        max_n: None,
        clip_r1: None,
        clip_r2: None,
        three_prime_clip_r1: None,
        three_prime_clip_r2: None,
        rename: false,
        nextseq: false,
        rrbs: false,
        non_directional: false,
        poly_a: false,
        poly_g: false,
        discard_untrimmed: false,
        gzip_level: trim_galore::fastq::DEFAULT_GZIP_LEVEL,
    }
}

/// Trims a pair of FASTQ files, writing `*_1_filtered.fastq.gz` and
/// `*_2_filtered.fastq.gz` beside them.
///
/// The output names follow upstream's, which the batch scripts and
/// `clean_files.pl` both match on:
///
/// ```text
/// first_filtered=`echo $first_read | sed 's/_1.fastq/_1_filtered.fastq/'`
/// ```
///
/// # Errors
/// Returns an error if either input cannot be read, an output cannot be
/// written, or adapter detection fails.
pub fn trim_pair(r1: &Path, r2: &Path, choice: &Adapters) -> Result<TrimSummary> {
    let (adapters, detection) = resolve_adapters(choice, r1)?;
    let adapter_name = adapters
        .first()
        .map_or_else(|| "none".to_string(), |(name, _)| name.clone());
    let config = config(adapters);

    let out_r1 = filtered_name(r1, 1)?;
    let out_r2 = filtered_name(r2, 2)?;

    let mut reader_r1 = trim_galore::fastq::FastqReader::open(r1)
        .with_context(|| format!("opening {}", r1.display()))?;
    let mut reader_r2 = trim_galore::fastq::FastqReader::open(r2)
        .with_context(|| format!("opening {}", r2.display()))?;

    // Compression follows the output name, as upstream's does: the batch
    // scripts write .fastq.gz and glob for it.
    let gzip = out_r1.extension().is_some_and(|e| e == "gz");
    let mut writer_r1 = trim_galore::fastq::FastqWriter::create(&out_r1, gzip, 1, config.gzip_level)
        .with_context(|| format!("creating {}", out_r1.display()))?;
    let mut writer_r2 = trim_galore::fastq::FastqWriter::create(&out_r2, gzip, 1, config.gzip_level)
        .with_context(|| format!("creating {}", out_r2.display()))?;

    let (stats_r1, _stats_r2, _pairs) = trimmer::run_paired_end(
        &mut reader_r1,
        &mut reader_r2,
        None,
        &mut writer_r1,
        &mut writer_r2,
        None,
        None,
        None,
        &config,
        // Neither mate is kept on its own: upstream discards the whole pair
        // when either read falls below the length floor, and the aligner is
        // run in paired mode.
        UnpairedLengths { r1: 0, r2: 0 },
    )?;

    Ok(TrimSummary {
        pairs_in: stats_r1.total_reads,
        pairs_out: stats_r1.reads_written,
        adapter: adapter_name,
        detection: detection.map(|d| d.message),
        output_r1: out_r1,
        output_r2: out_r2,
    })
}

/// `_1.fastq[.gz]` becomes `_1_filtered.fastq[.gz]`, as upstream's sed does.
///
/// # Errors
/// Returns an error if the name does not carry the expected mate marker, since
/// the batch scripts downstream glob for the filtered spelling and would
/// silently find nothing.
pub fn filtered_name(path: &Path, mate: u8) -> Result<PathBuf> {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .context("the input path has no file name")?;

    let marker = format!("_{mate}.fastq");
    let replacement = format!("_{mate}_filtered.fastq");
    anyhow::ensure!(
        name.contains(&marker),
        "expected {name} to contain {marker}, which is what the pipeline's \
         globs match on"
    );

    Ok(path.with_file_name(name.replacen(&marker, &replacement, 1)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_names_follow_upstreams_sed() {
        let p = Path::new("fastq/HG00250/ERR001_1.fastq.gz");
        assert_eq!(
            filtered_name(p, 1).unwrap(),
            Path::new("fastq/HG00250/ERR001_1_filtered.fastq.gz")
        );
        let p = Path::new("fastq/test_2.fastq");
        assert_eq!(
            filtered_name(p, 2).unwrap(),
            Path::new("fastq/test_2_filtered.fastq")
        );
    }

    /// A name without the mate marker is refused rather than silently passed
    /// through: the batch scripts glob for `*_1_filtered.fastq.gz` and would
    /// find nothing.
    #[test]
    fn a_name_without_the_mate_marker_is_refused() {
        assert!(filtered_name(Path::new("reads.fastq.gz"), 1).is_err());
    }

    /// The three parameters that come from upstream, and nothing else moved.
    #[test]
    fn the_config_keeps_upstreams_parameters() {
        let c = config(vec![("illumina".into(), b"AGATCGGAAGAGC".to_vec())]);
        assert_eq!(c.quality_cutoff, 10, "upstream's -q 10");
        assert_eq!(c.length_cutoff, 80, "upstream's --minimum-length 80");
        assert!(c.trim_n, "upstream's --trim-n");
        assert!(c.is_paired);
        assert_eq!(c.phred_offset, 33);
    }

    /// Reproducing upstream exactly means no adapter at all, since XXX cannot
    /// match a nucleotide.
    #[test]
    fn the_upstream_choice_trims_no_adapter() {
        let (adapters, detection) =
            resolve_adapters(&Adapters::None, Path::new("/nonexistent")).unwrap();
        assert!(adapters.is_empty(), "XXX matches nothing, so nothing is configured");
        assert!(detection.is_none());
    }

    #[test]
    fn presets_resolve_by_name() {
        assert!(preset_by_name("illumina").is_some());
        assert!(preset_by_name("Nextera").is_some(), "case-insensitive");
        assert!(preset_by_name("small_rna").is_some());
        assert!(preset_by_name("nonesuch").is_none());
    }

    #[test]
    fn an_explicit_sequence_is_used_as_given() {
        let (adapters, _) =
            resolve_adapters(&Adapters::Sequence(b"ACGTACGT".to_vec()), Path::new("/nonexistent"))
                .unwrap();
        assert_eq!(adapters.len(), 1);
        assert_eq!(adapters[0].1, b"ACGTACGT");
    }

    #[test]
    fn detection_is_the_default() {
        assert_eq!(Adapters::default(), Adapters::Detect);
    }
}
