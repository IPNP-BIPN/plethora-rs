//! GC content, the model built from it, and the correction it drives.
//!
//! Three upstream pieces: `gc_from_fasta.pl` measures GC per domain,
//! `build_gc_model.sh` extracts the sequences to measure, and
//! `gc_correction.R` turns coverage into copy number using the result.
//!
//! The whole calibration rests on regions assumed to sit at diploid copy
//! number. They define the correction curve and set the haploid unit, so the
//! reported copy numbers are relative to them.

pub mod correction;
pub mod from_fasta;
