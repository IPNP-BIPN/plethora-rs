//! The 1000 Genomes chain: everything upstream keeps in `code/1000genomes`.
//!
//! These scripts fetch and manage the public data the paper was validated on.
//! They are separate from single-sample processing, which is why the breakages
//! in them went unnoticed for years: a run on one's own cohort never touches
//! any of it.

pub mod align_report;
pub mod clean;
pub mod download;
pub mod sample_index;
