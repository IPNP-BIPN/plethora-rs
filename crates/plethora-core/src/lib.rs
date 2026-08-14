//! The plethora pipeline, stage by stage.
//!
//! Each module replaces one step of the upstream shell, Perl and R scripts.
//! Where a stage's output depends on the exact behaviour of a third-party tool,
//! that behaviour lives in `plethora-compat` and is pinned by golden vectors;
//! this crate holds the pipeline logic that calls it.

pub mod bam;
pub mod batch;
pub mod bed;
pub mod config;
pub mod coverage;
pub mod gc;
pub mod io;
pub mod merge_pairs;
pub mod onekg;
pub mod trim;
