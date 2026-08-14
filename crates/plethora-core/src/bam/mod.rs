//! Reading alignments and turning them into intervals.
//!
//! Replaces the two external tools `make_bed.sh` shells out to before any
//! interval arithmetic happens: `samtools sort -n` and `bedtools bamtobed`.
//! The ordering those tools impose reaches the final coverage numbers, so both
//! are reproduced rather than approximated; see [`plethora_compat::strnum`] for
//! the sort order and [`bamtobed`] for the conversion rules.

pub mod bamtobed;
pub mod namesort;
pub mod reader;
