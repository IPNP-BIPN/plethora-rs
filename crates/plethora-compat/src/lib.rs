//! Bit-exact ports of the third-party behaviour plethora's numbers depend on.
//!
//! The upstream pipeline is Perl and R calling out to GNU coreutils, samtools
//! and bedtools. Reproducing its output is not a matter of reimplementing the
//! documented algorithms: it means reproducing the specific implementations,
//! including the parts nobody would design on purpose. A tie in GNU `sort` is
//! broken by comparing the whole line; `phrtsd` drops the last character of its
//! phrase; R's `loess` returns interpolated values rather than fitted ones.
//!
//! Everything in this crate exists because one of those details moves a digit
//! in the final copy-number table. Each module names the exact source it was
//! transcribed from, and each is pinned by golden vectors generated from the
//! real tool (see `tests/oracle/`).

pub mod gnusort;
pub mod randlib;
pub mod rmath;
pub mod strnum;
