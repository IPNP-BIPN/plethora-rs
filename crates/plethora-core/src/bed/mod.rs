//! Interval arithmetic: the three bedtools calls in `make_bed.sh`, plus the
//! GNU sort between them.
//!
//! ```text
//! sort -k 1,1 -k 2,2n -T ./ ${output}_edited.bed > ${output}_sorted.bed
//! bedtools intersect -wao -sorted -a $reference_bed -b ${output}_sorted.bed > ${output}_temp.bed
//! awk 'OFS="\t" {print $4,$2,$3,$1,$13}' ${output}_temp.bed | bedtools merge -c 5 -o sum -i -
//! ```
//!
//! Written out rather than taken from an interval crate. The requirement is not
//! "compute overlaps" but "compute them the way bedtools 2.31 does", down to the
//! null record for an uncovered domain, the refusal to count book-ended
//! intervals as overlapping, and the error raised when a key reappears out of
//! order. No general-purpose library commits to those.

pub mod intersect;
pub mod merge;
pub mod sort;
