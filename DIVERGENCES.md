# Divergences from upstream

Plethora-rs reproduces [dpastling/plethora](https://github.com/dpastling/plethora),
the pipeline behind Astling et al., *"High resolution measurement of DUF1220
domain copy number from whole genome sequence data"*, BMC Genomics 18:614 (2017).

Reproducing it means reproducing its specific implementations, quirks included.
This file records every place the port does something else, and why. Anything not
listed here is intended to match upstream byte for byte.

## What "byte-identical" means here, and where it stops

Each stage is pinned by golden vectors or differential tests against the real
tool. Where a residual difference exists it is measured, not waved away.

| Stage | Oracle | Result |
|---|---|---|
| `merge_pairs.pl` | the Perl itself, with `Math::Random` | byte-identical, 3000 lines, 1058 random draws |
| `samtools sort -n` | samtools 1.24 | identical order, 1216 records |
| `bedtools bamtobed` | bedtools 2.31.1 | byte-identical, 800 records, BED and BEDPE |
| `sort -k1,1 -k2,2n` | GNU coreutils `sort` | byte-identical, 4000 lines |
| `bedtools intersect -wao -sorted` | bedtools 2.31.1 | byte-identical, 5185 lines |
| `bedtools merge -c 5 -o sum` | bedtools 2.31.1 | byte-identical, 360 rows |
| `awk` read depth | awk `OFMT` | byte-identical, 360 rows |
| `gc_from_fasta.pl` | the Perl itself | byte-identical, 200 sequences |
| `bedtools getfasta` | bedtools 2.31.1 | byte-identical, 60 domains |
| `gc_correction.R` | R 4.6.1 with dplyr | see below |

### The one numeric residual: BLAS inside `loess`

`gc_correction.R` fits an R `loess`, and `loess` calls LINPACK, which calls
BLAS. Which BLAS is part of the answer.

This port transcribes the reference BLAS that R ships in `src/extra/blas`, which
is what a stock R build uses. An R linked against an optimised BLAS gives
slightly different results: the R used for these tests links OpenBLAS 0.3.34,
whose NEON dot product accumulates into several partial sums and lands 8 ULP
from sequential summation on a 39-element product, which is exactly the
neighbourhood size `ehg127` uses.

Measured against a real 623,699-domain sample from a 394-genome cohort processed
with the upstream scripts:

- `percent.gc` agrees **as text on every one of 623,699 domains**
- the numeric columns agree to **2.0e-14** relative

For scale: running the upstream R script itself on this machine, against the
same sample's output produced on the cohort's machine, disagrees on 103,907 of
623,699 rows, always at the last digit. R does not reproduce itself across
machines any more closely than this port reproduces R.

## Deliberate divergences

### Trimming: Trim Galore instead of cutadapt

Upstream runs

```
cutadapt -a XXX -A XXX -q 10 --minimum-length 80 --trim-n ...
```

The adapter is the literal string `XXX`. X is not a nucleotide, so it cannot
match: on upstream's own test data, 0 of 135 read pairs had an adapter found and
cutadapt removed 11 bases in total. **Adapter contamination is left in the
reads.**

Plethora-rs keeps the quality cutoff of 10, the 80 bp length floor and
`--trim-n`, and detects adapters from the data instead of not trimming them.

Measured on those same 135 pairs:

- with `Adapters::None`, which reproduces the dummy `XXX`, output matches
  cutadapt exactly: 0 reads differ, 0 bases
- with `Adapters::Detect`, the detector finds no adapter, falls back to the
  Illumina preset, and removes 60 bases across 39 reads

Trim Galore's default stringency is 1, so one overlapping base with the fallback
adapter is enough to trim. On a real library that is right; on simulated reads
it removes signal for nothing. `Adapters::None` exists so an existing run can
still be reproduced.

### Alignment: BWA-MEM instead of bowtie2

Upstream aligns with `bowtie2 --very-sensitive`, end to end. Plethora-rs uses
[bwa-mem4](https://github.com/IPNP-BIPN/bwa-mem4).

BWA-MEM is a local aligner: it soft-clips, emits supplementary alignments, and
selects a primary among equal-scoring hits differently. On a locus with more
than 300 near-identical paralogues, that changes how multi-mapping reads are
distributed, so **copy numbers from this path are not the paper's**.

`plethora coverage --bam` accepts a BAM produced anywhere, so the bowtie2 recipe
remains available for paper parity:

```
bowtie2 -p 12 --very-sensitive --minins 0 --maxins 800 -x hg38 -1 R1.fq.gz -2 R2.fq.gz \
  | samtools view -Sb - > sample.bam
```

### Supplementary and secondary alignments are filtered

`bedtools bamtobed -bedpe` consumes records strictly two at a time and pairs
them only when their names match; on a mismatch it scans forward, warning per
skipped record. A third record under one name, which is what BWA-MEM emits for a
supplementary alignment, pushes every later record out of phase and makes
bedtools drop pairs in a cascade.

Plethora-rs filters flags `0x100` and `0x800` before pairing. bedtools does not.
Under bowtie2 the situation never arises, which is why the paper's pipeline
never met it; under BWA-MEM it always does.

### Name sorting follows samtools 1.20, and can follow the older rule

`samtools sort -n` orders on `strnum_cmp`, which reads runs of digits as
numbers. When two records carry the same name, the mates are ordered between
themselves, and **samtools changed how in 1.20**: before that release the
comparison fell back to the flags, after it the READ1/READ2 bits decide
directly. The two rules disagree on real files, so a name-sorted BAM is not one
thing to be identical to.

Plethora-rs implements both, as `TieBreak::Since1_20` and
`TieBreak::Before1_20`, and defaults to the current one. The differential tests
ask the installed `samtools --version` which rule to expect, which is why they
pass on this machine against 1.22 and on CI against 1.19.2. This surfaced as a
CI failure that looked like a port bug and was not one.

Downstream this only reorders records within a name, and `merge_pairs` reads a
pair as a unit, so the coverage figures do not move. It matters if you diff a
name-sorted BAM against one from another samtools.

### Two different counts, one index column

`READ_COUNT` is per file, so a paired sample's rows hold it twice. The two
consumers want different things from that, and getting either backwards makes a
whole stage silently unreachable rather than loudly wrong.

| Consumer | Expectation | Compared against |
|---|---|---|
| `clean_files.pl` | `sum(READ_COUNT)`, both mates | FASTQ records over both files, and BAM records |
| `trim_qc_report.R` | `sum(READ_COUNT) / 2`, pairs | one trimming-log row per pair |

Upstream says the first out loud where it does the sum: *"note both pairs will
be present in alignment file, so we need to count both"*. The second follows
from the R comparing `expected.files = n() / 2` against
`n.files = sum(type == "total")`. Only the BED is per fragment, and the cleanup
chain halves it there.

### Trimming says when nothing survives

Upstream prints cutadapt's two counts and stops there. A 1000 Genomes run whose
second mate is entirely Q2, which is Illumina's marker for a failed read
segment, trims to nothing under `-q 10`, and every pair then fails
`--minimum-length 80`. That is the right answer, and it looks exactly like a
broken trimmer. Plethora-rs names it, and warns below half survival too. Only
the message is new; the counts are the same.

### `trim_qc_report` does not delete files by default

`trim_qc_report.R` calls `cleanup_old_files()` unconditionally, so merely running
the report deletes FASTQ files. Here deletion requires an explicit flag.

### Output ordering of `gc_from_fasta`

The Perl iterates a hash, so its output order is whatever that run's hash
randomisation produces; the same input gives a differently ordered file each
time. Plethora-rs emits in FASTA order. Values are identical; only the order
differs.

## Upstream quirks kept on purpose

These look like bugs. They are, and they produced the published numbers, so they
are reproduced rather than fixed. Correcting any of them silently would rescale
results.

### The read-depth divisor is the domain length plus one

```
awk 'OFS="\t" { print $1, $4 / ($3 - $2 + 1)}'
```

A 1000 bp domain is divided by 1001. Every coverage figure is therefore 0.1%
low. Kept.

### Lowercase and ambiguous bases count as GC

`gc_from_fasta.pl` removes `[ATN]` and calls what is left GC. The substitution is
case-sensitive, so soft-masked repeat sequence, which `bedtools getfasta` writes
in lowercase, counts as GC in full: `acgt` contributes four GC bases, not two.
So do the IUPAC ambiguity codes. Uppercase `N` lowers the figure while lowercase
`n` raises it. Kept; `GcCounts` reports the soft-masked and ambiguous counts
separately so the effect is visible.

### Two different patterns select the conserved regions

`gc_correction.R` selects the domains that define the GC curve with
`^((baseline)|(uc))`, anchored, and the ones that set the haploid unit with
`((baseline)|(uc))`, unanchored. A domain merely containing `uc` therefore moves
the ploidy normalisation without touching the curve. Kept.

### `merge_pairs.pl` applies two different ceilings

`$max_inner_distance` is a file-scope global. The first pass classifies pairs
against 800; the second against `mean + 5 * sd`, measured during the first. The
two passes do not apply the same test, and the first ceiling is what selects
which pairs define the distribution at all. Kept.

### Fragment statistics stop at 50 million records

`last if ($#distance > $sufficient_number_of_reads)` compares the last index, not
the length, so collection stops at 50,000,002 entries. A whole-genome sample
reaches that, so the mean and standard deviation depend on which records came
first, and therefore on the exact `samtools sort -n` order. Kept, which is why
that sort order is reproduced exactly.

### Reads can be extended past the end of a chromosome

`merge_pairs.pl` clamps an extended read's start at zero and does not clamp its
end. The intersect stage drops the overhang, so upstream never noticed. Kept.

## Upstream breakages found, and reported

Five verified problems in upstream, each with a pull request. None of them
affects the core path a typical run takes, which is why the pipeline has worked
for years: the GC table ships prebuilt, so `build_gc_model.sh` is never run, and
the 1000 Genomes scripts are separate from single-sample processing.

| Problem | PR |
|---|---|
| `as.tbl()` is defunct from dplyr 1.2.0, so `gc_correction.R` stops before reading a row | [#13](https://github.com/dpastling/plethora/pull/13) |
| `getfasta -name` appends `::chrom:start-end` from bedtools 2.27, so the GC join silently matches nothing | [#14](https://github.com/dpastling/plethora/pull/14) |
| `preprocessing_1000genomes.R` does not parse: a misplaced parenthesis in `write.table` | [#15](https://github.com/dpastling/plethora/pull/15) |
| The README names three scripts that do not exist | [#16](https://github.com/dpastling/plethora/pull/16) |
| `trim_qc_report.R` reads `logs/trim_stats.txt`, which nothing in the repository writes | [#17](https://github.com/dpastling/plethora/pull/17) |

`logs/trim_stats.txt` is the first file `trim_qc_report.R` opens, and no script
produces it. `2_trim.sh` sends cutadapt's output to the job log and nothing
turns those logs into the three-column table the R reads, so the QC report
cannot run on a fresh checkout. Here `plethora trim` writes the file itself,
where the counts already exist, and `plethora qc-report` reads it. Pass
`--no-log` to suppress it.

`as.tbl()` was deprecated in dplyr 1.0.0 (May 2020) but only made **defunct** in
1.2.0. For five years it merely warned, so anyone who ran the pipeline before
1.2.0 got a warning and correct output. It is a forward-looking break.

Two more were noticed and not submitted. `download_sample.pl` reads `READ_COUNT`
from column 25 of the sequence index, where it is column 23, but the array it
fills is never used. `install_tools.sh` downloads two archives it never extracts
and installs into directories it does not create; anyone with the tools already
installed never runs it.

## Licence

Plethora-rs is GPL-3.0-only, because it links
[trim-galore](https://crates.io/crates/trim-galore), which is GPL-3.0-only. The
upstream plethora scripts are MIT, which is compatible in this direction: MIT
code may be incorporated into a GPL work.
