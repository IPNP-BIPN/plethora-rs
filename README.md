# plethora-rs

Native Rust reimplementation of [plethora](https://github.com/dpastling/plethora),
the Olduvai/DUF1220 copy-number pipeline, byte-identical target. Work in
progress. Not the official plethora.

The original is four shell scripts, four Perl scripts and two R scripts driving
cutadapt, bowtie2, samtools, bedtools, GNU sort, awk and R. Its own README warns
that "updates to samtools and bedtools may break plethora". This replaces the
lot with one binary.

> Astling DP, Heft IE, Jones KL, Sikela JM. "High resolution measurement of
> DUF1220 domain copy number from whole genome sequence data" (2017) BMC
> Genomics 18:614. https://doi.org/10.1186/s12864-017-3976-z

## Status

| Upstream stage | Replaced by | Verified against |
|---|---|---|
| `cutadapt` | [trim-galore](https://crates.io/crates/trim-galore) | cutadapt 5.2, exact on the no-adapter path |
| `bowtie2` | [bwa-mem4](https://github.com/IPNP-BIPN/bwa-mem4) | pending IPNP-BIPN/bwa-mem4#61 |
| `samtools sort -n` | `bam::namesort` | samtools 1.24, 1216 records |
| `bedtools bamtobed` | `bam::bamtobed` | bedtools 2.31.1, 800 records |
| `merge_pairs.pl` | `merge_pairs` | the Perl itself, 3000 lines |
| `sort -k1,1 -k2,2n` | `bed::sort` | GNU sort, 4000 lines |
| `bedtools intersect -wao` | `bed::intersect` | bedtools 2.31.1, 5185 lines |
| `bedtools merge -c -o sum` | `bed::merge` | bedtools 2.31.1, 360 rows |
| `awk` read depth | `coverage` | awk `OFMT`, 360 rows |
| `gc_from_fasta.pl` | `gc::from_fasta` | the Perl itself, 200 sequences |
| `build_gc_model.sh` | `gc::build_model` | bedtools getfasta, 60 domains |
| `gc_correction.R` | `gc::correction` | R 4.6.1, 623,699 real domains |
| `download_sample.pl` | `onekg::download` | the Perl's column indices |
| `preprocessing_1000genomes.R` | `onekg::preprocess` | the R, quota selection |
| `trim_qc_report.R` | `onekg::qc_report` | the R, minus the deletions |
| `clean_files.pl` | `clean` | the Perl itself |
| `zip.sh` | `onekg::align_report` | GNU sort and `uniq` |
| `config.sh` + `#BSUB` scripts | `config`, `batch` | LSF and Slurm arrays |
| the scripts, called by hand | the `plethora` binary | end to end, on both paths |

Still to come: `plethora align` over bwa-mem4, which waits on
[IPNP-BIPN/bwa-mem4#61](https://github.com/IPNP-BIPN/bwa-mem4/pull/61) being
merged and 4.3.2 published, and the vendored reference data.

## The interesting part

Reproducing a pipeline is not reimplementing the documented algorithms. It is
reproducing the specific implementations, including the parts nobody designed on
purpose. Some of what had to be ported:

- **RANDLIB.** `merge_pairs.pl` seeds Perl's `Math::Random` from an MD5 of every
  input line and draws a normal deviate, so every extended read's length is a
  function of RANDLIB's exact arithmetic. Two details decide it: `phrtsd` has two
  implementations and the one a plain `cpanm` compiles drops the last character
  of the phrase, so plethora hashes 31 characters of each 32-character MD5; and
  `ichr` is a C `char`, which is signed, so a byte ≥ 0x80 arrives negative.

- **R's `loess`.** `gc_correction.R` divides by `predict(loess(y ~ x))`. Those
  are not the local regression: with the default `surface = "interpolate"`,
  `loess` fits only at k-d tree vertices and interpolates between them. Measured
  on a real input, exactly one fitted value in 53 is bit-equal between
  `"interpolate"` and `"direct"`. Porting the easy one would have produced
  quietly wrong numbers, so the tree, the vertex fits, LINPACK's `dqrdc`,
  `dqrsl` and `dsvdc`, and the reference BLAS under them are all here.

- **GNU `sort`'s third key.** `sort -k1,1 -k2,2n` is not stable and breaks a tie
  on every named key by comparing the whole line. In a BED file of read
  intervals that last-resort comparison orders most of the file.

- **`awk`'s `OFMT`.** The published figures read `28.2794`, not `28.279401`,
  because awk formats numbers at `%.6g`. On top of that sits a rule easy to
  miss: an exactly integral value prints as an integer, which is how every
  uncovered domain writes `0` rather than `0.00000`.

- **bedtools' block ordering.** `bamtobed -bedpe` orders the two blocks by
  comparing the chromosome *name as a string*, not the reference id, so the
  first block is not read 1, and an unmapped mate always lands first because
  `.` sorts below every chromosome name.

Every one of these is pinned by golden vectors generated from the real tool, or
by a differential test that runs it. See `DIVERGENCES.md` for what does not
match and why.

## Using it

Two entry points. One sample at a time, which mirrors calling the upstream
scripts by hand:

```sh
plethora trim -1 fastq/S_1.fastq.gz -2 fastq/S_2.fastq.gz
# align however you like; see DIVERGENCES.md for the bowtie2 recipe that
# reproduces the paper
plethora coverage -r data/hg38_duf_full_domains_v2.3.bed -p paired \
                  -b alignments/S.bam -o results/S --gzip
plethora gc-correct results/S_read_depth.bed.gz data/hg38_duf_full_domains_v2.3_GC.txt
```

Or a whole cohort from one file, which `config.sh` never offered:

```sh
plethora init                      # writes plethora.toml
plethora run -j 8                  # every sample, locally, over rayon
plethora run --from coverage --to gc-correct   # or part of the chain
plethora emit-jobs --scheduler slurm           # or hand it to a cluster
plethora clean HG00250 --rm-fastq  # says what it would delete; --apply does it
plethora qc-report                 # how far each sample got, and what looks wrong
```

`--index` selects a single sample by position, which is what the emitted job
arrays use. Any output path ending in `.gz` is written compressed, and any input
is read compressed or not according to its own bytes rather than its name.

## Validation

Beyond the per-stage differential tests, the GC correction is checked against a
real upstream run: samples of a 394-genome cohort processed with the original
scripts, 623,699 domains each.

| sample | rows | `percent.gc` disagreeing as text | worst relative |
|---|---|---|---|
| S36742 | 623,699 | 0 | 2.0e-14 |
| S36958 | 623,699 | 0 | 8.2e-14 |
| S37111 | 623,699 | 0 | 9.0e-15 |

`percent.gc` is exact: `round()` reproduces R's on every one of 1.87 million
values. The residual is in `corrected.coverage`, and it is the BLAS inside
`loess`. S36958's 8.2e-14 is the relative error on a domain whose coverage is
1.27e-5, where an absolute difference of 1e-19 reads large; the absolute
agreement is the same everywhere.

For scale, running the upstream R script itself on a different machine disagrees
with that same output on 103,907 of 623,699 rows. R does not reproduce itself
across machines any more closely than this reproduces R.

The cohort's own copy of `hg38_duf_full_domains_v2.3_GC.txt` is byte-identical
to the one vendored here, so the comparison is against the same reference data,
and `data/` is confirmed by a second source.

`cargo xtask compare` re-runs that claim rather than quoting it. It clones
upstream, builds one corpus, hands it to `code/make_bed.sh` and
`code/gc_correction.R` with the real samtools, bedtools, GNU sort, awk, Perl and
R underneath, hands the same corpus to `plethora`, and diffs what each left on
disk:

```
file                      bytes  result
.bed                     346200  identical
_sorted.bed              217200  identical
_coverage.bed              4862  identical
_read_depth.bed            3290  identical
_gc_correct.txt            5703  equal to 3.4e-15 relative (52 rows differ in
                                 their last digits)
```

Everything from the BAM to the read depth is byte-identical. The last file is
the one that goes through `loess`, and it lands on the same residual. CI runs
this on every push.

`data/` holds upstream's reference files compressed, with `data/MANIFEST.toml`
recording the SHA-256 of each plaintext and the commit it was taken at.
`cargo xtask check-data` decompresses every one and re-checks it, so a
compressed copy stays traceable to the file it came from.

## Memory

Nothing holds the alignment. `coverage` reads the BAM, filters, name-sorts and
pairs as one stream, so its footprint is the sort's run buffer rather than the
file. Measured on the same input with the same 200,000-record run size:

```
collecting   390 MB
streaming    108 MB
```

That is the difference between running and not. A record costs about 181 bytes
once its QNAME and reference name are counted, so a 30x whole genome, upwards of
a billion records, would have needed 137 GB held at once. `run -j N` multiplies
the run buffer by N, which is the number to size a node against.

### Intermediates

`make_bed.sh` writes six files and deletes two of them on the way out. Those two
never exist here: each is a pipe into the stage that reads it, on its own
thread, so the two stages overlap instead of taking turns and nothing lands on
disk. It is the same resolution the shell reaches with
`awk ... | bedtools merge -i -`, and the four files upstream keeps are still
written, byte for byte.

On a million pairs, the two builds run alternately five times each:

```
                peak disk    wall    CPU
through files      201 MB   2.16 s  2.73 s
through pipes       89 MB   1.95 s  2.93 s
```

Less than half the disk, a tenth off the wall clock, and seven percent more CPU,
which is the trade a thread buys: the stages overlap rather than taking turns.

Both ends of the pipe are buffered, which is not an optimisation but the
difference between working and not: a pipe write is a syscall, these stages emit
a line at a time, and unbuffered the million-line intermediate cost a million
syscalls and ran three times slower than the file it replaced.

## Speed

The sort runs across every core. Measured on this machine, sixteen of them,
with the orders checked identical:

| lines | sequential | parallel | |
|---|---|---|---|
| 1,000,000 | 967 ms | 144 ms | 6.7x |
| 10,000,000 | 13.73 s | 1.99 s | 6.9x |

That is the stage worth parallelising, because it is the one that grows: a
whole-genome sample sorts a hundred million intervals, where the rest of the
chain is linear passes over the same file.

Three other things were tried and are not here, because measuring them
alternately against the unchanged binary showed nothing:

- **mimalloc and jemalloc.** No measurable difference, single-sample or under
  `-j 8`, where allocator contention would show if it existed: 35.1 s of CPU
  against 35.3 s. A first measurement suggested 3x, and was a cold-cache
  artefact of comparing a first run against warm ones.
- **Parallel BGZF decoding** of the BAM. Costs about 8% more CPU for no wall
  time, even on a 79 MB alignment of random sequence that deflate cannot help.
- **rapidgzip on ordinary gzip.** A single deflate stream has no block
  boundaries to split on, so the parallel decoder guesses and comes out slower:
  46 ms against flate2's 36 ms on the same 58 MB. It is used only where the
  input is BGZF, which the reader detects from the header rather than the name.

### Compression

`--gzip` writes BGZF rather than a single deflate stream. Every gzip reader
still takes it, `bgzip -t` accepts it, and its 64 KB blocks compress and
decompress on every core. Measured on a 58 MB intermediate:

| | write | size | read back |
|---|---|---|---|
| gzip, one stream | 210 ms | 13.8 MB | 46 ms |
| BGZF | 51 ms | 11.6 MB | 6.5 ms |

End to end on a million pairs that turns `--gzip` from a two-to-three times
penalty into no wall-clock cost at all, for a third of the disk:

```
plain   3.24 s    89 MB
gzip    2.90 s    28 MB
```

The CPU is higher, about 4.5 s against 3.5 s, because compressing is real work.
It is spread across cores, so the wall clock does not move.

### Sample-level

The axis that pays is samples, not stages, and it is now the default rather
than something to remember. Eight samples through `coverage` and `gc-correct`,
straight out of `plethora init`:

```
-j 1     23.05 s
default   4.40 s     5.2x
```

The default is the core count capped at eight, because concurrency costs memory
rather than being free: each sample in flight holds one name-sort run buffer of
roughly 900 MB. Raise it against memory, not against cores. `-j` overrides it,
and the emitted job arrays carry it through to the scheduler.

## Layout

```
crates/
  plethora-compat/   bit-exact ports of Perl, R, GNU coreutils and samtools behaviour
  plethora-core/     the pipeline stages
  plethora/          the binary
  xtask/             vendoring the reference data, and diffing against upstream
data/                upstream's reference files, compressed, with a manifest
```

`plethora-compat` is the fragile half, and it is separate on purpose: every
function in it exists because some undocumented detail of a third-party tool
moves a digit in the final table, and each is pinned by vectors generated from
that tool.

## Building and testing

```
cargo build --release
cargo test
```

The differential tests skip when their reference tool is missing, so the suite
runs without the bioinformatics stack. To run them all you need samtools,
bedtools, GNU coreutils `sort` (as `gsort` on macOS), cutadapt, R with dplyr,
and Perl with `Math::Random`:

```
env -u PERL5LIB cpanm --notest -l .oracle/perl5 Math::Random
```

## Upstream

Four verified problems in the original have been reported:
[#13](https://github.com/dpastling/plethora/pull/13),
[#14](https://github.com/dpastling/plethora/pull/14),
[#15](https://github.com/dpastling/plethora/pull/15),
[#16](https://github.com/dpastling/plethora/pull/16). None affects the core path
a typical single-sample run takes, which is why the pipeline has worked for
years.

## Licence

GPL-3.0-only, because it links `trim-galore`, which is GPL-3.0-only. The
upstream plethora scripts are MIT, which is compatible in this direction.
