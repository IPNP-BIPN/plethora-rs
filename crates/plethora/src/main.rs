//! The `plethora` binary.
//!
//! One subcommand per upstream script, keeping their names and their argument
//! letters where they had them, so a `make_bed.sh` invocation translates
//! directly:
//!
//! ```text
//! code/make_bed.sh -r data/domains.bed -p paired -b alignments/test.bam -o results/test
//! plethora coverage -r data/domains.bed -p paired -b alignments/test.bam -o results/test
//! ```
//!
//! `run` drives the whole pipeline from a `plethora.toml`, and `emit-jobs`
//! writes the array scripts for a cluster that already has a scheduler.

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand};

use plethora_core::batch::{self, Scheduler, Step};
use plethora_core::config::Config;
use plethora_core::coverage::{self, Pairing};
use plethora_core::gc;
use plethora_core::io as pio;
use plethora_core::onekg;
use plethora_core::trim;

/// Olduvai/DUF1220 copy number from whole genome sequence data.
#[derive(Debug, Parser)]
#[command(name = "plethora", version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Trim and filter a pair of FASTQ files.
    Trim {
        /// First mate.
        #[arg(short = '1')]
        read1: PathBuf,
        /// Second mate.
        #[arg(short = '2')]
        read2: PathBuf,
        /// How adapters are chosen: detect, none, or a preset name.
        #[arg(long, default_value = "detect")]
        adapters: String,
        /// Append the read counts here, for `qc-report` to read later.
        #[arg(long, default_value = "logs/trim_stats.txt")]
        log: PathBuf,
        /// Do not write the trimming log.
        #[arg(long)]
        no_log: bool,
    },

    /// Coverage per domain, from an alignment. `make_bed.sh`.
    Coverage {
        /// The BED of domains to measure. Upstream's -r.
        #[arg(short = 'r')]
        reference: PathBuf,
        /// paired or single. Upstream's -p.
        #[arg(short = 'p', default_value = "paired")]
        pairing: String,
        /// The alignment. Upstream's -b.
        #[arg(short = 'b')]
        bam: PathBuf,
        /// Output prefix, not a file. Upstream's -o.
        #[arg(short = 'o')]
        output: PathBuf,
        /// Where the sorts spill.
        #[arg(long, default_value = ".")]
        temp: PathBuf,
        /// Write the outputs gzipped.
        #[arg(long)]
        gzip: bool,
    },

    /// Correct coverage for GC bias and normalise for ploidy.
    GcCorrect {
        /// A `_read_depth.bed` from `coverage`.
        read_depth: PathBuf,
        /// The per-domain GC table.
        gc_table: PathBuf,
        /// Where to write. Defaults to upstream's `_gc_correct.txt` name.
        #[arg(short = 'o')]
        output: Option<PathBuf>,
        /// Write the output gzipped.
        #[arg(long)]
        gzip: bool,
    },

    /// Percent GC per sequence in a FASTA. `gc_from_fasta.pl`.
    GcFromFasta { fasta: PathBuf },

    /// Build a GC table from a BED and an indexed genome. `build_gc_model.sh`.
    BuildGcModel {
        /// The BED of domains.
        #[arg(short = 'b')]
        bed: PathBuf,
        /// The genome FASTA. Its `.fai` index must exist.
        #[arg(short = 'f')]
        fasta: PathBuf,
        /// Where to write the table.
        #[arg(short = 'o')]
        output: PathBuf,
    },

    /// Report which intermediates can be removed. `clean_files.pl`.
    Clean {
        /// The sample name.
        sample: String,
        /// Also remove the FASTQ files. Upstream's --rm-fastq.
        #[arg(long)]
        rm_fastq: bool,
        /// Remove the files rather than only listing them.
        #[arg(long)]
        apply: bool,
        /// The 1000 Genomes sequence index, which gives the expected read
        /// count. Without it the first stage present anchors the chain.
        #[arg(long)]
        index: Option<PathBuf>,
        #[arg(long, default_value = "plethora.toml")]
        config: PathBuf,
    },

    /// Fetch a sample's FASTQ files. `download_sample.pl`.
    Download {
        sample: String,
        /// The 1000 Genomes sequence index.
        #[arg(long)]
        index: PathBuf,
        /// Where `fastq/<sample>/` goes.
        #[arg(long, default_value = ".")]
        root: PathBuf,
    },

    /// Choose which public samples to process. `preprocessing_1000genomes.R`.
    SelectSamples {
        /// The 1000 Genomes sequence index.
        #[arg(long)]
        index: PathBuf,
        /// Samples that failed QC before, one per line.
        #[arg(long)]
        failed: Option<PathBuf>,
        /// Samples to take regardless of quota, one per line. Repeatable.
        #[arg(long)]
        of_interest: Vec<PathBuf>,
    },

    /// Count distinct aligned fragments in a sorted BED. `zip.sh`.
    AlignReport {
        /// A `_sorted.bed` from `coverage`.
        #[arg(short = 'b')]
        bed: PathBuf,
        /// The sample name to record.
        #[arg(long)]
        sample: String,
        /// Append to this report instead of writing to stdout.
        #[arg(long)]
        report: Option<PathBuf>,
    },

    /// Run the pipeline over every sample in a configuration.
    Run {
        #[arg(long, default_value = "plethora.toml")]
        config: PathBuf,
        /// First step to run.
        #[arg(long, default_value = "trim")]
        from: String,
        /// Last step to run.
        #[arg(long, default_value = "gc-correct")]
        to: String,
        /// Run only the sample at this one-based index, as a job array does.
        #[arg(long)]
        index: Option<usize>,
        /// How many samples to hold in flight, overriding `options.jobs`.
        #[arg(short = 'j', long)]
        jobs: Option<usize>,
    },

    /// Write scheduler job scripts for a configuration.
    EmitJobs {
        #[arg(long, default_value = "plethora.toml")]
        config: PathBuf,
        /// lsf or slurm.
        #[arg(long, default_value = "slurm")]
        scheduler: String,
        /// Where to write the scripts.
        #[arg(long, default_value = "jobs")]
        output: PathBuf,
    },

    /// Report how far each sample got, and what looks wrong. `trim_qc_report.R`.
    QcReport {
        #[arg(long, default_value = "plethora.toml")]
        config: PathBuf,
        /// The trimming log.
        #[arg(long, default_value = "logs/trim_stats.txt")]
        trim_stats: PathBuf,
        /// The alignment report, from `align-report --report`.
        #[arg(long, default_value = "align_report.txt")]
        align_report: PathBuf,
        /// The 1000 Genomes sequence index, which gives the expected counts.
        #[arg(long)]
        index: Option<PathBuf>,
    },

    /// Write a starting `plethora.toml`.
    Init {
        #[arg(long, default_value = "plethora.toml")]
        output: PathBuf,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Trim {
            read1,
            read2,
            adapters,
            log,
            no_log,
        } => run_trim(
            &read1,
            &read2,
            &adapters,
            (!no_log).then_some(log).as_deref(),
        ),
        Command::Coverage {
            reference,
            pairing,
            bam,
            output,
            temp,
            gzip,
        } => run_coverage(&reference, &pairing, &bam, &output, &temp, gzip),
        Command::GcCorrect {
            read_depth,
            gc_table,
            output,
            gzip,
        } => run_gc_correct(&read_depth, &gc_table, output.as_deref(), gzip),
        Command::GcFromFasta { fasta } => run_gc_from_fasta(&fasta),
        Command::BuildGcModel { bed, fasta, output } => run_build_gc_model(&bed, &fasta, &output),
        Command::Clean {
            sample,
            rm_fastq,
            apply,
            index,
            config,
        } => run_clean(&sample, rm_fastq, apply, &config, index.as_deref()),
        Command::Download {
            sample,
            index,
            root,
        } => run_download(&sample, &index, &root),
        Command::SelectSamples {
            index,
            failed,
            of_interest,
        } => run_select(&index, failed.as_deref(), &of_interest),
        Command::AlignReport {
            bed,
            sample,
            report,
        } => run_align_report(&bed, &sample, report.as_deref()),
        Command::Run {
            config,
            from,
            to,
            index,
            jobs,
        } => run_pipeline(&config, &from, &to, index, jobs),
        Command::EmitJobs {
            config,
            scheduler,
            output,
        } => run_emit(&config, &scheduler, &output),
        Command::QcReport {
            config,
            trim_stats,
            align_report,
            index,
        } => run_qc_report(&config, &trim_stats, Some(&align_report), index.as_deref()),
        Command::Init { output } => run_init(&output),
    }
}

fn run_trim(read1: &Path, read2: &Path, adapters: &str, log: Option<&Path>) -> Result<()> {
    let choice = match adapters {
        "detect" => trim::Adapters::Detect,
        "none" => trim::Adapters::None,
        name => {
            let preset = trim::preset_by_name(name)
                .with_context(|| format!("unknown adapter preset: {name}"))?;
            trim::Adapters::Preset(preset.name)
        }
    };

    let summary = trim::trim_pair(read1, read2, &choice)?;
    if let Some(message) = &summary.detection {
        eprintln!("{}", message.trim_end());
    }
    println!(
        "{} pairs in, {} out, adapter {}",
        summary.pairs_in, summary.pairs_out, summary.adapter
    );
    warn_if_nothing_survived(&summary);
    println!("{}", summary.output_r1.display());
    println!("{}", summary.output_r2.display());

    // `qc-report` reads this. Upstream never writes it, so its own QC script
    // opens a file nothing produces; here the counts are already in hand.
    if let Some(log) = log {
        let discarded = (summary.pairs_in - summary.pairs_out) as u64;
        // One entry for the pair, keyed by the first mate: see append_trim_stats
        // for why a row per mate doubles every count the report checks.
        let entries = [(logged_path(read1), summary.pairs_in as u64, discarded)];
        onekg::qc_report::append_trim_stats(log, &entries)
            .with_context(|| format!("appending to {}", log.display()))?;
    }
    Ok(())
}

/// Puts words to [`trim::TrimSummary::survival`]. The judgement is in the core,
/// where it is tested; only the wording is here.
fn warn_if_nothing_survived(summary: &trim::TrimSummary) {
    match summary.survival() {
        trim::Survival::NoInput => eprintln!("warning: the input held no read pairs"),
        trim::Survival::None => eprintln!(
            "warning: no pair survived trimming. Both mates must stay at least {} bases \
             after quality trimming at Q{}, so a run whose second mate is entirely low \
             quality loses every pair. Check the per-mate quality before assuming the \
             trimmer is at fault.",
            trim::MIN_LENGTH,
            trim::QUALITY_CUTOFF
        ),
        trim::Survival::Low => eprintln!(
            "warning: only {:.1}% of pairs survived trimming",
            summary.survival_rate() * 100.0
        ),
        trim::Survival::Ordinary => {}
    }
}

/// The path as the trimming log records it, `fastq/<sample>/<file>`.
///
/// The sample is read back out of that prefix, so a path given from somewhere
/// else keeps only its last two components.
fn logged_path(path: &Path) -> String {
    let mut parts: Vec<String> = path
        .components()
        .rev()
        .take(2)
        .map(|c| c.as_os_str().to_string_lossy().into_owned())
        .collect();
    parts.reverse();
    format!("fastq/{}", parts.join("/"))
}

fn run_coverage(
    reference: &Path,
    pairing: &str,
    bam: &Path,
    output: &Path,
    temp: &Path,
    gzip: bool,
) -> Result<()> {
    let pairing: Pairing = pairing.parse().map_err(|e: String| anyhow::anyhow!(e))?;
    let compress = if gzip {
        pio::Compress::Yes
    } else {
        pio::Compress::No
    };

    let outputs = coverage::make_bed(bam, reference, pairing, output, temp, compress)
        .map_err(|e| anyhow::anyhow!("{e}"))?;
    println!("{}", outputs.read_depth.display());
    Ok(())
}

fn run_gc_correct(
    read_depth: &Path,
    gc_table: &Path,
    output: Option<&Path>,
    gzip: bool,
) -> Result<()> {
    let depth = read_pairs(read_depth)?;
    let gc: HashMap<String, f64> = read_pairs(gc_table)?.into_iter().collect();

    let rows = gc::correction::correct(&depth, &gc)?;

    // Upstream derives the name from the input, and so does this by default.
    let destination = output.map_or_else(
        || PathBuf::from(gc::correction::output_name(&read_depth.to_string_lossy())),
        Path::to_path_buf,
    );
    let destination = if gzip {
        pio::with_gzip(&destination)
    } else {
        destination
    };

    gc::correction::write_table_to(&rows, &destination)?;
    println!("{}", destination.display());
    Ok(())
}

/// Reads a two-column table of name and number, decompressing if needed.
fn read_pairs(path: &Path) -> Result<Vec<(String, f64)>> {
    use std::io::BufRead as _;
    let mut out = Vec::new();
    for line in pio::open(path)
        .with_context(|| format!("reading {}", path.display()))?
        .lines()
    {
        let line = line?;
        if let Some((name, value)) = line.split_once('\t')
            && let Ok(number) = value.trim().parse::<f64>()
        {
            out.push((name.to_string(), number));
        }
    }
    Ok(out)
}

fn run_gc_from_fasta(fasta: &Path) -> Result<()> {
    let rows = gc::from_fasta::gc_from_fasta(pio::open(fasta)?)?;
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    for row in &rows {
        writeln!(out, "{}", row.to_line())?;
    }
    Ok(())
}

fn run_build_gc_model(bed: &Path, fasta: &Path, output: &Path) -> Result<()> {
    let mut genome = gc::build_model::Genome::open(fasta)
        .with_context(|| format!("opening {}", fasta.display()))?;
    let rows = gc::build_model::build(pio::open(bed)?, &mut genome)?;

    let mut out = pio::create(output)?;
    for row in &rows {
        writeln!(out, "{}", row.to_line())?;
    }
    out.flush()?;
    println!("{} domains -> {}", rows.len(), output.display());
    Ok(())
}

fn run_clean(
    sample: &str,
    rm_fastq: bool,
    apply: bool,
    config: &Path,
    index: Option<&Path>,
) -> Result<()> {
    let config = Config::load(config)?;
    let paths = onekg::clean::Paths {
        fastq: &config.paths.fastq,
        alignments: &config.paths.alignments,
        results: &config.paths.results,
    };

    let gathered = onekg::clean::gather(sample, &paths)?;
    if gathered.counts.present() == 0 {
        println!("nothing found on disk for sample {sample}");
        return Ok(());
    }

    let pairing = if config.options.paired {
        onekg::clean::Pairing::Paired
    } else {
        onekg::clean::Pairing::Single
    };
    if let Some(named) = gathered.named_pairing
        && named != pairing
    {
        eprintln!(
            "warning: the FASTQ names look {named:?} but the configuration says {pairing:?}; \
             going with the configuration"
        );
    }

    // The sequence index is what upstream calls the manifest: without it the
    // first stage present anchors the chain instead.
    let expected = match index {
        Some(path) => {
            let records = onekg::sample_index::read(pio::open(path)?)?;
            let wanted = onekg::sample_index::for_sample(&records, sample);
            if wanted.is_empty() {
                bail!("sample {sample} is not in {}", path.display());
            }
            // The index lists both mates, so its read counts are per file.
            Some(wanted.iter().filter_map(|r| r.read_count).sum::<u64>() / 2)
        }
        None => None,
    };

    let plan = onekg::clean::plan_with_reasons(&gathered.counts, expected, pairing, rm_fastq)
        .map_err(|e| anyhow::anyhow!("{e}"))?;

    print_counts(&gathered.counts);
    if plan.is_empty() {
        println!("nothing can be removed yet for sample {sample}");
        return Ok(());
    }

    for step in &plan {
        let files = onekg::clean::files_for(step.removal, sample, &paths, &gathered);
        let present: Vec<_> = files.into_iter().filter(|f| f.exists()).collect();
        if present.is_empty() {
            continue;
        }
        // Upstream's exact line, which is what `qc-report` greps back out of
        // the clean-up logs. Only on --apply, since it asserts the removal
        // happened.
        if apply {
            println!("{}", onekg::clean::log_line(sample, *step));
        }
        for file in present {
            if apply {
                std::fs::remove_file(&file)
                    .with_context(|| format!("removing {}", file.display()))?;
                println!("  removed {}", file.display());
            } else {
                println!(
                    "would remove {} (the {} matched)",
                    file.display(),
                    step.verified
                );
            }
        }
    }
    if !apply {
        println!("\nnothing was removed; pass --apply to do it");
    }
    Ok(())
}

/// `trim_qc_report.R`, minus the deletions it does on the way.
fn run_qc_report(
    config: &Path,
    trim_stats: &Path,
    align_report: Option<&Path>,
    index: Option<&Path>,
) -> Result<()> {
    let config = Config::load(config)?;

    let stats = onekg::qc_report::read_trim_stats(pio::open(trim_stats)?)
        .with_context(|| format!("reading {}", trim_stats.display()))?;
    let summaries = onekg::qc_report::summarise(&stats);

    // Expected counts come from the sequence index, halved because it lists
    // both mates. Without one, no sample can reach the trimmed stage, which is
    // upstream's behaviour too.
    let mut expected: HashMap<String, (u64, usize)> = HashMap::new();
    if let Some(path) = index {
        let records = onekg::sample_index::read(pio::open(path)?)?;
        for record in &records {
            let entry = expected.entry(record.sample_name.clone()).or_insert((0, 0));
            entry.0 += record.read_count.unwrap_or(0);
            entry.1 += 1;
        }
        for value in expected.values_mut() {
            value.0 /= 2;
            value.1 /= 2;
        }
    }

    let mut aligned: HashMap<String, usize> = HashMap::new();
    if let Some(path) = align_report
        && path.exists()
    {
        for entry in onekg::align_report::read_report(path)? {
            aligned.insert(entry.sample, entry.aligned_fragments);
        }
    }

    // `grep "correct number of reads" logs/clean_*.out`, without the shell.
    let mut bams = HashSet::new();
    let mut beds = HashSet::new();
    if let Ok(entries) = std::fs::read_dir(&config.paths.logs) {
        let mut logs: Vec<PathBuf> = entries
            .filter_map(Result::ok)
            .map(|e| e.path())
            .filter(|p| {
                p.file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(|n| n.starts_with("clean_"))
            })
            .collect();
        logs.sort();
        for log in logs {
            let (b, d) = onekg::qc_report::read_clean_report(pio::open(&log)?)?;
            bams.extend(b);
            beds.extend(d);
        }
    }

    let reports = onekg::qc_report::build(
        &config.samples,
        &summaries,
        &expected,
        &aligned,
        &bams,
        &beds,
    );

    println!(
        "{:<5} {:<16} {:<8} {:>14} {:>14} {:>10} {:>12}",
        "index", "sample", "stage", "total reads", "after trim", "filtered", "aligned"
    );
    let mut problems = Vec::new();
    for report in &reports {
        let index = report
            .index
            .map_or_else(|| "-".to_string(), |i| i.to_string());
        let aligned = report
            .percent_aligned()
            .map_or_else(|| "-".to_string(), |p| format!("{:.1}%", p * 100.0));
        println!(
            "{index:<5} {:<16} {:<8} {:>14} {:>14} {:>9.1}% {aligned:>12}",
            report.sample,
            format!("{:?}", report.stage),
            report.total_reads,
            report.remaining_reads,
            report.percent_filtered * 100.0,
        );
        if report.has_quality_problem() {
            problems.push(report);
        }
    }

    if !problems.is_empty() {
        println!(
            "\n{} sample(s) with fewer than 100M reads left or more than 10% filtered:",
            problems.len()
        );
        for report in problems {
            println!("  {}", report.sample);
        }
    }
    println!("\nNo files were removed. Upstream's script deletes FASTQ files as a side");
    println!("effect of being run; `plethora clean --apply` is where that lives here.");
    Ok(())
}

fn print_counts(counts: &plethora_core::onekg::clean::Counts) {
    let show = |name: &str, value: Option<u64>| match value {
        Some(n) => println!("  {name:<12} {n}"),
        None => println!("  {name:<12} absent"),
    };
    show("fastq", counts.fastq);
    show("bam", counts.bam);
    show("sorted bam", counts.sorted_bam);
    show("bed", counts.bed);
}

fn run_download(sample: &str, index: &Path, root: &Path) -> Result<()> {
    let records = onekg::sample_index::read(pio::open(index)?)?;
    let wanted = onekg::sample_index::for_sample(&records, sample);
    if wanted.is_empty() {
        bail!("sample {sample} is not in {}", index.display());
    }

    let outcomes =
        onekg::download::fetch_sample(&onekg::download::HttpFetcher, root, sample, &wanted)?;
    for (file, outcome) in &outcomes {
        println!("{file}: {outcome:?}");
    }
    Ok(())
}

fn run_select(index: &Path, failed: Option<&Path>, of_interest: &[PathBuf]) -> Result<()> {
    let records = onekg::sample_index::read(pio::open(index)?)?;

    let mut lists = onekg::preprocess::Lists::default();
    if let Some(path) = failed {
        lists.failed = read_names(path)?;
    }
    for path in of_interest {
        lists.sudmant.extend(read_names(path)?);
    }

    for sample in onekg::preprocess::select(&records, &lists) {
        println!("{sample}");
    }
    Ok(())
}

/// Reads a list of sample names, one per line.
fn read_names(path: &Path) -> Result<HashSet<String>> {
    use std::io::BufRead as _;
    Ok(pio::open(path)
        .with_context(|| format!("reading {}", path.display()))?
        .lines()
        .map_while(Result::ok)
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect())
}

fn run_align_report(bed: &Path, sample: &str, report: Option<&Path>) -> Result<()> {
    let entry = onekg::align_report::Entry {
        sample: sample.to_string(),
        aligned_fragments: onekg::align_report::count_for_sample(bed)?,
    };
    match report {
        Some(path) => onekg::align_report::append_line(path, &entry)?,
        None => println!("{}", entry.to_line()),
    }
    Ok(())
}

fn run_pipeline(
    config_path: &Path,
    from: &str,
    to: &str,
    index: Option<usize>,
    jobs: Option<usize>,
) -> Result<()> {
    let mut config = Config::load(config_path)?;
    if let Some(jobs) = jobs {
        if jobs == 0 {
            bail!("-j must be at least one");
        }
        config.options.jobs = jobs;
    }
    config.check_reference()?;
    config.create_directories()?;

    let first: Step = from.parse().map_err(|e: String| anyhow::anyhow!(e))?;
    let last: Step = to.parse().map_err(|e: String| anyhow::anyhow!(e))?;
    let steps = Step::range(first, last);
    if steps.is_empty() {
        bail!("no steps between {from} and {to}; they may be the wrong way round");
    }

    // A job array runs one sample per task, selected by its index.
    if let Some(index) = index {
        let sample = batch::sample_at(&config, index)
            .with_context(|| format!("no sample at index {index}"))?
            .clone();
        config.samples = vec![sample];
    }

    let outcomes = batch::run(&config, &steps, |sample, step| {
        step_for_sample(&config, sample, step).map_err(|e| format!("{e:#}"))
    });

    let failed = outcomes.iter().filter(|o| !o.succeeded()).count();
    for outcome in &outcomes {
        match &outcome.error {
            Some(error) => eprintln!("{}: {} failed: {error}", outcome.sample, outcome.step),
            None => println!("{}: through {}", outcome.sample, outcome.step),
        }
    }
    if failed > 0 {
        bail!("{failed} of {} samples failed", outcomes.len());
    }
    Ok(())
}

/// Runs one step for one sample.
///
/// Every stage that exists is driven here. Alignment is the exception: it needs
/// bwa-mem4's command implementations, which are only reachable by spawning the
/// binary until IPNP-BIPN/bwa-mem4#61 lands a library target. Until then the
/// step says so and points at the BAM it expects, rather than pretending.
fn step_for_sample(config: &Config, sample: &str, step: Step) -> Result<()> {
    match step {
        Step::Download => {
            let index = config
                .reference
                .sample_index
                .as_ref()
                .context("downloading needs reference.sample_index in the configuration")?;
            let records = onekg::sample_index::read(pio::open(index)?)?;
            let wanted = onekg::sample_index::for_sample(&records, sample);
            if wanted.is_empty() {
                bail!("sample {sample} is not in {}", index.display());
            }
            onekg::download::fetch_sample(
                &onekg::download::HttpFetcher,
                Path::new("."),
                sample,
                &wanted,
            )?;
            Ok(())
        }

        Step::Trim => {
            let dir = config.fastq_dir(sample);
            let (read1, read2) = mate_pair(&dir)
                .with_context(|| format!("looking for reads in {}", dir.display()))?;
            let summary = trim::trim_pair(&read1, &read2, &trim::Adapters::Detect)?;
            warn_if_nothing_survived(&summary);
            // The counts go to the log `qc-report` reads. Every sample appends
            // to the same file, which is why the reader keeps the last entry
            // per file and kind rather than the first.
            let discarded = (summary.pairs_in - summary.pairs_out) as u64;
            let entries = [(logged_path(&read1), summary.pairs_in as u64, discarded)];
            onekg::qc_report::append_trim_stats(
                &config.paths.logs.join("trim_stats.txt"),
                &entries,
            )?;
            Ok(())
        }

        Step::Align => {
            let bam = config.alignment(sample);
            if bam.exists() {
                return Ok(());
            }
            bail!(
                "no aligner is wired in yet, and {} does not exist. Produce it \
                 with bowtie2 for parity with the paper, or with bwa-mem4, then \
                 rerun from the coverage step",
                bam.display()
            )
        }

        Step::Coverage => {
            let pairing = if config.options.paired {
                Pairing::Paired
            } else {
                Pairing::Single
            };
            coverage::make_bed(
                &config.alignment(sample),
                &config.reference.domains,
                pairing,
                &config.result_prefix(sample),
                &config.paths.temp,
                config.compress(),
            )
            .map_err(|e| anyhow::anyhow!("{e}"))?;
            Ok(())
        }

        Step::Clean => {
            // Reporting only: removing files is an explicit choice, not a step
            // a pipeline run takes on its own.
            Ok(())
        }

        Step::GcCorrect => {
            let prefix = config.result_prefix(sample);
            let read_depth = config.compress().apply(&PathBuf::from(format!(
                "{}_read_depth.bed",
                prefix.display()
            )));
            let depth = read_pairs(&read_depth)?;
            let gc: HashMap<String, f64> = read_pairs(&config.reference.gc_table)?
                .into_iter()
                .collect();
            let rows = gc::correction::correct(&depth, &gc)?;
            let output = config.compress().apply(&PathBuf::from(format!(
                "{}_gc_correct.txt",
                prefix.display()
            )));
            gc::correction::write_table_to(&rows, &output)?;
            Ok(())
        }
    }
}

/// Finds a sample's two mate files in its FASTQ directory.
///
/// The names follow upstream's convention, `*_1.fastq*` and `*_2.fastq*`, which
/// is what its batch scripts glob for.
fn mate_pair(dir: &Path) -> Result<(PathBuf, PathBuf)> {
    let mut first: Option<PathBuf> = None;
    let mut second: Option<PathBuf> = None;

    for entry in std::fs::read_dir(dir)? {
        let path = entry?.path();
        let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        // A file already trimmed is not an input.
        if name.contains("_filtered.") {
            continue;
        }
        if name.contains("_1.fastq") && first.is_none() {
            first = Some(path);
        } else if name.contains("_2.fastq") && second.is_none() {
            second = Some(path);
        }
    }

    match (first, second) {
        (Some(a), Some(b)) => Ok((a, b)),
        _ => bail!(
            "expected a *_1.fastq* and a *_2.fastq* in {}",
            dir.display()
        ),
    }
}

fn run_emit(config_path: &Path, scheduler: &str, output: &Path) -> Result<()> {
    let config = Config::load(config_path)?;
    let scheduler: Scheduler = scheduler.parse().map_err(|e: String| anyhow::anyhow!(e))?;

    // An empty list yields the array range 1-0, which both schedulers reject at
    // submission. Better to say so now than to hand over scripts that fail.
    if config.samples.is_empty() {
        bail!(
            "{} lists no samples, so the job arrays would be empty",
            config_path.display()
        );
    }
    std::fs::create_dir_all(output)?;

    for (i, step) in Step::ALL.iter().enumerate() {
        let script = batch::emit(&config, *step, scheduler, &config_path.to_string_lossy());
        // The scheduler is in the name so that emitting both does not leave one
        // silently overwritten by the other.
        let path = output.join(format!("{}_{}.{scheduler}.sh", i + 1, step.name()));
        std::fs::write(&path, script)?;
        println!("{}", path.display());
    }
    Ok(())
}

fn run_init(output: &Path) -> Result<()> {
    if output.exists() {
        bail!("{} already exists", output.display());
    }
    let template = r#"# plethora.toml, which replaces config.sh.

# Samples to process, in order. A job array indexes into this list.
samples = []

[reference]
domains = "data/hg38_duf_full_domains_v2.3.bed"
gc_table = "data/hg38_duf_full_domains_v2.3_GC.txt"
# index = "genomes/hg38"
# fasta = "genomes/hg38.fa"
# sample_index = "data/1000Genomes_samples.txt"

[paths]
fastq = "fastq"
alignments = "alignments"
results = "results"
logs = "logs"
temp = "."

[options]
paired = true
# Upstream's cutadapt -q 10 and --minimum-length 80.
quality_cutoff = 10
min_length = 80
compress = false
threads = 12
jobs = 1
aligner_args = []
"#;
    std::fs::write(output, template)?;
    println!("{}", output.display());
    Ok(())
}
