//! Running many samples, locally or through a scheduler.
//!
//! Upstream is seven LSF scripts chained by dependency:
//!
//! ```text
//! bsub < code/1000genomes/1_download.sh
//! bsub -w "done(download[*])" < code/1000genomes/2_trim.sh
//! bsub -w "done(trim[*])" < code/1000genomes/3_batch_bowtie.sh
//! ...
//! ```
//!
//! Each is a job array indexed by `LSB_JOBINDEX` into the `SAMPLES` array in
//! `config.sh`. That is a good design on a cluster and unusable anywhere else,
//! so this offers both: [`run`] does the same work locally with a bounded number
//! of samples in flight, and [`emit`] writes the job scripts for a scheduler
//! that already exists.
//!
//! One behaviour is worth keeping from the array model: a sample that fails does
//! not stop the others. The scripts are independent by construction, and a
//! cohort where one download times out should not lose the other three hundred.

use std::fmt;
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::config::Config;

/// One stage of the pipeline, in the order they run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Step {
    /// Fetch a sample's reads.
    Download,
    /// Trim and filter them.
    Trim,
    /// Align to the reference.
    Align,
    /// Turn the alignment into per-domain coverage.
    Coverage,
    /// Remove intermediates whose successor checks out.
    Clean,
    /// Correct for GC bias and normalise for ploidy.
    GcCorrect,
}

impl Step {
    /// Every step, in order.
    pub const ALL: [Self; 6] = [
        Self::Download,
        Self::Trim,
        Self::Align,
        Self::Coverage,
        Self::Clean,
        Self::GcCorrect,
    ];

    /// The name used on the command line and in job scripts.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Download => "download",
            Self::Trim => "trim",
            Self::Align => "align",
            Self::Coverage => "coverage",
            Self::Clean => "clean",
            Self::GcCorrect => "gc-correct",
        }
    }

    /// The steps from `first` to `last` inclusive.
    ///
    /// An inverted range yields nothing rather than the whole pipeline, which
    /// is what asking for `--from coverage --to trim` deserves.
    #[must_use]
    pub fn range(first: Self, last: Self) -> Vec<Self> {
        Self::ALL
            .into_iter()
            .filter(|s| *s >= first && *s <= last)
            .collect()
    }
}

impl fmt::Display for Step {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

impl std::str::FromStr for Step {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, String> {
        Self::ALL
            .into_iter()
            .find(|step| step.name() == s)
            .ok_or_else(|| {
                let names: Vec<&str> = Self::ALL.iter().map(|s| s.name()).collect();
                format!("unknown step {s}, expected one of {}", names.join(", "))
            })
    }
}

/// What happened to one sample.
#[derive(Debug, Clone)]
pub struct Outcome {
    pub sample: String,
    /// The step that failed, or the last one that ran.
    pub step: Step,
    /// The failure, if there was one.
    pub error: Option<String>,
    pub elapsed: Duration,
}

impl Outcome {
    /// Whether every step ran.
    #[must_use]
    pub const fn succeeded(&self) -> bool {
        self.error.is_none()
    }
}

/// Runs `steps` over every sample, `jobs` at a time.
///
/// `work` is called once per sample and step, in order; the first failure for a
/// sample stops that sample and leaves the rest running. Outcomes come back in
/// the sample list's order, not completion order, so a run is reproducible to
/// read.
///
/// # Panics
/// Panics if a rayon thread pool of the requested size cannot be built.
pub fn run<F>(config: &Config, steps: &[Step], work: F) -> Vec<Outcome>
where
    F: Fn(&str, Step) -> Result<(), String> + Sync,
{
    let jobs = config.options.jobs.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(jobs)
        .build()
        .expect("a thread pool of the requested size");

    pool.install(|| {
        config
            .samples
            .par_iter()
            .map(|sample| {
                let started = Instant::now();
                let mut last = steps.first().copied().unwrap_or(Step::Download);
                for step in steps {
                    last = *step;
                    if let Err(error) = work(sample, *step) {
                        return Outcome {
                            sample: sample.clone(),
                            step: *step,
                            error: Some(error),
                            elapsed: started.elapsed(),
                        };
                    }
                }
                Outcome {
                    sample: sample.clone(),
                    step: last,
                    error: None,
                    elapsed: started.elapsed(),
                }
            })
            .collect()
    })
}

/// Which scheduler to write job scripts for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scheduler {
    /// LSF, which is what upstream targets.
    Lsf,
    /// Slurm.
    Slurm,
}

impl std::fmt::Display for Scheduler {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Lsf => "lsf",
            Self::Slurm => "slurm",
        })
    }
}

impl std::str::FromStr for Scheduler {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, String> {
        match s.to_ascii_lowercase().as_str() {
            "lsf" | "bsub" => Ok(Self::Lsf),
            "slurm" | "sbatch" => Ok(Self::Slurm),
            other => Err(format!("unknown scheduler {other}, expected lsf or slurm")),
        }
    }
}

/// Writes a job array script for one step.
///
/// The array is indexed into the sample list in the configuration, so the
/// script and the configuration must stay together; that is the same coupling
/// `LSB_JOBINDEX` and `SAMPLES` already have, made explicit.
#[must_use]
pub fn emit(config: &Config, step: Step, scheduler: Scheduler, config_path: &str) -> String {
    let n = config.samples.len();
    let name = step.name();
    let threads = config.options.threads;
    let logs = config.paths.logs.display();

    match scheduler {
        Scheduler::Lsf => format!(
            "#!/usr/bin/env bash\n\
             #BSUB -J {name}[1-{n}]%{concurrency}\n\
             #BSUB -e {logs}/{name}_%J.log\n\
             #BSUB -o {logs}/{name}_%J.out\n\
             #BSUB -n {threads}\n\
             #BSUB -q normal\n\
             \n\
             set -o nounset -o pipefail -o errexit\n\
             \n\
             # The array index selects a sample from the list in {config_path}.\n\
             plethora {name} --config {config_path} --index \"$LSB_JOBINDEX\"\n",
            concurrency = config.options.jobs.max(1),
        ),
        Scheduler::Slurm => format!(
            "#!/usr/bin/env bash\n\
             #SBATCH --job-name={name}\n\
             #SBATCH --array=1-{n}%{concurrency}\n\
             #SBATCH --error={logs}/{name}_%A_%a.log\n\
             #SBATCH --output={logs}/{name}_%A_%a.out\n\
             #SBATCH --cpus-per-task={threads}\n\
             \n\
             set -o nounset -o pipefail -o errexit\n\
             \n\
             # The array index selects a sample from the list in {config_path}.\n\
             plethora {name} --config {config_path} --index \"$SLURM_ARRAY_TASK_ID\"\n",
            concurrency = config.options.jobs.max(1),
        ),
    }
}

/// The sample a job array index selects.
///
/// One-based, because both `LSB_JOBINDEX` and `SLURM_ARRAY_TASK_ID` are.
#[must_use]
pub fn sample_at(config: &Config, index: usize) -> Option<&String> {
    config.samples.get(index.checked_sub(1)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    fn config(samples: &[&str], jobs: usize) -> Config {
        let list = samples
            .iter()
            .map(|s| format!("\"{s}\""))
            .collect::<Vec<_>>()
            .join(", ");
        let text = format!(
            "samples = [{list}]\n\
             [reference]\n\
             domains = \"d.bed\"\n\
             gc_table = \"d_GC.txt\"\n\
             [options]\n\
             jobs = {jobs}\n"
        );
        toml::from_str(&text).expect("a valid config")
    }

    #[test]
    fn steps_parse_and_print_by_name() {
        assert_eq!("coverage".parse::<Step>().unwrap(), Step::Coverage);
        assert_eq!(Step::GcCorrect.to_string(), "gc-correct");
        let err = "nonsense".parse::<Step>().unwrap_err();
        assert!(
            err.contains("download"),
            "the message should list them: {err}"
        );
    }

    #[test]
    fn a_range_of_steps_is_inclusive() {
        assert_eq!(
            Step::range(Step::Trim, Step::Coverage),
            [Step::Trim, Step::Align, Step::Coverage]
        );
        assert_eq!(Step::range(Step::Trim, Step::Trim), [Step::Trim]);
    }

    /// Asking for the steps backwards yields nothing rather than everything.
    #[test]
    fn an_inverted_range_is_empty() {
        assert!(Step::range(Step::Coverage, Step::Trim).is_empty());
    }

    #[test]
    fn every_sample_runs_every_step_in_order() {
        let config = config(&["A", "B"], 1);
        let seen: Mutex<Vec<String>> = Mutex::new(Vec::new());
        let outcomes = run(
            &config,
            &Step::range(Step::Trim, Step::Align),
            |sample, step| {
                seen.lock().unwrap().push(format!("{sample}:{step}"));
                Ok(())
            },
        );

        assert_eq!(outcomes.len(), 2);
        assert!(outcomes.iter().all(Outcome::succeeded));
        let seen = seen.lock().unwrap();
        assert_eq!(seen.len(), 4);
        // Within a sample the order holds, whatever order the samples ran in.
        let a: Vec<&String> = seen.iter().filter(|s| s.starts_with("A:")).collect();
        assert_eq!(a, ["A:trim", "A:align"]);
    }

    /// The array model's useful property: one sample failing leaves the rest.
    #[test]
    fn a_failure_stops_its_sample_and_no_other() {
        let config = config(&["good", "bad", "also_good"], 2);
        let outcomes = run(&config, &Step::ALL, |sample, step| {
            if sample == "bad" && step == Step::Align {
                return Err("the aligner fell over".into());
            }
            Ok(())
        });

        assert_eq!(outcomes.len(), 3);
        let bad = outcomes.iter().find(|o| o.sample == "bad").unwrap();
        assert!(!bad.succeeded());
        assert_eq!(bad.step, Step::Align, "it stopped where it failed");
        assert!(bad.error.as_deref().unwrap().contains("fell over"));

        assert!(
            outcomes
                .iter()
                .filter(|o| o.sample != "bad")
                .all(Outcome::succeeded)
        );
    }

    /// A failing sample runs no step after the one that failed.
    #[test]
    fn a_failure_stops_the_later_steps() {
        let config = config(&["S"], 1);
        let ran: Mutex<Vec<Step>> = Mutex::new(Vec::new());
        run(&config, &Step::ALL, |_, step| {
            ran.lock().unwrap().push(step);
            if step == Step::Trim {
                return Err("nope".into());
            }
            Ok(())
        });
        assert_eq!(*ran.lock().unwrap(), [Step::Download, Step::Trim]);
    }

    /// Outcomes come back in the sample list's order, not completion order.
    #[test]
    fn outcomes_are_reported_in_list_order() {
        let config = config(&["z", "a", "m"], 3);
        let outcomes = run(&config, &[Step::Trim], |_, _| Ok(()));
        let names: Vec<&str> = outcomes.iter().map(|o| o.sample.as_str()).collect();
        assert_eq!(names, ["z", "a", "m"]);
    }

    #[test]
    fn an_empty_sample_list_runs_nothing() {
        assert!(run(&config(&[], 4), &Step::ALL, |_, _| Ok(())).is_empty());
    }

    /// The index is one-based, as both schedulers make it.
    #[test]
    fn the_array_index_is_one_based() {
        let config = config(&["A", "B", "C"], 1);
        assert_eq!(sample_at(&config, 1).unwrap(), "A");
        assert_eq!(sample_at(&config, 3).unwrap(), "C");
        assert_eq!(sample_at(&config, 0), None, "there is no job zero");
        assert_eq!(sample_at(&config, 4), None);
    }

    #[test]
    fn schedulers_parse_by_name_or_command() {
        assert_eq!("lsf".parse::<Scheduler>().unwrap(), Scheduler::Lsf);
        assert_eq!("bsub".parse::<Scheduler>().unwrap(), Scheduler::Lsf);
        assert_eq!("SLURM".parse::<Scheduler>().unwrap(), Scheduler::Slurm);
        assert!("pbs".parse::<Scheduler>().is_err());
    }

    /// The emitted array must span the sample list and read the right index
    /// variable, or every job processes the wrong sample.
    #[test]
    fn an_lsf_script_indexes_the_whole_list() {
        let config = config(&["A", "B", "C"], 2);
        let script = emit(&config, Step::Coverage, Scheduler::Lsf, "plethora.toml");
        assert!(script.contains("#BSUB -J coverage[1-3]%2"), "{script}");
        assert!(script.contains("$LSB_JOBINDEX"));
        assert!(script.contains("plethora coverage --config plethora.toml"));
    }

    #[test]
    fn a_slurm_script_indexes_the_whole_list() {
        let config = config(&["A", "B", "C"], 2);
        let script = emit(&config, Step::Trim, Scheduler::Slurm, "run.toml");
        assert!(script.contains("#SBATCH --array=1-3%2"), "{script}");
        assert!(script.contains("$SLURM_ARRAY_TASK_ID"));
        assert!(script.contains("--config run.toml"));
    }

    #[test]
    fn the_thread_count_reaches_the_scheduler() {
        let mut config = config(&["A"], 1);
        config.options.threads = 24;
        assert!(emit(&config, Step::Align, Scheduler::Lsf, "c.toml").contains("#BSUB -n 24"));
        assert!(
            emit(&config, Step::Align, Scheduler::Slurm, "c.toml").contains("--cpus-per-task=24")
        );
    }
}
