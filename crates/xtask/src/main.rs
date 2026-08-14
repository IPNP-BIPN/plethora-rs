//! Repository maintenance, run as `cargo xtask <task>`.
//!
//! Two things live here that do not belong in the binary a user installs:
//! vendoring upstream's reference data into `data/`, and running both pipelines
//! side by side to diff their intermediates. The second is what turns
//! "byte-identical" from a claim into something a reader can re-run.

mod compare;
mod data;
mod manifest;

use std::path::PathBuf;

use anyhow::Result;
use clap::{Parser, Subcommand};

/// The commit the vendored copies were taken at. Upstream has not moved since
/// 2017, but pinning means a future push cannot silently change what `data/`
/// means.
const UPSTREAM_COMMIT: &str = "4f4a5c734c81baf058dc7fcbf5d6cebc6155a211";

#[derive(Debug, Parser)]
#[command(name = "xtask", about = "Repository maintenance for plethora-rs")]
struct Cli {
    #[command(subcommand)]
    task: Task,
}

#[derive(Debug, Subcommand)]
enum Task {
    /// Download upstream's reference data and vendor it into `data/`.
    FetchData {
        /// The upstream commit to take the files from.
        #[arg(long, default_value = UPSTREAM_COMMIT)]
        commit: String,
    },

    /// Verify `data/` against its manifest. This is what CI runs.
    CheckData,

    /// Run both pipelines on the same input and diff every intermediate.
    Compare {
        /// A checkout of dpastling/plethora. Cloned to a temporary directory
        /// if not given.
        #[arg(long)]
        upstream: Option<PathBuf>,
        /// Keep the working directory instead of removing it, to inspect a
        /// disagreement.
        #[arg(long)]
        keep: bool,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    // Tasks run from the workspace root whatever directory cargo was invoked
    // from, since `data/` and `crates/` are named relative to it.
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()?;

    match cli.task {
        Task::FetchData { commit } => data::fetch(&root, &commit),
        Task::CheckData => data::check(&root),
        Task::Compare { upstream, keep } => compare::run(&root, upstream.as_deref(), keep),
    }
}
