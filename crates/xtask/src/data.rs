//! Vendoring the reference data, and checking the vendored copies.

use std::io::{BufRead, Read, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::manifest::{Digested, Entry, Manifest, Source, digest};

/// Which upstream files are vendored, and which get compressed.
///
/// The threshold is not a rule, it is a judgement: the three large files are
/// 46 MB between them and compress to 7 MB, while the four lists are under a
/// kilobyte each and are worth keeping greppable in the tree.
const FILES: &[(&str, bool, &str)] = &[
    (
        "hg38_duf_full_domains_v2.3.bed",
        true,
        "The Olduvai/DUF1220 domain coordinates on hg38, which every run measures against",
    ),
    (
        "hg38_duf_full_domains_v2.3_GC.txt",
        true,
        "Percent GC per domain, the input to the loess correction",
    ),
    (
        "1000Genomes_samples.txt",
        true,
        "The 1000 Genomes sequence index, as of the paper",
    ),
    (
        "failed_samples.txt",
        false,
        "Samples excluded from selection because they failed QC",
    ),
    (
        "sudmant_samples.txt",
        false,
        "Samples from Sudmant et al., taken regardless of quota",
    ),
    (
        "CLM_DNA_sample_names.txt",
        false,
        "Colombian samples of interest, taken regardless of quota",
    ),
    (
        "rd-irys-samples.txt",
        false,
        "Irys samples of interest, taken regardless of quota",
    ),
];

const REPOSITORY: &str = "https://github.com/dpastling/plethora";

/// Downloads each file, checks it, and writes the vendored copy.
pub fn fetch(root: &Path, commit: &str) -> Result<()> {
    let data = root.join("data");
    std::fs::create_dir_all(&data)?;

    let mut entries = Vec::new();
    for &(name, compress, description) in FILES {
        let url = format!("{REPOSITORY}/raw/{commit}/data/{name}");
        println!("fetching {name}");
        let body = download(&url).with_context(|| format!("downloading {url}"))?;

        let got = digest(&body[..])?;
        let vendored = if compress {
            let path = data.join(format!("{name}.gz"));
            write_bgzf(&path, &body)?;
            // Read it straight back, so a vendored file that cannot be
            // decompressed is caught here rather than by whoever clones next.
            let round_trip = digest(open(&path)?)?;
            if round_trip.sha256 != got.sha256 {
                bail!("{name}: the compressed copy does not decompress to the original");
            }
            format!("{name}.gz")
        } else {
            std::fs::write(data.join(name), &body)?;
            name.to_string()
        };

        println!(
            "  {} bytes, {} lines -> data/{vendored}",
            got.bytes, got.lines
        );
        entries.push(Entry {
            name: name.to_string(),
            vendored,
            sha256: got.sha256,
            bytes: got.bytes,
            lines: got.lines,
            description: description.to_string(),
        });
    }

    let manifest = Manifest {
        source: Source {
            repository: REPOSITORY.to_string(),
            commit: commit.to_string(),
        },
        files: entries,
    };
    manifest.store(root)?;
    println!("wrote {}", crate::manifest::MANIFEST);
    Ok(())
}

/// Re-reads every vendored file and compares it against the manifest.
pub fn check(root: &Path) -> Result<()> {
    let manifest = Manifest::load(root)?;
    let mut failures = 0;

    for entry in &manifest.files {
        let path = root.join("data").join(&entry.vendored);
        if !path.exists() {
            println!("MISSING {}", entry.vendored);
            failures += 1;
            continue;
        }
        let got: Digested = match open(&path).map_err(anyhow::Error::from).and_then(digest) {
            Ok(got) => got,
            Err(e) => {
                println!("UNREADABLE {}: {e}", entry.vendored);
                failures += 1;
                continue;
            }
        };
        match entry.verify(&got) {
            Ok(()) => println!("ok {:<40} {} lines", entry.vendored, got.lines),
            Err(e) => {
                println!("MISMATCH {e}");
                failures += 1;
            }
        }
    }

    println!(
        "\n{} of {} files match {} at {}",
        manifest.files.len() - failures,
        manifest.files.len(),
        manifest.source.repository,
        &manifest.source.commit[..manifest.source.commit.len().min(12)]
    );
    if failures > 0 {
        bail!("{failures} vendored files do not match the manifest");
    }
    Ok(())
}

/// Opens a vendored file, compressed or not, the way plethora itself does.
fn open(path: &Path) -> std::io::Result<Box<dyn Read>> {
    let mut file = std::io::BufReader::new(std::fs::File::open(path)?);
    // Sniff rather than trust the name, which is what `plethora_core::io` does.
    let gzipped = matches!(file.fill_buf()?, [0x1f, 0x8b, ..]);
    if gzipped {
        Ok(Box::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

/// Writes BGZF rather than plain gzip. Every gzip reader takes it, including
/// the one in `plethora_core::io`, which is why nothing downstream had to
/// change; the reason to prefer it is that its 64 KB block framing leaves
/// indexed and parallel reads open later, which plain gzip forecloses. Nothing
/// uses that yet.
///
/// At the maximum level, because these files are written once and then live in
/// everybody's clone forever. htslib's `bgzip -l 9` still beats this by about a
/// fifth, since it deflates through libdeflate; matching it would mean linking
/// C, which is the thing this project is for not doing.
fn write_bgzf(path: &Path, body: &[u8]) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut writer = noodles_bgzf::io::writer::Builder::default()
        .set_compression_level(noodles_bgzf::io::writer::CompressionLevel::BEST)
        .build_from_writer(file);
    writer.write_all(body)?;
    writer.finish()?;
    Ok(())
}

fn download(url: &str) -> Result<Vec<u8>> {
    let mut response = ureq::get(url).call()?;
    if response.status() != 200 {
        bail!("HTTP {}", response.status());
    }
    let mut body = Vec::new();
    response
        .body_mut()
        .as_reader()
        .take(256 << 20)
        .read_to_end(&mut body)?;
    Ok(body)
}
