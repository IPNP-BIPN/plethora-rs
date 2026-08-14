//! `gc_correction.R`: normalise coverage for GC bias, then for ploidy.
//!
//! ```text
//! code/gc_correction.R results/HG00250_read_depth.bed data/hg38_duf_full_domains_v2.3_GC.txt
//! ```
//!
//! Conserved regions assumed to sit at diploid copy number carry the whole
//! calibration. They define the GC correction curve, and their median corrected
//! coverage sets the haploid unit that every domain is finally divided by. So
//! the copy numbers the pipeline reports are relative to those regions, and any
//! change in how they are selected moves every result.
//!
//! The fit is R's `loess`, which is why [`plethora_compat::loess`] exists.

use std::collections::HashMap;
use std::hash::BuildHasher;
use std::io::{self, Write};

use plethora_compat::loess::eval::Loess;
use plethora_compat::rmath::{format_as_r, fround};

/// The GC window the model is fitted over.
///
/// Below and above it the correction is held flat at the nearest fitted value
/// rather than extrapolated.
pub const MIN_GC: f64 = 0.2;
/// Exclusive upper bound of the fitted window.
pub const MAX_GC: f64 = 0.73;
/// Domains quieter than this are left out of the fit.
pub const MIN_COVERAGE: f64 = 5e-2;

/// One output row.
#[derive(Debug, Clone, PartialEq)]
pub struct Row {
    pub domain: String,
    pub coverage: f64,
    pub percent_gc: f64,
    pub k_gc: f64,
    pub corrected_coverage: f64,
}

/// True for the conserved regions the calibration rests on.
///
/// The R script uses two different patterns for this, and the difference is
/// deliberate enough to be worth keeping: `^((baseline)|(uc))` is anchored when
/// selecting the domains that define the GC curve, and `((baseline)|(uc))` is
/// not when selecting the ones that set the haploid unit. A domain whose name
/// merely contains "uc" therefore contributes to the ploidy normalisation but
/// not to the curve.
#[must_use]
pub fn is_conserved(domain: &str, anchored: bool) -> bool {
    if anchored {
        domain.starts_with("baseline") || domain.starts_with("uc")
    } else {
        domain.contains("baseline") || domain.contains("uc")
    }
}

/// Why a GC model could not be fitted.
#[derive(Debug, Clone, thiserror::Error, PartialEq, Eq)]
pub enum ModelError {
    /// The inner join produced nothing.
    ///
    /// The likely cause is worth naming: `bedtools getfasta -name` wrote plain
    /// domain names up to version 2.26 and `name::chrom:start-end` from 2.27
    /// on, so a GC file built with a modern bedtools has keys that cannot match
    /// the domain names in the read-depth file. `build_gc_model.sh` still
    /// passes `-name`; the fix upstream is `-nameOnly`.
    #[error(
        "no domain survived the join between the read-depth and GC files; if the GC file was \
         built with bedtools 2.27 or later, its keys carry a `::chrom:start-end` suffix that \
         cannot match the domain names (use `getfasta -nameOnly`)"
    )]
    EmptyJoin,

    /// Too few GC bins for a local quadratic.
    ///
    /// R warns "span too small. fewer data values than degrees of freedom" and
    /// carries on into an out-of-bounds read its own source marks with a
    /// FIXME. Refusing is the only defensible reading of that.
    #[error(
        "only {bins} GC bin(s) qualified, which is fewer than the {needed} a local quadratic needs"
    )]
    TooFewBins { bins: usize, needed: usize },
}

/// The fitted correction factor per GC bin.
#[derive(Debug, Clone)]
pub struct GcModel {
    /// Bins in ascending GC order, each with its correction factor.
    pub bins: Vec<(f64, f64)>,
}

impl GcModel {
    /// The factor for a GC value, held flat outside the fitted window.
    ///
    /// A value inside the window that matched no bin gets 1, which is R's
    /// `ifelse(is.na(k.gc), 1, k.gc)`: no correction rather than no result.
    #[must_use]
    pub fn factor_for(&self, percent_gc: f64) -> f64 {
        if self.bins.is_empty() {
            return 1.0;
        }
        if percent_gc < MIN_GC {
            return self.bins[0].1;
        }
        if percent_gc >= MAX_GC {
            return self.bins[self.bins.len() - 1].1;
        }
        self.bins
            .iter()
            .find(|(gc, _)| *gc == percent_gc)
            .map_or(1.0, |(_, k)| *k)
    }
}

/// Fits the GC model from the conserved domains.
///
/// The pipeline: keep the anchored conserved domains inside the GC window and
/// above the coverage floor, take logs, average per GC bin, smooth with
/// `loess`, and divide the overall mean by the smoothed value. Bins come out in
/// ascending GC order because `group_by` sorts them, and that order is what the
/// loess sees.
///
/// # Errors
/// Returns [`ModelError::TooFewBins`] when too few bins qualify for the local
/// quadratic to be determined.
pub fn fit_model(rows: &[(String, f64, f64)]) -> Result<GcModel, ModelError> {
    let mut by_bin: Vec<(f64, Vec<f64>)> = Vec::new();

    for (domain, coverage, percent_gc) in rows {
        if !is_conserved(domain, true) {
            continue;
        }
        if *percent_gc < MIN_GC || *percent_gc >= MAX_GC {
            continue;
        }
        if *coverage <= MIN_COVERAGE {
            continue;
        }
        let log_coverage = coverage.ln();
        match by_bin.iter_mut().find(|(gc, _)| *gc == *percent_gc) {
            Some((_, values)) => values.push(log_coverage),
            None => by_bin.push((*percent_gc, vec![log_coverage])),
        }
    }

    // loess with span 0.75 and degree 2 takes floor(n * 0.75) neighbours and
    // needs three of them to determine a quadratic, so eight bins is the
    // smallest input R fits without complaint.
    let needed = 8;
    if by_bin.len() < needed {
        return Err(ModelError::TooFewBins {
            bins: by_bin.len(),
            needed,
        });
    }

    // group_by orders the bins ascending, and the loess reads them in that order.
    by_bin.sort_by(|a, b| a.0.total_cmp(&b.0));

    let x: Vec<f64> = by_bin.iter().map(|(gc, _)| *gc).collect();
    let y: Vec<f64> = by_bin
        .iter()
        .map(|(_, values)| values.iter().sum::<f64>() / values.len() as f64)
        .collect();

    let overall_mean: f64 = y.iter().sum::<f64>() / y.len() as f64;
    let fitted = Loess::fit(&x, &y).fitted(&x);

    Ok(GcModel {
        bins: x
            .iter()
            .zip(&fitted)
            .map(|(gc, y_hat)| (*gc, overall_mean / y_hat))
            .collect(),
    })
}

/// Runs the correction over a read-depth table and a GC table.
///
/// The join is an inner one on domain name, keeping the read-depth order: a
/// domain missing from either file simply does not appear. That is the point
/// where a GC file built with a modern `bedtools getfasta -name` produces
/// nothing at all, because its keys carry a `::chrom:start-end` suffix; see
/// [`ModelError::EmptyJoin`].
///
/// # Errors
/// Returns an error if the join is empty or the GC model cannot be fitted.
pub fn correct<S: BuildHasher>(
    read_depth: &[(String, f64)],
    gc: &HashMap<String, f64, S>,
) -> Result<Vec<Row>, ModelError> {
    // Inner join, then round the GC to two decimals as R does. The rounding is
    // R's, not a naive one: the GC file is full of three-decimal values, which
    // land on rounding boundaries constantly.
    let joined: Vec<(String, f64, f64)> = read_depth
        .iter()
        .filter_map(|(domain, coverage)| {
            gc.get(domain)
                .map(|percent_gc| (domain.clone(), *coverage, fround(*percent_gc, 2.0)))
        })
        .collect();

    if joined.is_empty() {
        return Err(ModelError::EmptyJoin);
    }

    let model = fit_model(&joined)?;

    let corrected: Vec<Row> = joined
        .iter()
        .map(|(domain, coverage, percent_gc)| {
            let k_gc = model.factor_for(*percent_gc);
            Row {
                domain: domain.clone(),
                coverage: *coverage,
                percent_gc: *percent_gc,
                k_gc,
                // Correcting in log space, then back.
                corrected_coverage: (coverage.ln() * k_gc).exp(),
            }
        })
        .collect();

    // The haploid unit: half the median corrected coverage of the conserved
    // regions, selected with the unanchored pattern.
    let mut conserved: Vec<f64> = corrected
        .iter()
        .filter(|r| is_conserved(&r.domain, false))
        .map(|r| r.corrected_coverage)
        .collect();
    let haploid = median(&mut conserved) / 2.0;

    Ok(corrected
        .into_iter()
        .map(|mut r| {
            r.corrected_coverage /= haploid;
            r
        })
        .collect())
}

/// R's `median`: the middle value, or the mean of the two middle ones.
fn median(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    values.sort_by(f64::total_cmp);
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2]
    } else {
        // R averages the two middle values as (a + b) / 2. f64::midpoint
        // computes the same quantity by a different route and can land on a
        // different last bit, so the plain form is kept.
        #[allow(clippy::manual_midpoint)]
        {
            (values[n / 2 - 1] + values[n / 2]) / 2.0
        }
    }
}

/// Writes the table to a path, compressing when the name ends in `.gz`.
///
/// A 623,699-domain table is about 40 MB written plain and 10 MB gzipped, and
/// there is one per sample, which is why cohorts keep them compressed.
///
/// # Errors
/// Returns an error if the file cannot be created or written.
pub fn write_table_to(rows: &[Row], path: &std::path::Path) -> io::Result<()> {
    write_table(rows, crate::io::create(path)?)
}

/// Writes the table as `write.table(sep = "\t", row.names = FALSE, quote = FALSE)`.
///
/// # Errors
/// Returns an error if writing fails.
pub fn write_table<W: Write>(rows: &[Row], mut out: W) -> io::Result<()> {
    writeln!(
        out,
        "domain\tcoverage\tpercent.gc\tk.gc\tcorrected.coverage"
    )?;
    for r in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            r.domain,
            format_as_r(r.coverage),
            format_as_r(r.percent_gc),
            format_as_r(r.k_gc),
            format_as_r(r.corrected_coverage)
        )?;
    }
    out.flush()
}

/// The output name the R script derives from its input.
///
/// `gsub("_read.depth.bed", "_gc_correct.txt", f)` is a regex, and the dots
/// match any character, which is the only reason it matches the
/// `_read_depth.bed` that `make_bed.sh` actually writes. Reproduced with the
/// same latitude so a caller passing either spelling gets the same answer.
#[must_use]
pub fn output_name(input: &str) -> String {
    let pattern = "_read.depth.bed";
    let bytes = input.as_bytes();
    for start in 0..bytes.len().saturating_sub(pattern.len() - 1) {
        let matches = pattern
            .as_bytes()
            .iter()
            .enumerate()
            .all(|(i, &p)| p == b'.' || bytes[start + i] == p);
        if matches {
            return format!(
                "{}_gc_correct.txt{}",
                &input[..start],
                &input[start + pattern.len()..]
            );
        }
    }
    input.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The two patterns really do differ, and a name can satisfy one but not
    /// the other.
    #[test]
    fn the_two_conserved_patterns_differ() {
        assert!(is_conserved("baseline_1_1", true));
        assert!(is_conserved("uc001abc", true));
        assert!(
            !is_conserved("NBPF1_uc_CON1", true),
            "anchored: must start with it"
        );
        assert!(
            is_conserved("NBPF1_uc_CON1", false),
            "unanchored: contains it"
        );
        assert!(!is_conserved("NBPF1_CON1_1", false));
    }

    #[test]
    fn the_factor_is_flat_outside_the_window() {
        let model = GcModel {
            bins: vec![(0.3, 1.1), (0.4, 1.2), (0.5, 1.3)],
        };
        assert_eq!(
            model.factor_for(0.1),
            1.1,
            "below the window: the first factor"
        );
        assert_eq!(
            model.factor_for(0.9),
            1.3,
            "above the window: the last factor"
        );
        assert_eq!(model.factor_for(0.4), 1.2);
    }

    /// A GC value inside the window with no matching bin gets no correction
    /// rather than no answer.
    #[test]
    fn an_unmatched_bin_inside_the_window_gets_one() {
        let model = GcModel {
            bins: vec![(0.3, 1.1), (0.5, 1.3)],
        };
        assert_eq!(model.factor_for(0.4), 1.0);
    }

    #[test]
    fn an_empty_model_corrects_nothing() {
        let model = GcModel { bins: Vec::new() };
        assert_eq!(model.factor_for(0.4), 1.0);
    }

    /// Only anchored conserved domains inside the window and above the floor
    /// reach the fit. Padded to eight bins, the smallest input a local
    /// quadratic can be determined from.
    #[test]
    fn the_fit_selects_on_three_criteria() {
        let mut rows = vec![
            // Excluded: not conserved.
            ("NBPF1_CON1_1".to_string(), 30.0, 0.40),
            // Excluded: below the GC window.
            ("baseline_lo".to_string(), 30.0, 0.10),
            // Excluded: at the exclusive upper bound.
            ("baseline_hi".to_string(), 30.0, 0.73),
            // Excluded: coverage at the floor.
            ("baseline_quiet".to_string(), 0.05, 0.41),
        ];
        // Eight qualifying bins.
        for i in 0..8 {
            rows.push((
                format!("baseline_{i}"),
                30.0 + f64::from(i),
                0.30 + f64::from(i) * 0.01,
            ));
        }

        let model = fit_model(&rows).expect("eight bins is enough");
        assert_eq!(
            model.bins.len(),
            8,
            "only the qualifying bins reach the fit"
        );
        assert!(model.bins.iter().all(|(gc, _)| *gc >= 0.30 && *gc <= 0.37));
    }

    /// Domains sharing a GC bin are averaged in log space before the fit, so
    /// the bin count follows the distinct GC values rather than the row count.
    #[test]
    fn domains_sharing_a_bin_are_averaged() {
        let mut rows = Vec::new();
        for i in 0..8 {
            let gc = 0.30 + f64::from(i) * 0.01;
            rows.push((format!("baseline_{i}a"), 10.0, gc));
            rows.push((format!("baseline_{i}b"), 40.0, gc));
        }
        let model = fit_model(&rows).expect("eight bins");
        assert_eq!(model.bins.len(), 8, "eight bins from sixteen rows");
    }

    #[test]
    fn bins_come_out_in_ascending_gc_order() {
        // Deliberately supplied in descending order.
        let rows: Vec<(String, f64, f64)> = (0..8)
            .rev()
            .map(|i| {
                (
                    format!("baseline_{i}"),
                    30.0 + f64::from(i),
                    0.30 + f64::from(i) * 0.01,
                )
            })
            .collect();
        let model = fit_model(&rows).expect("eight bins");
        let order: Vec<f64> = model.bins.iter().map(|(gc, _)| *gc).collect();
        let mut sorted = order.clone();
        sorted.sort_by(f64::total_cmp);
        assert_eq!(order, sorted);
    }

    /// Too few bins is refused rather than fitted, because R's own code reads
    /// out of bounds there.
    #[test]
    fn too_few_bins_is_an_error() {
        let rows = vec![
            ("baseline_a".to_string(), 10.0, 0.30),
            ("baseline_b".to_string(), 20.0, 0.31),
        ];
        assert_eq!(
            fit_model(&rows).unwrap_err(),
            ModelError::TooFewBins { bins: 2, needed: 8 }
        );
    }

    /// An empty join names the modern-bedtools trap rather than failing
    /// obscurely later.
    #[test]
    fn an_empty_join_says_why() {
        let read_depth = vec![("domA".to_string(), 1.0)];
        let gc: HashMap<String, f64> = [("domA::chr1:0-12".to_string(), 0.4)].into_iter().collect();
        let err = correct(&read_depth, &gc).unwrap_err();
        assert_eq!(err, ModelError::EmptyJoin);
        assert!(
            err.to_string().contains("nameOnly"),
            "the message must name the fix"
        );
    }

    #[test]
    fn median_handles_both_parities() {
        assert_eq!(median(&mut [3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&mut [4.0, 1.0, 3.0, 2.0]), 2.5);
        assert!(median(&mut []).is_nan());
    }

    /// Conserved regions are defined as diploid, so their median corrected
    /// coverage must come out at 2.
    #[test]
    fn the_conserved_regions_normalise_to_two_copies() {
        let read_depth: Vec<(String, f64)> = (0..40)
            .map(|i| (format!("baseline_{i}"), 30.0 + f64::from(i % 5)))
            .collect();
        let gc: HashMap<String, f64> = read_depth
            .iter()
            .enumerate()
            .map(|(i, (d, _))| (d.clone(), 0.30 + (i % 20) as f64 * 0.01))
            .collect();

        let rows = correct(&read_depth, &gc).expect("enough bins");
        let mut conserved: Vec<f64> = rows.iter().map(|r| r.corrected_coverage).collect();
        let m = median(&mut conserved);
        assert!(
            (m - 2.0).abs() < 1e-9,
            "the median conserved region should sit at two copies, got {m}"
        );
    }

    /// A domain missing from either table drops out entirely.
    #[test]
    fn the_join_is_an_inner_one() {
        let read_depth = vec![("a".to_string(), 1.0), ("b".to_string(), 2.0)];
        let gc: HashMap<String, f64> = [("a".to_string(), 0.4)].into_iter().collect();
        // Only one row survives, which is too few to fit: the join is still
        // what decided that, and the error says so.
        assert_eq!(
            correct(&read_depth, &gc).unwrap_err(),
            ModelError::TooFewBins { bins: 0, needed: 8 }
        );
    }

    /// The dots in the R pattern match any character, which is the only reason
    /// it matches the file the pipeline actually writes.
    #[test]
    fn the_output_name_follows_the_r_substitution() {
        assert_eq!(
            output_name("results/HG00250_read_depth.bed"),
            "results/HG00250_gc_correct.txt"
        );
        assert_eq!(
            output_name("results/HG00250_read.depth.bed"),
            "results/HG00250_gc_correct.txt"
        );
        assert_eq!(output_name("results/other.txt"), "results/other.txt");
    }

    /// The same table, written compressed, reads back identically.
    #[test]
    fn the_table_round_trips_through_gzip() {
        use std::io::Read as _;

        let rows = vec![Row {
            domain: "baseline_1_1".into(),
            coverage: 33.2787,
            percent_gc: 0.41,
            k_gc: 0.990_282_562_851_831,
            corrected_coverage: 2.222_236_689_446_68,
        }];

        let dir = tempfile::tempdir().unwrap();
        let plain = dir.path().join("s_gc_correct.txt");
        let gz = dir.path().join("s_gc_correct.txt.gz");
        write_table_to(&rows, &plain).unwrap();
        write_table_to(&rows, &gz).unwrap();

        assert!(std::fs::metadata(&gz).unwrap().len() > 0);
        let mut compressed = String::new();
        crate::io::open(&gz)
            .unwrap()
            .read_to_string(&mut compressed)
            .unwrap();
        assert_eq!(compressed, std::fs::read_to_string(&plain).unwrap());
        assert!(compressed.contains("baseline_1_1\t33.2787\t0.41\t"));
    }

    #[test]
    fn the_table_carries_a_header_and_no_quotes() {
        let rows = vec![Row {
            domain: "d".into(),
            coverage: 28.2794,
            percent_gc: 0.34,
            k_gc: 1.0,
            corrected_coverage: 2.0,
        }];
        let mut out = Vec::new();
        write_table(&rows, &mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "domain\tcoverage\tpercent.gc\tk.gc\tcorrected.coverage\nd\t28.2794\t0.34\t1\t2\n"
        );
    }
}
