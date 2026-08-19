//! Cell-barcode whitelisting and error correction
//! (Python source: `umi_tools/whitelist_methods.py`).
//!
//! Given the barcodes observed in a single-cell library and how often each
//! occurred, decide which correspond to real cells. Two "knee" estimators are
//! available, and either can be followed by error correction that folds
//! near-miss barcodes into their whitelisted neighbours.
//!
//! # Plots are not reproduced; their data is
//!
//! Python's `--plot-prefix` writes three matplotlib PNGs plus a
//! `_cell_thresholds.tsv` sidecar. The PNGs are not library semantics and no
//! Rust plotting stack reproduces matplotlib's output bytes, so this module
//! returns every quantity the plots draw — the cumulative curve, the
//! distance-to-line curve, the density curve, the candidate thresholds and
//! which was selected — inside [`KneeResult`], and
//! [`write_cell_thresholds`] writes the TSV byte for byte. See
//! `port_reports/umi_tools/05_design.md` deviation D5.

use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use std::path::Path;

use crate::bktree::BkTree;
use crate::error::{PortError, Result};
use crate::fuzzy::{matches_within_errors, matches_within_substitutions};
use crate::numeric;
use crate::utilities;

/// A cell barcode, as raw bytes.
pub type Barcode = Vec<u8>;

/// Barcode to observed count, remembering the order barcodes were first seen
/// (Python source: the `collections.Counter` the `whitelist` command builds).
///
/// The order is load-bearing. `Counter.most_common()` sorts by count with a
/// **stable** sort, so equal-count barcodes come out in the counter's
/// insertion order — and that order decides which of two mutually-erroneous
/// barcodes `error_detect_above_threshold` discards. A plain `BTreeMap` would
/// order ties by barcode instead and silently pick the other one.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CellCounts {
    order: Vec<Barcode>,
    counts: BTreeMap<Barcode, u64>,
}

impl CellCounts {
    /// Creates an empty counter.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds one to `barcode`'s count, recording it if new.
    pub fn observe(&mut self, barcode: &[u8]) {
        if let Some(count) = self.counts.get_mut(barcode) {
            *count += 1;
        } else {
            self.order.push(barcode.to_vec());
            self.counts.insert(barcode.to_vec(), 1);
        }
    }

    /// Sets `barcode`'s count outright, recording it if new.
    ///
    /// `whitelist --method=umis` replaces each read count with the number of
    /// distinct UMIs after the scan, which is what this is for.
    pub fn set(&mut self, barcode: &[u8], count: u64) {
        if !self.counts.contains_key(barcode) {
            self.order.push(barcode.to_vec());
        }
        self.counts.insert(barcode.to_vec(), count);
    }

    /// `barcode`'s count, or `None` if it was never seen.
    #[must_use]
    pub fn get(&self, barcode: &[u8]) -> Option<u64> {
        self.counts.get(barcode).copied()
    }

    /// How many distinct barcodes were seen.
    #[must_use]
    pub fn len(&self) -> usize {
        self.counts.len()
    }

    /// Whether no barcodes were seen.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    /// The barcodes in first-seen order.
    #[must_use]
    pub fn barcodes(&self) -> &[Barcode] {
        &self.order
    }

    /// `(barcode, count)` pairs in first-seen order.
    pub fn iter(&self) -> impl Iterator<Item = (&Barcode, u64)> + '_ {
        self.order
            .iter()
            .map(move |b| (b, self.counts.get(b).copied().unwrap_or(0)))
    }
}

impl FromIterator<(Barcode, u64)> for CellCounts {
    fn from_iter<I: IntoIterator<Item = (Barcode, u64)>>(iter: I) -> Self {
        let mut out = Self::new();
        for (barcode, count) in iter {
            out.set(&barcode, count);
        }
        out
    }
}

/// Whitelisted barcode to the erroneous barcodes that correct to it
/// (Python source: `true_to_false_map`).
pub type TrueToFalseMap = BTreeMap<Barcode, BTreeSet<Barcode>>;

/// Erroneous barcode to its whitelisted parent, or `None` when it is within
/// threshold of more than one and so cannot be corrected
/// (Python source: `false_to_true_map`).
///
/// The `None` is load-bearing: `getUserDefinedBarcodes` writes it for
/// ambiguous barcodes and `ExtractFilterAndUpdate.filterUMIBarcode` branches
/// on it, distinguishing "not correctable because ambiguous" from "not in the
/// map at all".
pub type FalseToTrueMap = BTreeMap<Barcode, Option<Barcode>>;

/// How many candidate thresholds the density estimator evaluates
/// (Python source: `xx_values = 10000`).
pub const DENSITY_XX_VALUES: usize = 10_000;

/// The kernel bandwidth the density estimator uses
/// (Python source: `gaussian_kde(log_counts, bw_method=0.1)` — the comment
/// there reads "guassian density with hardcoded bw").
pub const DENSITY_BANDWIDTH: f64 = 0.1;

/// Counts below this fraction of the most abundant barcode are excluded from
/// the density estimate (Python source: `threshold = 0.001 *
/// cell_barcode_counts.most_common(1)[0][1]`).
pub const DENSITY_MIN_COUNT_FRACTION: f64 = 0.001;

/// Iteration cap on the distance estimator's refinement loop
/// (Python source: `max_iterations = 100`).
pub const KNEE_MAX_ITERATIONS: usize = 100;

/// What a knee estimator concluded, plus everything its plots would have drawn.
#[derive(Debug, Clone, PartialEq)]
pub struct KneeResult {
    /// The selected barcodes, or `None` when no threshold was accepted.
    ///
    /// Python returns `None` from the density estimator when no local minimum
    /// passes its heuristics; `whitelist` then either warns or exits depending
    /// on `--allow-threshold-error`.
    pub barcodes: Option<BTreeSet<Barcode>>,
    /// The same barcodes in **selection order** — most abundant first.
    ///
    /// The distance estimator returns a Python *list* in `most_common` order,
    /// not a set, and [`error_detect_above_threshold`] sorts that list by
    /// ascending count with a stable sort — so equal-count barcodes keep this
    /// order, and that decides which of two mutually-erroneous barcodes is
    /// discarded. Sorting them by barcode instead silently picks the other
    /// one. See `port_reports/umi_tools/05_design.md` § Determinism decisions.
    pub ordered_barcodes: Vec<Barcode>,
    /// Barcodes ordered by descending count, with their counts.
    pub ranked: Vec<(Barcode, u64)>,
    /// Cumulative counts over `ranked`, the curve the knee plot draws.
    pub cumulative: Vec<i64>,
    /// Distance from each point to the first-to-last chord, for the distance
    /// estimator; empty for the density estimator.
    pub dist_to_line: Vec<f64>,
    /// Index of the chosen knee point in `ranked`, for the distance estimator.
    pub knee_index: Option<usize>,
    /// How many barcodes each candidate density threshold would accept, in the
    /// order Python evaluates them (highest count index first).
    ///
    /// This is Python's `local_mins_counts`, and it is exactly what the
    /// `_cell_thresholds.tsv` sidecar reports.
    pub local_mins_counts: Vec<usize>,
    /// The candidate that was selected, as an index into `local_mins_counts`.
    pub selected_local_min: Option<usize>,
}

impl KneeResult {
    /// The number of selected barcodes, or 0 when none were.
    #[must_use]
    pub fn len(&self) -> usize {
        self.barcodes.as_ref().map_or(0, BTreeSet::len)
    }

    /// Whether no barcodes were selected.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Orders barcodes by descending count (Python source: `Counter.most_common`).
///
/// A **stable** sort over the counter's first-seen order, exactly as
/// `Counter.most_common` is — so equal-count barcodes keep the order they were
/// observed in. That tie order is observable: it decides which barcodes fall
/// inside the knee cut, and which of two mutually-erroneous barcodes
/// [`error_detect_above_threshold`] discards.
#[must_use]
pub fn most_common(counts: &CellCounts) -> Vec<(Barcode, u64)> {
    let mut ranked: Vec<(Barcode, u64)> = counts.iter().map(|(k, v)| (k.clone(), v)).collect();
    // Stable and descending: `sort_by_key` with a reversed key keeps ties in
    // first-seen order, which is what `Counter.most_common` does.
    ranked.sort_by_key(|(_, count)| std::cmp::Reverse(*count));
    ranked
}

/// Distance from each cumulative-count point to the chord joining the first
/// and last, and the index of the farthest
/// (Python source: the nested `getKneeDistance` in `getKneeEstimateDistance`).
///
/// The construction: treat `(rank, cumulative_count)` as points, draw the line
/// from the first point to the last, and find the point standing furthest off
/// it. That point is the "elbow" of the cumulative curve.
///
/// Returns `None` for fewer than two points, where the chord is degenerate and
/// Python would divide by a zero norm.
#[must_use]
pub fn knee_distance(values: &[i64]) -> Option<(Vec<f64>, usize)> {
    let n = values.len();
    if n < 2 {
        return None;
    }
    let first = (0.0_f64, values[0] as f64);
    let last = ((n - 1) as f64, values[n - 1] as f64);
    let line = (last.0 - first.0, last.1 - first.1);
    let norm = (line.0 * line.0 + line.1 * line.1).sqrt();
    if norm == 0.0 {
        return None;
    }
    let unit = (line.0 / norm, line.1 / norm);

    let dist: Vec<f64> = (0..n)
        .map(|i| {
            let from_first = (i as f64 - first.0, values[i] as f64 - first.1);
            // Project onto the chord, then take the norm of what is left.
            let scalar = from_first.0 * unit.0 + from_first.1 * unit.1;
            let parallel = (scalar * unit.0, scalar * unit.1);
            let perpendicular = (from_first.0 - parallel.0, from_first.1 - parallel.1);
            (perpendicular.0 * perpendicular.0 + perpendicular.1 * perpendicular.1).sqrt()
        })
        .collect();

    let idx = numeric::argmax(&dist)?;
    Some((dist, idx))
}

/// Selects cell barcodes by the maximum-distance knee method
/// (Python source: `whitelist_methods.getKneeEstimateDistance`).
///
/// The estimator is sensitive to how many barcodes it is shown and
/// overestimates when given too many, so Python re-runs it on the leading
/// `3 * knee_index` barcodes until the index stops moving or 100 iterations
/// pass. Both the refinement and the cap are reproduced.
///
/// `cell_number`, when given, overrides the knee entirely: the threshold
/// becomes the count of the `cell_number`-th most abundant barcode and every
/// barcode with a **strictly greater** count is selected. `Some(0)` is treated
/// as absent, because Python tests the value for truthiness and `0` is falsy —
/// deviation-free, but easy to get wrong.
///
/// # Errors
///
/// - [`PortError::Value`] with Python's message `"Something's gone wrong
///   here!!"` when the first pass puts the knee at index 0, which means the
///   cumulative curve has no elbow.
/// - [`PortError::Value`] when `cell_number` indexes past the barcodes
///   available, where Python raises `IndexError`.
/// - [`PortError::Value`] for fewer than two barcodes, where the chord is
///   degenerate.
///
/// # Examples
///
/// ```
/// use umi_tools_core::whitelist::{knee_estimate_distance, CellCounts};
///
/// let counts: CellCounts = [
///     (b"AAAA".to_vec(), 1000), (b"AAAT".to_vec(), 800),
///     (b"AATT".to_vec(), 600),  (b"ATTT".to_vec(), 400),
///     (b"TTTT".to_vec(), 200),  (b"CCCC".to_vec(), 10),
///     (b"CCCG".to_vec(), 8),    (b"CCGG".to_vec(), 6),
///     (b"CGGG".to_vec(), 4),    (b"GGGG".to_vec(), 2),
/// ].into_iter().collect();
///
/// let result = knee_estimate_distance(&counts, None)?;
/// assert_eq!(result.len(), 4);
/// # Ok::<(), umi_tools_core::PortError>(())
/// ```
pub fn knee_estimate_distance(
    counts: &CellCounts,
    cell_number: Option<usize>,
) -> Result<KneeResult> {
    let cell_number = cell_number.filter(|n| *n != 0);
    let ranked = most_common(counts);
    let raw_counts: Vec<i64> = ranked.iter().map(|(_, c)| *c as i64).collect();
    let cumulative = numeric::cumsum_i64(&raw_counts);

    let Some((mut dist_to_line, mut idx)) = knee_distance(&cumulative) else {
        return Err(PortError::value(
            "cannot estimate a knee from fewer than two barcodes",
        ));
    };
    if idx == 0 {
        return Err(PortError::value("Something's gone wrong here!!"));
    }

    let mut previous = 0usize;
    let mut iterations = 0usize;
    while idx != previous {
        previous = idx;
        iterations += 1;
        if iterations > KNEE_MAX_ITERATIONS {
            break;
        }
        let take = (idx * 3).min(cumulative.len());
        let Some((d, i)) = knee_distance(&cumulative[..take]) else {
            break;
        };
        dist_to_line = d;
        idx = i;
    }

    let knee_barcodes: BTreeSet<Barcode> = ranked
        .iter()
        .take(idx + 1)
        .map(|(b, _)| b.clone())
        .collect();

    let barcodes: BTreeSet<Barcode> = match cell_number {
        Some(n) => {
            let threshold = raw_counts.get(n).copied().ok_or_else(|| {
                PortError::value(format!(
                    "--set-cell-number={n} exceeds the {} observed cell barcodes",
                    raw_counts.len()
                ))
            })?;
            counts
                .iter()
                .filter(|(_, c)| *c as i64 > threshold)
                .map(|(b, _)| b.clone())
                .collect()
        }
        None => knee_barcodes,
    };
    // Most-abundant-first, which is the order Python's list carries.
    let ordered_barcodes: Vec<Barcode> = ranked
        .iter()
        .filter(|(b, _)| barcodes.contains(b))
        .map(|(b, _)| b.clone())
        .collect();

    Ok(KneeResult {
        barcodes: Some(barcodes),
        ordered_barcodes,
        ranked,
        cumulative,
        dist_to_line,
        knee_index: Some(idx),
        local_mins_counts: Vec::new(),
        selected_local_min: None,
    })
}

/// A one-dimensional Gaussian kernel density estimate with a fixed bandwidth
/// factor (Python dependency: `scipy.stats.gaussian_kde`).
///
/// Reproduces `SciPy`'s scalar-`bw_method` convention exactly: the kernel
/// covariance is the **unbiased** sample variance (`ddof = 1`) times the
/// square of the bandwidth factor, and the estimate at `x` is
///
/// ```text
/// sum_i exp(-0.5 * (x - x_i)^2 / covariance) / (n * sqrt(2*pi*covariance))
/// ```
///
/// Confirmed against `SciPy` 1.17.1 to within 1.1e-16 over the whole evaluation
/// grid, not derived from its documentation.
#[derive(Debug, Clone)]
pub struct GaussianKde {
    dataset: Vec<f64>,
    covariance: f64,
    norm: f64,
}

impl GaussianKde {
    /// Fits an estimate to `dataset` with bandwidth factor `bw_method`.
    ///
    /// # Errors
    ///
    /// Returns [`PortError::Value`] when the dataset has fewer than two points
    /// or zero variance, where `SciPy` raises
    /// `LinAlgError: singular matrix` because the covariance cannot be
    /// inverted.
    pub fn new(dataset: &[f64], bw_method: f64) -> Result<Self> {
        if dataset.len() < 2 {
            return Err(PortError::value(
                "gaussian_kde requires at least two data points",
            ));
        }
        let covariance = numeric::variance(dataset, 1) * bw_method * bw_method;
        if covariance.is_nan() || covariance <= 0.0 {
            return Err(PortError::value(
                "gaussian_kde: dataset has zero variance (singular covariance)",
            ));
        }
        let norm = (dataset.len() as f64) * (2.0 * std::f64::consts::PI * covariance).sqrt();
        Ok(Self {
            dataset: dataset.to_vec(),
            covariance,
            norm,
        })
    }

    /// The kernel covariance, `SciPy`'s `kde.covariance[0][0]`.
    #[must_use]
    pub fn covariance(&self) -> f64 {
        self.covariance
    }

    /// Evaluates the density at one point.
    #[must_use]
    pub fn eval(&self, x: f64) -> f64 {
        let acc: f64 = self
            .dataset
            .iter()
            .map(|xi| {
                let d = x - xi;
                (-0.5 * d * d / self.covariance).exp()
            })
            .sum();
        acc / self.norm
    }

    /// Evaluates the density at each point of `xs`.
    ///
    /// Summation is sequential and single-threaded, deliberately: a parallel
    /// reduction would reorder the sum and push this comparison from Tier 2 to
    /// Tier 3 for no useful gain at these sizes.
    #[must_use]
    pub fn eval_many(&self, xs: &[f64]) -> Vec<f64> {
        xs.iter().map(|&x| self.eval(x)).collect()
    }
}

/// Selects cell barcodes by the density method
/// (Python source: `whitelist_methods.getKneeEstimateDensity`).
///
/// Fits a kernel density estimate to the base-10 log of the counts, finds the
/// local minima of that density, and picks one as the cut between real cells
/// and background. Which minimum is picked depends on what the caller knows:
///
/// - With `cell_number`, no minimum is used at all: the threshold is the count
///   of the `cell_number`-th most abundant barcode.
/// - With `expect_cells`, the highest-count minimum accepting between
///   `0.1 * expect_cells` and `expect_cells` barcodes is chosen.
/// - With neither, Python applies a heuristic its own comment describes as
///   "very unlikely to be the best way to achieve this": the minimum must lie
///   in the upper 80% of the evaluation grid, and must either be more than 0.5
///   log-units below the largest count or below half of it.
///
/// When no minimum qualifies, [`KneeResult::barcodes`] is `None` — Python's
/// `final_barcodes = None` — and `whitelist` reports that no threshold was
/// accepted.
///
/// `Some(0)` for either optional argument is treated as absent, mirroring
/// Python's truthiness test.
///
/// # Errors
///
/// - [`PortError::Value`] when fewer than two distinct counts survive the
///   0.1%-of-maximum filter, where `SciPy`'s `gaussian_kde` fails on a singular
///   covariance.
/// - [`PortError::Value`] when `cell_number` indexes past the barcodes
///   available.
pub fn knee_estimate_density(
    counts: &CellCounts,
    expect_cells: Option<u64>,
    cell_number: Option<usize>,
) -> Result<KneeResult> {
    let expect_cells = expect_cells.filter(|n| *n != 0);
    let cell_number = cell_number.filter(|n| *n != 0);

    let ranked = most_common(counts);
    let raw_counts: Vec<i64> = ranked.iter().map(|(_, c)| *c as i64).collect();
    let cumulative = numeric::cumsum_i64(&raw_counts);

    let max_count = raw_counts.first().copied().unwrap_or(0) as f64;
    let min_count = DENSITY_MIN_COUNT_FRACTION * max_count;
    let thresholded: Vec<f64> = raw_counts
        .iter()
        .map(|&c| c as f64)
        .filter(|&c| c > min_count)
        .collect();
    let log_counts = numeric::log10(&thresholded);

    let (Some(log_min), Some(log_max)) = (numeric::min(&log_counts), numeric::max(&log_counts))
    else {
        return Err(PortError::value(
            "no cell barcodes survive the density estimator's abundance filter",
        ));
    };

    let mut local_mins_counts: Vec<usize> = Vec::new();
    let mut selected_local_min: Option<usize> = None;
    let threshold: Option<f64>;

    if let Some(n) = cell_number {
        threshold = Some(raw_counts.get(n).copied().ok_or_else(|| {
            PortError::value(format!(
                "--set-cell-number={n} exceeds the {} observed cell barcodes",
                raw_counts.len()
            ))
        })? as f64);
    } else {
        let kde = GaussianKde::new(&log_counts, DENSITY_BANDWIDTH)?;
        let xx = numeric::linspace(log_min, log_max, DENSITY_XX_VALUES);
        let density = kde.eval_many(&xx);
        let local_mins = numeric::argrelmin(&density);

        // Python walks the minima from the highest count downwards
        // (`local_mins[::-1]`) and keeps the first that qualifies, while still
        // recording the accept-count of every one for the sidecar TSV.
        let mut chosen: Option<f64> = None;
        for (position, &candidate) in local_mins.iter().rev().enumerate() {
            let cut = 10.0_f64.powf(xx[candidate]);
            let passing = counts.iter().filter(|(_, c)| *c as f64 > cut).count();
            local_mins_counts.push(passing);

            if selected_local_min.is_some() {
                continue;
            }
            let qualifies = match expect_cells {
                Some(expected) => {
                    let expected = expected as f64;
                    passing as f64 > expected * 0.1 && passing as f64 <= expected
                }
                None => {
                    // The heuristic, verbatim from the Python source.
                    candidate as f64 >= 0.2 * DENSITY_XX_VALUES as f64
                        && (log_max - xx[candidate] > 0.5 || xx[candidate] < log_max / 2.0)
                }
            };
            if qualifies {
                selected_local_min = Some(position);
                chosen = Some(cut);
            }
        }
        threshold = chosen;
    }

    let barcodes = threshold.map(|cut| {
        counts
            .iter()
            .filter(|(_, c)| *c as f64 > cut)
            .map(|(b, _)| b.clone())
            .collect::<BTreeSet<Barcode>>()
    });
    // Python's density estimator returns a `set`, whose order is a hash
    // artefact; most-abundant-first is deterministic and matches the order the
    // distance estimator's list already carries.
    let ordered_barcodes: Vec<Barcode> = match barcodes.as_ref() {
        Some(selected) => ranked
            .iter()
            .filter(|(b, _)| selected.contains(b))
            .map(|(b, _)| b.clone())
            .collect(),
        None => Vec::new(),
    };

    Ok(KneeResult {
        barcodes,
        ordered_barcodes,
        ranked,
        cumulative,
        dist_to_line: Vec::new(),
        knee_index: None,
        local_mins_counts,
        selected_local_min,
    })
}

/// Which knee estimator to use (Python source: the `--knee-method` choice).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KneeMethod {
    /// Maximum distance from the cumulative curve's chord.
    Distance,
    /// Local minimum of a kernel density estimate.
    Density,
}

impl std::str::FromStr for KneeMethod {
    type Err = PortError;

    /// # Errors
    ///
    /// Returns [`PortError::Value`] with Python's message
    /// `"knee_method must be 'distance' or 'density'"`.
    fn from_str(s: &str) -> Result<Self> {
        match s {
            "distance" => Ok(KneeMethod::Distance),
            "density" => Ok(KneeMethod::Density),
            _ => Err(PortError::value(
                "knee_method must be 'distance' or 'density'",
            )),
        }
    }
}

/// Writes the `_cell_thresholds.tsv` sidecar
/// (Python source: the `U.openFile("%s_cell_thresholds.tsv" ...)` blocks in
/// both estimators).
///
/// The two estimators write different files: the density one lists every
/// candidate threshold with `Selected`/`Rejected`, the distance one writes a
/// single `count` column holding the knee index. Both are reproduced byte for
/// byte, because the test suite compares them.
///
/// Neither writes anything when `cell_number` was supplied, matching the
/// `if not cell_number:` guard in both.
///
/// # Errors
///
/// Returns [`PortError::Io`] if the file cannot be written.
pub fn write_cell_thresholds(
    result: &KneeResult,
    method: KneeMethod,
    plotfile_prefix: &Path,
    cell_number: Option<usize>,
) -> Result<()> {
    if cell_number.as_ref().is_some_and(|n| *n != 0) {
        return Ok(());
    }
    let path = plotfile_prefix.with_file_name(format!(
        "{}_cell_thresholds.tsv",
        plotfile_prefix
            .file_name()
            .map_or_else(String::new, |n| n.to_string_lossy().into_owned())
    ));
    let mut out = utilities::create_file(&path, false)?;

    match method {
        KneeMethod::Distance => {
            writeln!(out, "count").map_err(|e| PortError::io(&path, e))?;
            writeln!(out, "{}", result.knee_index.unwrap_or(0))
                .map_err(|e| PortError::io(&path, e))?;
        }
        KneeMethod::Density => {
            writeln!(out, "count\taction").map_err(|e| PortError::io(&path, e))?;
            let selected_count = result
                .selected_local_min
                .map(|i| result.local_mins_counts[i]);
            for &count in &result.local_mins_counts {
                // Python compares each candidate's accept-count against
                // `len(final_barcodes)`, so a rejected candidate that happens
                // to accept the same number is also labelled "Selected".
                let action = if Some(count) == selected_count {
                    "Selected"
                } else {
                    "Rejected"
                };
                writeln!(out, "{count}\t{action}").map_err(|e| PortError::io(&path, e))?;
            }
        }
    }
    out.finish()
}

/// Maps whitelisted barcodes to the observed barcodes that correct to them
/// (Python source: `whitelist_methods.getErrorCorrectMapping`).
///
/// An observed barcode is corrected only when it is within `threshold` of
/// **exactly one** whitelisted barcode; ambiguous ones are dropped, as are
/// barcodes with no whitelisted neighbour. Barcodes already on the whitelist
/// are skipped.
///
/// # Deviation
///
/// Python normalises the whitelist with `set([str(x) for x in whitelist])`,
/// which turns a `bytes` barcode into the literal text `"b'ACGT'"`. That is
/// defensive coercion of an unknown input type, not intended semantics — every
/// in-package caller passes `str`, where it is a no-op. This function takes
/// bytes directly. See `05_design.md` deviation D7.
#[must_use]
pub fn get_error_correct_mapping(
    cell_barcodes: &[Barcode],
    whitelist: &BTreeSet<Barcode>,
    threshold: u32,
) -> TrueToFalseMap {
    utilities::info("building bktree");
    let tree = BkTree::from_items(whitelist.iter().cloned());
    utilities::info("done building bktree");

    let mut true_to_false: TrueToFalseMap = TrueToFalseMap::new();
    for barcode in cell_barcodes {
        if whitelist.contains(barcode) {
            continue;
        }
        // `d > 0` excludes the barcode itself, which the whitelist check has
        // already ruled out but which Python filters again.
        let candidates: Vec<Barcode> = tree
            .find(barcode, threshold)
            .into_iter()
            .filter(|(d, _)| *d > 0)
            .map(|(_, item)| item)
            .collect();

        if candidates.len() == 1 {
            true_to_false
                .entry(candidates[0].clone())
                .or_default()
                .insert(barcode.clone());
        }
    }
    true_to_false
}

/// Runs the chosen knee estimator and, optionally, error correction
/// (Python source: `whitelist_methods.getCellWhitelist`).
///
/// Error correction runs only when a whitelist was selected **and**
/// `error_correct_threshold` is greater than zero, matching Python's
/// `if cell_whitelist and error_correct_threshold > 0:`.
///
/// # Errors
///
/// Propagates the chosen estimator's errors.
pub fn get_cell_whitelist(
    counts: &CellCounts,
    knee_method: KneeMethod,
    expect_cells: Option<u64>,
    cell_number: Option<usize>,
    error_correct_threshold: u32,
) -> Result<(KneeResult, Option<TrueToFalseMap>)> {
    let result = match knee_method {
        KneeMethod::Distance => knee_estimate_distance(counts, cell_number)?,
        KneeMethod::Density => knee_estimate_density(counts, expect_cells, cell_number)?,
    };
    utilities::info("Finished - whitelist determination");

    let mut true_to_false = None;
    if let Some(whitelist) = result.barcodes.as_ref() {
        if !whitelist.is_empty() && error_correct_threshold > 0 {
            utilities::info("Starting - finding putative error cell barcodes");
            let observed: Vec<Barcode> = counts.barcodes().to_vec();
            true_to_false = Some(get_error_correct_mapping(
                &observed,
                whitelist,
                error_correct_threshold,
            ));
            utilities::info("Finished - finding putative error cell barcodes");
        }
    }
    Ok((result, true_to_false))
}

/// How a user-supplied whitelist file should be turned into a correction map
/// (Python source: the mutually exclusive `getErrorCorrection` and
/// `deriveErrorCorrection` flags of `getUserDefinedBarcodes`).
///
/// Python takes them as two independent booleans even though only three of the
/// four combinations are meaningful, and silently prefers `derive` when both
/// are set. An enum makes the choice total.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorCorrectionMode {
    /// Read the whitelist only; no correction map.
    None,
    /// Read a comma-separated second column as the barcodes correcting to the
    /// first.
    FromFile,
    /// Enumerate every barcode within `threshold` substitutions of a
    /// whitelisted one.
    Derive,
}

/// The substitutions the derived correction map considers at each position
/// (Python source: `base2errors` in `getUserDefinedBarcodes`).
///
/// Only `A`, `T`, `C` and `G` are keys, so a whitelist containing an `N`
/// makes Python raise `KeyError`; this port returns
/// [`PortError::Format`] instead, naming the offending base.
const BASE_ERRORS: [(u8, [u8; 4]); 4] = [
    (b'A', *b"TCGN"),
    (b'T', *b"ACGN"),
    (b'C', *b"TAGN"),
    (b'G', *b"TCAN"),
];

fn substitutions_for(base: u8) -> Option<[u8; 4]> {
    BASE_ERRORS
        .iter()
        .find(|(b, _)| *b == base)
        .map(|(_, s)| *s)
}

fn read_barcode_column(path: &Path, column: usize) -> Result<Vec<Vec<Vec<u8>>>> {
    use std::io::BufRead as _;
    let reader = utilities::open_file(path)?;
    let mut rows = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| PortError::io(path, e))?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<Vec<u8>> = line
            .trim_end_matches(['\n', '\r'])
            .trim()
            .split('\t')
            .map(|f| f.as_bytes().to_vec())
            .collect();
        if fields.len() <= column && column > 0 {
            return Err(PortError::format(format!(
                "{}: expected at least {} tab-separated columns",
                path.display(),
                column + 1
            )));
        }
        rows.push(fields);
    }
    Ok(rows)
}

/// Reads one or two whitelist files and optionally builds a correction map
/// (Python source: `whitelist_methods.getUserDefinedBarcodes`).
///
/// With `tsv2`, the whitelist is the Cartesian product of the two files'
/// first columns, concatenated — the layout paired-end cell barcodes use.
///
/// # Errors
///
/// - [`PortError::Assertion`] when `tsv2` is combined with
///   [`ErrorCorrectionMode::FromFile`], matching Python's `assert not
///   whitelist_tsv2, "Can only extract errors from the whitelist if a single
///   whitelist is given"`.
/// - [`PortError::Format`] when [`ErrorCorrectionMode::FromFile`] meets a line
///   with no second column (Python raises `IndexError`), or when
///   [`ErrorCorrectionMode::Derive`] meets a base outside `ATCG` (Python
///   raises `KeyError`).
/// - [`PortError::Io`] if a file cannot be read.
pub fn get_user_defined_barcodes(
    tsv: &Path,
    tsv2: Option<&Path>,
    mode: ErrorCorrectionMode,
    threshold: u32,
) -> Result<(BTreeSet<Barcode>, Option<FalseToTrueMap>)> {
    let barcodes = |rows: Vec<Vec<Vec<u8>>>| -> Vec<Barcode> {
        rows.into_iter()
            .filter_map(|mut f| {
                if f.is_empty() {
                    None
                } else {
                    Some(std::mem::take(&mut f[0]))
                }
            })
            .collect()
    };

    let combined = |first: Vec<Barcode>, second: Option<Vec<Barcode>>| -> Vec<Barcode> {
        match second {
            None => first,
            Some(second) => {
                let mut out = Vec::with_capacity(first.len() * second.len());
                for a in &first {
                    for b in &second {
                        let mut joined = a.clone();
                        joined.extend_from_slice(b);
                        out.push(joined);
                    }
                }
                out
            }
        }
    };

    match mode {
        ErrorCorrectionMode::None => {
            let first = barcodes(read_barcode_column(tsv, 0)?);
            let second = tsv2
                .map(|p| read_barcode_column(p, 0).map(barcodes))
                .transpose()?;
            Ok((combined(first, second).into_iter().collect(), None))
        }

        ErrorCorrectionMode::FromFile => {
            if tsv2.is_some() {
                return Err(PortError::assertion(
                    "Can only extract errors from the whitelist if a single whitelist is given",
                ));
            }
            let rows = read_barcode_column(tsv, 1)?;
            let mut whitelist = BTreeSet::new();
            let mut map: FalseToTrueMap = FalseToTrueMap::new();
            for row in rows {
                let parent = row[0].clone();
                whitelist.insert(parent.clone());
                for error in row[1].split(|b| *b == b',') {
                    map.insert(error.to_vec(), Some(parent.clone()));
                }
            }
            Ok((whitelist, Some(map)))
        }

        ErrorCorrectionMode::Derive => {
            let first = barcodes(read_barcode_column(tsv, 0)?);
            let second = tsv2
                .map(|p| read_barcode_column(p, 0).map(barcodes))
                .transpose()?;
            let all = combined(first, second);

            let mut whitelist = BTreeSet::new();
            let mut map: FalseToTrueMap = FalseToTrueMap::new();

            for parent in &all {
                whitelist.insert(parent.clone());
                for positions in position_tuples(parent.len(), threshold as usize) {
                    let mut choices: Vec<[u8; 4]> = Vec::with_capacity(positions.len());
                    for &p in &positions {
                        choices.push(substitutions_for(parent[p]).ok_or_else(|| {
                            PortError::format(format!(
                                "whitelist barcode {} contains base '{}', \
                                 which has no defined substitutions",
                                String::from_utf8_lossy(parent),
                                parent[p] as char
                            ))
                        })?);
                    }
                    for combo in cartesian(&choices) {
                        let mut error = parent.clone();
                        for (&pos, &base) in positions.iter().zip(combo.iter()) {
                            error[pos] = base;
                        }
                        match map.get(&error) {
                            Some(existing) => {
                                if let Some(other) = existing {
                                    utilities::info(&format!(
                                        "Error barcode {} can be assigned to more than \
                                         one possible true barcode: {} or {}",
                                        String::from_utf8_lossy(&error),
                                        String::from_utf8_lossy(other),
                                        String::from_utf8_lossy(parent)
                                    ));
                                }
                                map.insert(error, None);
                            }
                            None => {
                                map.insert(error, Some(parent.clone()));
                            }
                        }
                    }
                }
            }
            Ok((whitelist, Some(map)))
        }
    }
}

/// Every `repeat`-length tuple of positions in `0..len`, with repetition
/// (Python source: `itertools.product(range(0, len(barcode)), repeat=threshold)`).
///
/// Positions may repeat, so a threshold of 2 also enumerates every
/// single-position substitution — which is why the derived map contains
/// distance-1 barcodes even at threshold 2.
fn position_tuples(len: usize, repeat: usize) -> Vec<Vec<usize>> {
    if repeat == 0 || len == 0 {
        return if repeat == 0 {
            vec![Vec::new()]
        } else {
            Vec::new()
        };
    }
    let mut out: Vec<Vec<usize>> = vec![Vec::new()];
    for _ in 0..repeat {
        let mut next = Vec::with_capacity(out.len() * len);
        for prefix in &out {
            for p in 0..len {
                let mut extended = prefix.clone();
                extended.push(p);
                next.push(extended);
            }
        }
        out = next;
    }
    out
}

/// The Cartesian product of per-position substitution choices
/// (Python source: `itertools.product(*m_bases)`).
fn cartesian(choices: &[[u8; 4]]) -> Vec<Vec<u8>> {
    let mut out: Vec<Vec<u8>> = vec![Vec::new()];
    for choice in choices {
        let mut next = Vec::with_capacity(out.len() * choice.len());
        for prefix in &out {
            for &base in choice {
                let mut extended = prefix.clone();
                extended.push(base);
                next.push(extended);
            }
        }
        out = next;
    }
    out
}

/// Whitelisted barcodes that `barcode` could be an error of
/// (Python source: `whitelist_methods.checkError`).
///
/// Returns as soon as a second match is found, because every caller only
/// distinguishes "none", "exactly one" and "more than one". Barcodes whose
/// length differs from `barcode` by more than `errors` are skipped without
/// comparison, and `barcode` never matches itself.
///
/// # Examples
///
/// ```
/// use umi_tools_core::whitelist::check_error;
/// let whitelist = vec![b"AAAT".to_vec(), b"AATT".to_vec(),
///                      b"AAAA".to_vec(), b"AAAG".to_vec()];
/// assert_eq!(check_error(b"AAAA", &whitelist, 1),
///            vec![b"AAAT".to_vec(), b"AAAG".to_vec()]);
/// ```
#[must_use]
pub fn check_error(barcode: &[u8], whitelist: &[Barcode], errors: u32) -> Vec<Barcode> {
    let mut near_matches = Vec::new();
    for candidate in whitelist {
        if candidate.as_slice() == barcode {
            continue;
        }
        let (long, short) = if barcode.len() > candidate.len() {
            (barcode.len(), candidate.len())
        } else {
            (candidate.len(), barcode.len())
        };
        if long > short + errors as usize {
            continue;
        }
        if matches_within_errors(barcode, candidate, errors) {
            near_matches.push(candidate.clone());
            if near_matches.len() > 1 {
                return near_matches;
            }
        }
    }
    near_matches
}

/// What to do with a whitelisted barcode that looks like an error of another
/// (Python source: the `--ed-above-threshold` choice).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResolutionMethod {
    /// Drop it from the whitelist.
    Discard,
    /// Fold it into its neighbour when the difference is a pure substitution;
    /// drop it when an indel is involved.
    Correct,
}

impl std::str::FromStr for ResolutionMethod {
    type Err = PortError;

    /// # Errors
    ///
    /// Returns [`PortError::Assertion`] with Python's message
    /// `"resolution method must be discard or correct"`.
    fn from_str(s: &str) -> Result<Self> {
        match s {
            "discard" => Ok(ResolutionMethod::Discard),
            "correct" => Ok(ResolutionMethod::Correct),
            _ => Err(PortError::assertion(
                "resolution method must be discard or correct",
            )),
        }
    }
}

/// Removes or corrects whitelisted barcodes that look like errors of other
/// whitelisted barcodes
/// (Python source: `whitelist_methods.errorDetectAboveThreshold`).
///
/// Barcodes are examined in ascending count order, each compared only against
/// those with a higher count, so a barcode is only ever folded into a more
/// abundant one.
///
/// # A mirrored upstream defect
///
/// Python's `error_counter["error_discarded_mt_1"]` is a bare subscript with
/// no increment. Because the counter is a `defaultdict(int)` the expression
/// merely creates the key, so the "discarded due to possible errors from
/// multiple other CBs" line always reports 0 however many barcodes were
/// dropped for that reason. The counter is reproduced with the same defect,
/// because the log line it feeds is compared by the test suite. See
/// `05_design.md` § Intentional deviations.
///
/// # Errors
///
/// This function is infallible in practice; the [`Result`] exists so callers
/// can propagate a resolution-method parse failure through the same path.
pub fn error_detect_above_threshold(
    counts: &CellCounts,
    cell_whitelist: &[Barcode],
    true_to_false_map: Option<&TrueToFalseMap>,
    errors: u32,
    resolution_method: ResolutionMethod,
) -> Result<(BTreeSet<Barcode>, TrueToFalseMap)> {
    let mut new_map: TrueToFalseMap = true_to_false_map.cloned().unwrap_or_default();
    let mut discard: BTreeSet<Barcode> = BTreeSet::new();

    let mut substitution_corrected = 0u64;
    let mut indel_discarded = 0u64;
    // Mirrors Python's never-incremented counter; see the note above.
    let error_discarded_mt_1 = 0u64;

    // Stable, so equal counts keep the caller's order — most-abundant-first,
    // which is what Python's list carries.
    let mut ordered: Vec<Barcode> = cell_whitelist.to_vec();
    ordered.sort_by_key(|barcode| counts.get(barcode).unwrap_or(0));

    for ix in 0..ordered.len() {
        let cb = &ordered[ix];
        let near_misses = check_error(cb, &ordered[ix + 1..], errors);

        if !near_misses.is_empty() {
            discard.insert(cb.clone());
        }

        if resolution_method == ResolutionMethod::Correct && near_misses.len() == 1 {
            // Only pure substitutions are corrected: an indel shifts the UMI
            // as well, so fixing the cell barcode alone would not repair the
            // read.
            if matches_within_substitutions(cb, &near_misses[0], errors) {
                new_map
                    .entry(near_misses[0].clone())
                    .or_default()
                    .insert(cb.clone());
                substitution_corrected += 1;
            } else {
                discard.insert(cb.clone());
                indel_discarded += 1;
            }
        }
    }

    if resolution_method == ResolutionMethod::Correct {
        utilities::info(&format!(
            "CBs above the knee corrected due to possible substitutions: {substitution_corrected}"
        ));
        utilities::info(&format!(
            "CBs above the knee discarded due to possible INDELs: {indel_discarded}"
        ));
        utilities::info(&format!(
            "CBs above the knee discarded due to possible errors from \
             multiple other CBs: {error_discarded_mt_1}"
        ));
    } else {
        utilities::info(&format!(
            "CBs above the knee discarded due to possible errors: {}",
            discard.len()
        ));
    }

    let remaining: BTreeSet<Barcode> = cell_whitelist
        .iter()
        .filter(|b| !discard.contains(*b))
        .cloned()
        .collect();
    Ok((remaining, new_map))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bc(s: &str) -> Barcode {
        s.as_bytes().to_vec()
    }

    /// Builds a counter in the order given, so tie behaviour is exercised.
    fn counts_of(pairs: &[(&str, u64)]) -> CellCounts {
        let mut out = CellCounts::new();
        for (barcode, count) in pairs {
            out.set(barcode.as_bytes(), *count);
        }
        out
    }

    fn bimodal() -> CellCounts {
        counts_of(&[
            ("AAAA", 1000),
            ("AAAT", 800),
            ("AATT", 600),
            ("ATTT", 400),
            ("TTTT", 200),
            ("CCCC", 10),
            ("CCCG", 8),
            ("CCGG", 6),
            ("CGGG", 4),
            ("GGGG", 2),
        ])
    }

    #[test]
    fn most_common_orders_by_descending_count() {
        let ranked = most_common(&counts_of(&[("A", 1), ("B", 3), ("C", 2)]));
        assert_eq!(ranked, vec![(bc("B"), 3), (bc("C"), 2), (bc("A"), 1)]);
    }

    /// `Counter.most_common` is stable over insertion order, so the barcode
    /// seen first wins a tie — not the alphabetically smaller one.
    #[test]
    fn most_common_breaks_ties_on_first_seen_order() {
        let ranked = most_common(&counts_of(&[("B", 5), ("A", 5)]));
        assert_eq!(ranked[0].0, bc("B"));
        let ranked = most_common(&counts_of(&[("A", 5), ("B", 5)]));
        assert_eq!(ranked[0].0, bc("A"));
    }

    #[test]
    fn observing_accumulates_and_records_order() {
        let mut counts = CellCounts::new();
        counts.observe(b"TTTT");
        counts.observe(b"AAAA");
        counts.observe(b"TTTT");
        assert_eq!(counts.get(b"TTTT"), Some(2));
        assert_eq!(counts.get(b"AAAA"), Some(1));
        assert_eq!(counts.get(b"GGGG"), None);
        assert_eq!(counts.barcodes(), [bc("TTTT"), bc("AAAA")]);
        assert_eq!(counts.len(), 2);
    }

    #[test]
    fn knee_distance_needs_two_points() {
        assert!(knee_distance(&[]).is_none());
        assert!(knee_distance(&[1]).is_none());
        assert!(knee_distance(&[1, 1]).is_some());
    }

    #[test]
    fn distance_knee_matches_python() {
        let result = knee_estimate_distance(&bimodal(), None).unwrap();
        let selected: Vec<Barcode> = result.barcodes.unwrap().into_iter().collect();
        assert_eq!(
            selected,
            vec![bc("AAAA"), bc("AAAT"), bc("AATT"), bc("ATTT")]
        );
        assert_eq!(result.knee_index, Some(3));
        assert_eq!(result.cumulative[..5], [1000, 1800, 2400, 2800, 3000]);
    }

    #[test]
    fn distance_knee_honours_an_explicit_cell_number() {
        let result = knee_estimate_distance(&bimodal(), Some(3)).unwrap();
        // The 4th most abundant count is 400, so barcodes above 400 are kept.
        assert_eq!(result.len(), 3);
    }

    #[test]
    fn cell_number_zero_is_treated_as_absent() {
        let with_zero = knee_estimate_distance(&bimodal(), Some(0)).unwrap();
        let without = knee_estimate_distance(&bimodal(), None).unwrap();
        assert_eq!(with_zero.barcodes, without.barcodes);
    }

    #[test]
    fn cell_number_past_the_end_is_an_error() {
        assert!(knee_estimate_distance(&bimodal(), Some(999)).is_err());
    }

    /// On a cumulative curve that is (almost) a straight line, every point is
    /// nominally zero distance from the chord and the winner is decided by
    /// floating-point residuals around 1e-16. Python behaves the same way:
    /// `getKneeEstimateDistance` on `{AAAA: 1000000, CCCC: 1, GGGG: 1}`
    /// returns all three barcodes rather than raising, because the residuals
    /// grow along the line and `argmax` lands on the last index, not the
    /// first. The `"Something's gone wrong here!!"` branch is therefore
    /// effectively unreachable for real inputs; this test pins the behaviour
    /// that actually occurs.
    #[test]
    fn straight_cumulative_curve_is_decided_by_float_residuals() {
        let counts = counts_of(&[("AAAA", 1_000_000), ("CCCC", 1), ("GGGG", 1)]);
        let result = knee_estimate_distance(&counts, None).unwrap();
        assert_eq!(result.len(), 3);
        assert_eq!(result.knee_index, Some(2));
    }

    /// The residual ordering must match `NumPy`'s, since it is what selects the
    /// knee on a degenerate curve. These are the values `NumPy` produces for the
    /// same input.
    #[test]
    fn knee_distance_residuals_match_numpy_on_a_straight_line() {
        let (dist, idx) = knee_distance(&[1_000_000, 1_000_001, 1_000_002]).unwrap();
        assert_eq!(idx, 2);
        assert_eq!(dist[0], 0.0);
        assert!(dist[1] > 0.0 && dist[1] < 1e-15);
        assert!(dist[2] > dist[1]);
    }

    #[test]
    fn distance_knee_needs_at_least_two_barcodes() {
        let err = knee_estimate_distance(&counts_of(&[("A", 1)]), None).unwrap_err();
        assert!(matches!(err, PortError::Value(_)));
    }

    #[test]
    fn kde_covariance_matches_scipy_convention() {
        let data = [1.0, 2.0, 3.0, 4.0, 10.0];
        let kde = GaussianKde::new(&data, 0.1).unwrap();
        let expected = numeric::variance(&data, 1) * 0.01;
        assert!((kde.covariance() - expected).abs() < 1e-15);
    }

    #[test]
    fn kde_rejects_degenerate_datasets() {
        assert!(GaussianKde::new(&[1.0], 0.1).is_err());
        assert!(GaussianKde::new(&[2.0, 2.0, 2.0], 0.1).is_err());
    }

    #[test]
    fn kde_integrates_to_about_one() {
        let data = [0.0, 1.0, 2.0, 3.0, 4.0];
        let kde = GaussianKde::new(&data, 1.0).unwrap();
        let xs = numeric::linspace(-10.0, 14.0, 5000);
        let step = xs[1] - xs[0];
        let mass: f64 = kde.eval_many(&xs).iter().sum::<f64>() * step;
        assert!((mass - 1.0).abs() < 1e-3, "mass = {mass}");
    }

    #[test]
    fn density_knee_matches_python_on_a_bimodal_distribution() {
        let result = knee_estimate_density(&bimodal(), None, None).unwrap();
        let selected: Vec<Barcode> = result.barcodes.unwrap().into_iter().collect();
        assert_eq!(
            selected,
            vec![bc("AAAA"), bc("AAAT"), bc("AATT"), bc("ATTT")]
        );
    }

    /// The 0.1%-of-maximum filter can leave a single count, at which point
    /// `SciPy`'s `gaussian_kde` raises ``ValueError: `dataset` input should have
    /// multiple elements.`` — confirmed by running this exact distribution
    /// through the Python original.
    #[test]
    fn density_knee_errors_when_the_abundance_filter_leaves_one_count() {
        let counts = counts_of(&[
            ("AAAA", 10000),
            ("AAAT", 5),
            ("AATT", 4),
            ("ATTT", 3),
            ("TTTT", 2),
            ("CCCC", 1),
        ]);
        let err = knee_estimate_density(&counts, None, None).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.contains("at least two data points")));
    }

    /// A distribution with several distinct abundances but no local minimum
    /// passing the heuristics yields no whitelist at all, which is what
    /// `whitelist` reports as "No local minima was accepted".
    #[test]
    fn density_knee_returns_none_when_no_minimum_qualifies() {
        let counts = counts_of(&[
            ("AAAA", 100),
            ("AAAT", 99),
            ("AATT", 98),
            ("ATTT", 97),
            ("TTTT", 96),
            ("CCCC", 95),
        ]);
        let result = knee_estimate_density(&counts, None, None).unwrap();
        assert!(result.barcodes.is_none());
        assert!(result.is_empty());
    }

    #[test]
    fn density_knee_rejects_a_zero_variance_distribution() {
        let counts = counts_of(&[("AAAA", 100), ("CCCC", 100), ("GGGG", 100)]);
        assert!(knee_estimate_density(&counts, None, None).is_err());
    }

    #[test]
    fn knee_method_parsing() {
        assert_eq!(
            "distance".parse::<KneeMethod>().unwrap(),
            KneeMethod::Distance
        );
        assert_eq!(
            "density".parse::<KneeMethod>().unwrap(),
            KneeMethod::Density
        );
        let err = "other".parse::<KneeMethod>().unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m == "knee_method must be 'distance' or 'density'"));
    }

    #[test]
    fn error_correct_mapping_matches_python() {
        let counts = counts_of(&[
            ("AAAA", 1000),
            ("AAAT", 800),
            ("AATT", 600),
            ("ATTT", 400),
            ("TTTT", 200),
            ("CCCC", 10),
            ("CCCG", 8),
        ]);
        let observed: Vec<Barcode> = counts.barcodes().to_vec();
        let whitelist: BTreeSet<Barcode> = [bc("AAAA"), bc("CCCC")].into_iter().collect();
        let map = get_error_correct_mapping(&observed, &whitelist, 1);
        assert_eq!(map[&bc("AAAA")], [bc("AAAT")].into_iter().collect());
        assert_eq!(map[&bc("CCCC")], [bc("CCCG")].into_iter().collect());
        assert_eq!(map.len(), 2);
    }

    #[test]
    fn error_correct_mapping_drops_ambiguous_barcodes() {
        // AAAT is one substitution from both AAAA and AAAG.
        let whitelist: BTreeSet<Barcode> = [bc("AAAA"), bc("AAAG")].into_iter().collect();
        let map = get_error_correct_mapping(&[bc("AAAT")], &whitelist, 1);
        assert!(map.is_empty());
    }

    #[test]
    fn check_error_matches_python() {
        let whitelist = vec![bc("AAAT"), bc("AATT"), bc("AAAA"), bc("AAAG")];
        assert_eq!(
            check_error(b"AAAA", &whitelist, 1),
            vec![bc("AAAT"), bc("AAAG")]
        );
        assert!(check_error(b"AAAA", &[bc("CCCC"), bc("GGGG")], 1).is_empty());
        assert!(check_error(b"ACGT", &[], 1).is_empty());
    }

    #[test]
    fn check_error_returns_early_after_two_matches() {
        let whitelist = vec![bc("ACGA"), bc("ACGC"), bc("ACGG")];
        assert_eq!(check_error(b"ACGT", &whitelist, 1).len(), 2);
    }

    #[test]
    fn check_error_skips_incompatible_lengths() {
        assert!(check_error(b"AAAA", &[bc("AA")], 1).is_empty());
        assert_eq!(check_error(b"AAAA", &[bc("AAA")], 1), vec![bc("AAA")]);
    }

    #[test]
    fn error_detect_discard_matches_python() {
        let counts = counts_of(&[
            ("AAAA", 1000),
            ("AAAT", 800),
            ("AATT", 600),
            ("ATTT", 400),
            ("TTTT", 200),
        ]);
        // Most-abundant-first, as the estimator returns them.
        let whitelist = vec![bc("AAAA"), bc("AAAT"), bc("AATT"), bc("ATTT")];
        let (out, _) =
            error_detect_above_threshold(&counts, &whitelist, None, 1, ResolutionMethod::Discard)
                .unwrap();
        assert_eq!(out.into_iter().collect::<Vec<_>>(), vec![bc("AAAA")]);
    }

    #[test]
    fn error_detect_correct_builds_the_correction_chain() {
        let counts = counts_of(&[
            ("AAAA", 1000),
            ("AAAT", 800),
            ("AATT", 600),
            ("ATTT", 400),
            ("TTTT", 200),
        ]);
        let whitelist = vec![bc("AAAA"), bc("AAAT"), bc("AATT"), bc("ATTT")];
        let mut seed: TrueToFalseMap = TrueToFalseMap::new();
        seed.insert(bc("ATTT"), [bc("TTTT")].into_iter().collect());

        let (out, map) = error_detect_above_threshold(
            &counts,
            &whitelist,
            Some(&seed),
            1,
            ResolutionMethod::Correct,
        )
        .unwrap();
        assert_eq!(out.into_iter().collect::<Vec<_>>(), vec![bc("AAAA")]);
        assert_eq!(map[&bc("AAAA")], [bc("AAAT")].into_iter().collect());
        assert_eq!(map[&bc("AAAT")], [bc("AATT")].into_iter().collect());
        assert_eq!(map[&bc("AATT")], [bc("ATTT")].into_iter().collect());
        assert_eq!(map[&bc("ATTT")], [bc("TTTT")].into_iter().collect());
    }

    #[test]
    fn resolution_method_parsing() {
        assert_eq!(
            "discard".parse::<ResolutionMethod>().unwrap(),
            ResolutionMethod::Discard
        );
        assert_eq!(
            "correct".parse::<ResolutionMethod>().unwrap(),
            ResolutionMethod::Correct
        );
        assert!("other".parse::<ResolutionMethod>().is_err());
    }

    #[test]
    fn position_tuples_enumerate_with_repetition() {
        assert_eq!(position_tuples(2, 1), vec![vec![0], vec![1]]);
        assert_eq!(
            position_tuples(2, 2),
            vec![vec![0, 0], vec![0, 1], vec![1, 0], vec![1, 1]]
        );
        assert_eq!(position_tuples(3, 0), vec![Vec::<usize>::new()]);
        assert!(position_tuples(0, 1).is_empty());
    }

    fn write_tsv(dir: &Path, name: &str, body: &str) -> std::path::PathBuf {
        let path = dir.join(name);
        std::fs::write(&path, body).unwrap();
        path
    }

    #[test]
    fn user_defined_plain_whitelist() {
        let dir = tempfile::tempdir().unwrap();
        let p = write_tsv(
            dir.path(),
            "wl.tsv",
            "# comment\nAAAA\tAAAT,AATT\nCCCC\tCCCG\n",
        );
        let (wl, map) = get_user_defined_barcodes(&p, None, ErrorCorrectionMode::None, 1).unwrap();
        assert_eq!(
            wl.into_iter().collect::<Vec<_>>(),
            vec![bc("AAAA"), bc("CCCC")]
        );
        assert!(map.is_none());
    }

    #[test]
    fn user_defined_error_correction_from_file() {
        let dir = tempfile::tempdir().unwrap();
        let p = write_tsv(
            dir.path(),
            "wl.tsv",
            "# comment\nAAAA\tAAAT,AATT\nCCCC\tCCCG\n",
        );
        let (wl, map) =
            get_user_defined_barcodes(&p, None, ErrorCorrectionMode::FromFile, 1).unwrap();
        assert_eq!(wl.len(), 2);
        let map = map.unwrap();
        assert_eq!(map[&bc("AAAT")], Some(bc("AAAA")));
        assert_eq!(map[&bc("AATT")], Some(bc("AAAA")));
        assert_eq!(map[&bc("CCCG")], Some(bc("CCCC")));
        assert_eq!(map.len(), 3);
    }

    #[test]
    fn user_defined_paired_whitelists_take_the_product() {
        let dir = tempfile::tempdir().unwrap();
        let p1 = write_tsv(dir.path(), "wl1.tsv", "AAAA\tAAAT\nCCCC\tCCCG\n");
        let p2 = write_tsv(dir.path(), "wl2.tsv", "GG\nTT\n");
        let (wl, _) =
            get_user_defined_barcodes(&p1, Some(&p2), ErrorCorrectionMode::None, 1).unwrap();
        assert_eq!(
            wl.into_iter().collect::<Vec<_>>(),
            vec![bc("AAAAGG"), bc("AAAATT"), bc("CCCCGG"), bc("CCCCTT")]
        );
    }

    #[test]
    fn user_defined_paired_with_file_correction_is_an_assertion_error() {
        let dir = tempfile::tempdir().unwrap();
        let p1 = write_tsv(dir.path(), "wl1.tsv", "AAAA\tAAAT\n");
        let p2 = write_tsv(dir.path(), "wl2.tsv", "GG\n");
        let err = get_user_defined_barcodes(&p1, Some(&p2), ErrorCorrectionMode::FromFile, 1)
            .unwrap_err();
        assert!(matches!(err, PortError::Assertion(_)));
    }

    #[test]
    fn user_defined_derived_correction_matches_python() {
        let dir = tempfile::tempdir().unwrap();
        let p = write_tsv(dir.path(), "wl.tsv", "AC\nGT\n");
        let (wl, map) =
            get_user_defined_barcodes(&p, None, ErrorCorrectionMode::Derive, 1).unwrap();
        assert_eq!(wl.into_iter().collect::<Vec<_>>(), vec![bc("AC"), bc("GT")]);

        let map = map.unwrap();
        // Reachable from AC only.
        assert_eq!(map[&bc("AA")], Some(bc("AC")));
        assert_eq!(map[&bc("TC")], Some(bc("AC")));
        // Reachable from both, so ambiguous.
        assert_eq!(map[&bc("AT")], None);
        assert_eq!(map[&bc("GC")], None);
        // The whitelist barcodes themselves are not keys.
        assert!(!map.contains_key(&bc("AC")));
        assert!(!map.contains_key(&bc("GT")));
        assert_eq!(map.len(), 14);
    }

    #[test]
    fn derived_correction_rejects_bases_outside_atcg() {
        let dir = tempfile::tempdir().unwrap();
        let p = write_tsv(dir.path(), "wl.tsv", "AN\n");
        let err = get_user_defined_barcodes(&p, None, ErrorCorrectionMode::Derive, 1).unwrap_err();
        assert!(matches!(err, PortError::Format(ref m) if m.contains("'N'")));
    }

    #[test]
    fn user_defined_missing_second_column_is_a_format_error() {
        let dir = tempfile::tempdir().unwrap();
        let p = write_tsv(dir.path(), "wl.tsv", "AAAA\n");
        let err =
            get_user_defined_barcodes(&p, None, ErrorCorrectionMode::FromFile, 1).unwrap_err();
        assert!(matches!(err, PortError::Format(_)));
    }

    #[test]
    fn get_cell_whitelist_runs_error_correction_only_when_asked() {
        let counts = bimodal();
        let (_, map) = get_cell_whitelist(&counts, KneeMethod::Distance, None, None, 0).unwrap();
        assert!(map.is_none());

        let (result, map) =
            get_cell_whitelist(&counts, KneeMethod::Distance, None, None, 1).unwrap();
        assert_eq!(result.len(), 4);
        let map = map.unwrap();
        assert_eq!(map[&bc("ATTT")], [bc("TTTT")].into_iter().collect());
    }

    #[test]
    fn cell_thresholds_tsv_matches_pythons_layout() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("plots");
        let result = knee_estimate_distance(&bimodal(), None).unwrap();
        write_cell_thresholds(&result, KneeMethod::Distance, &prefix, None).unwrap();
        let body = std::fs::read_to_string(dir.path().join("plots_cell_thresholds.tsv")).unwrap();
        assert_eq!(body, "count\n3\n");

        let result = knee_estimate_density(&bimodal(), None, None).unwrap();
        write_cell_thresholds(&result, KneeMethod::Density, &prefix, None).unwrap();
        let body = std::fs::read_to_string(dir.path().join("plots_cell_thresholds.tsv")).unwrap();
        assert!(body.starts_with("count\taction\n"));
        assert!(body.contains("\tSelected\n"));
    }

    #[test]
    fn cell_thresholds_tsv_is_skipped_for_an_explicit_cell_number() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("plots");
        let result = knee_estimate_distance(&bimodal(), Some(3)).unwrap();
        write_cell_thresholds(&result, KneeMethod::Distance, &prefix, Some(3)).unwrap();
        assert!(!dir.path().join("plots_cell_thresholds.tsv").exists());
    }
}
