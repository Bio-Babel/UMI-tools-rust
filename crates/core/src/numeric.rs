//! The NumPy/SciPy subset the Python package actually uses.
//!
//! Python sources: `network.py`, `whitelist_methods.py`, `dedup.py`,
//! `umi_methods.py`.
//!
//! Every `NumPy` call in `umi_tools` is a one-dimensional reduction or an
//! elementwise operation over a `Vec<f64>`/`Vec<i64>`; none of them touches
//! BLAS, LAPACK, broadcasting across more than two ranks, or `dtype`
//! promotion between differently-sized integer types. Reimplementing that
//! subset here is smaller and easier to validate than adopting `ndarray`,
//! and it keeps summation order fixed and single-threaded so the divergence
//! against `NumPy` stays at Tier 1 rather than Tier 3.
//!
//! Each function documents the exact `NumPy` definition it reproduces,
//! including the boundary conventions that differ between plausible
//! implementations (median tie handling, `digitize(right=True)`'s
//! half-open interval, `bincount`'s `minlength`).

/// `numpy.median` over a slice: the middle element for odd lengths, the mean
/// of the two middle elements for even lengths.
///
/// Returns `f64::NAN` for an empty input, matching `NumPy`, which emits a
/// `RuntimeWarning` and produces `nan` rather than raising.
///
/// Used by `UMIClusterer._get_best_percentile` (`np.median(list(counts.values()))/100`)
/// and by the per-UMI statistics aggregation in `dedup.py`.
///
/// # Examples
///
/// ```
/// use umi_tools_core::numeric::median;
/// assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
/// assert_eq!(median(&[1.0, 2.0, 3.0, 4.0]), 2.5);
/// ```
#[must_use]
pub fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        // NumPy averages the two central order statistics; it does not pick
        // the lower one.
        f64::midpoint(sorted[n / 2 - 1], sorted[n / 2])
    }
}

/// `numpy.cumsum`: the running total, left to right.
///
/// Summation order is fixed and sequential, exactly as `NumPy`'s default
/// (non-pairwise for `cumsum`), so the float result matches bit for bit for
/// inputs of the magnitudes this package produces.
#[must_use]
pub fn cumsum(values: &[f64]) -> Vec<f64> {
    let mut acc = 0.0;
    values
        .iter()
        .map(|v| {
            acc += v;
            acc
        })
        .collect()
}

/// `numpy.cumsum` over integers, used for cumulative barcode counts.
#[must_use]
pub fn cumsum_i64(values: &[i64]) -> Vec<i64> {
    let mut acc = 0i64;
    values
        .iter()
        .map(|v| {
            acc += v;
            acc
        })
        .collect()
}

/// `numpy.argmax`: the index of the **first** maximal element.
///
/// The "first" is load-bearing — `getKneeEstimateDistance` uses the returned
/// index as a cut point, so a last-wins implementation would select a
/// different barcode set on a plateau.
///
/// Returns `None` for an empty input, where `NumPy` raises `ValueError`.
#[must_use]
pub fn argmax(values: &[f64]) -> Option<usize> {
    let mut best: Option<(usize, f64)> = None;
    for (i, &v) in values.iter().enumerate() {
        // Strictly greater, so the *first* maximum wins — NumPy's rule, and
        // load-bearing for the knee estimator, which uses the index as a cut
        // point.
        let replace = match best {
            None => true,
            Some((_, current)) => v > current,
        };
        if replace {
            best = Some((i, v));
        }
    }
    best.map(|(i, _)| i)
}

/// `numpy.linspace(start, stop, num)`: `num` evenly spaced samples,
/// inclusive of both endpoints.
///
/// Reproduces `NumPy`'s endpoint handling exactly: the last element is set to
/// `stop` rather than accumulated, so `linspace(a, b, n).last() == b` holds
/// without rounding drift. `num == 1` yields `[start]`, and `num == 0` an
/// empty vector.
///
/// Used by `getKneeEstimateDensity` with `num = 10_000`.
#[must_use]
pub fn linspace(start: f64, stop: f64, num: usize) -> Vec<f64> {
    match num {
        0 => Vec::new(),
        1 => vec![start],
        _ => {
            let step = (stop - start) / ((num - 1) as f64);
            let mut out: Vec<f64> = (0..num).map(|i| start + step * (i as f64)).collect();
            // NumPy assigns the endpoint directly rather than accumulating to it.
            out[num - 1] = stop;
            out
        }
    }
}

/// `numpy.digitize(values, bins, right=True)` for monotonically increasing
/// `bins`.
///
/// With `right=True` `NumPy` returns the index `i` satisfying
/// `bins[i-1] < x <= bins[i]`, which is exactly a binary search for the
/// leftmost insertion point (`numpy.searchsorted(bins, x, side='left')`);
/// verified against the live oracle rather than inferred from the docs.
///
/// Used by `dedup.py`'s `bin_clusters` to bucket average edit distances.
#[must_use]
pub fn digitize_right(values: &[f64], bins: &[f64]) -> Vec<usize> {
    values
        .iter()
        .map(|&x| bins.partition_point(|&b| b < x))
        .collect()
}

/// `numpy.bincount(values, minlength)`: counts of each non-negative integer.
///
/// The output length is `max(minlength, max(values) + 1)`, and is
/// `minlength` for an empty input — the same rule `NumPy` applies.
///
/// Used by `dedup.py`'s `tallyCounts` with `minlength = max_edit_distance + 3`.
#[must_use]
pub fn bincount(values: &[usize], minlength: usize) -> Vec<usize> {
    let needed = values.iter().copied().max().map_or(0, |m| m + 1);
    let mut out = vec![0usize; needed.max(minlength)];
    for &v in values {
        out[v] += 1;
    }
    out
}

/// Elementwise base-10 logarithm, as `numpy.log10`.
///
/// Non-positive inputs yield `-inf` (for zero) or `NaN` (for negatives),
/// matching `NumPy`, which warns rather than raising.
#[must_use]
pub fn log10(values: &[f64]) -> Vec<f64> {
    values.iter().map(|v| v.log10()).collect()
}

/// Sum of a slice, left to right — `numpy.sum` for a 1-D array.
///
/// `NumPy` uses pairwise summation for large arrays, which can differ from a
/// sequential sum in the last ulp. The vectors this package sums are counts
/// and probabilities of at most a few thousand elements, where the two agree;
/// the difference is bounded by the Tier-1 threshold and is recorded in
/// `08_validation.md`.
#[must_use]
pub fn sum(values: &[f64]) -> f64 {
    values.iter().sum()
}

/// Maximum of a slice, or `None` when empty (`NumPy` raises `ValueError`).
#[must_use]
pub fn max(values: &[f64]) -> Option<f64> {
    values
        .iter()
        .copied()
        .reduce(|a, b| if b > a { b } else { a })
}

/// Minimum of a slice, or `None` when empty (`NumPy` raises `ValueError`).
#[must_use]
pub fn min(values: &[f64]) -> Option<f64> {
    values
        .iter()
        .copied()
        .reduce(|a, b| if b < a { b } else { a })
}

/// `scipy.signal.argrelextrema(values, numpy.less)` with the default
/// `order=1`, `mode='clip'`.
///
/// Returns the indices of strict local minima: positions `i` where
/// `values[i] < values[i-1]` **and** `values[i] < values[i+1]`. Under
/// `mode='clip'` the boundary comparisons clamp the neighbour index, so
/// index 0 compares against itself on the left and can never satisfy the
/// strict inequality; the first and last indices are therefore never
/// returned. A plateau contains no strict minimum and is likewise never
/// returned — which is why `getKneeEstimateDensity` can find no local minimum
/// at all and fall through to `final_barcodes = None`.
///
/// # Examples
///
/// ```
/// use umi_tools_core::numeric::argrelmin;
/// assert_eq!(argrelmin(&[3.0, 1.0, 2.0, 0.5, 4.0]), vec![1, 3]);
/// assert_eq!(argrelmin(&[1.0, 1.0, 1.0]), Vec::<usize>::new());
/// ```
#[must_use]
pub fn argrelmin(values: &[f64]) -> Vec<usize> {
    if values.len() < 3 {
        return Vec::new();
    }
    (1..values.len() - 1)
        .filter(|&i| values[i] < values[i - 1] && values[i] < values[i + 1])
        .collect()
}

/// Population (biased) or sample (unbiased) variance of a slice.
///
/// `ddof` is the delta degrees of freedom: `0` reproduces
/// `numpy.var(x)`, `1` reproduces `numpy.cov(x, bias=False)` for a single
/// variable — the convention `scipy.stats.gaussian_kde` uses when computing
/// its data covariance.
///
/// Returns `f64::NAN` when `values.len() <= ddof`.
#[must_use]
pub fn variance(values: &[f64], ddof: usize) -> f64 {
    let n = values.len();
    if n <= ddof {
        return f64::NAN;
    }
    let mean = sum(values) / (n as f64);
    let ss: f64 = values.iter().map(|v| (v - mean) * (v - mean)).sum();
    ss / ((n - ddof) as f64)
}

/// Pearson correlation coefficient between two equal-length slices.
///
/// Used by the differential-test harness to score Tier-2 comparisons; kept in
/// the library so `cargo test` can apply the same criterion without the
/// Python side.
///
/// Returns `None` when the slices differ in length, are shorter than two
/// elements, or either has zero variance.
#[must_use]
pub fn pearson_r(a: &[f64], b: &[f64]) -> Option<f64> {
    if a.len() != b.len() || a.len() < 2 {
        return None;
    }
    let n = a.len() as f64;
    let ma = sum(a) / n;
    let mb = sum(b) / n;
    let mut num = 0.0;
    let mut da = 0.0;
    let mut db = 0.0;
    for (x, y) in a.iter().zip(b.iter()) {
        let dx = x - ma;
        let dy = y - mb;
        num += dx * dy;
        da += dx * dx;
        db += dy * dy;
    }
    if da == 0.0 || db == 0.0 {
        return None;
    }
    Some(num / (da.sqrt() * db.sqrt()))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn median_odd_and_even() {
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&[1.0, 2.0, 3.0, 4.0]), 2.5);
        assert_eq!(median(&[1.0]), 1.0);
    }

    #[test]
    fn median_of_empty_is_nan() {
        assert!(median(&[]).is_nan());
    }

    #[test]
    fn argmax_returns_first_maximum() {
        assert_eq!(argmax(&[1.0, 5.0, 5.0, 2.0]), Some(1));
        assert_eq!(argmax(&[]), None);
        assert_eq!(argmax(&[-3.0, -1.0, -2.0]), Some(1));
    }

    #[test]
    fn linspace_endpoints_are_exact() {
        let xs = linspace(1.0, 3.0, 10_000);
        assert_eq!(xs.len(), 10_000);
        assert_eq!(xs[0], 1.0);
        assert_eq!(xs[9_999], 3.0);
        assert_eq!(linspace(0.0, 1.0, 1), vec![0.0]);
        assert!(linspace(0.0, 1.0, 0).is_empty());
    }

    #[test]
    fn digitize_matches_left_insertion_point() {
        let bins: Vec<f64> = (-1..5).map(f64::from).collect();
        assert_eq!(
            digitize_right(&[-1.0, 0.0, 0.5, 1.0, 2.0, 3.0], &bins),
            vec![0, 1, 2, 2, 3, 4]
        );
    }

    #[test]
    fn bincount_honours_minlength_and_max() {
        assert_eq!(bincount(&[0, 1, 1, 3], 6), vec![1, 2, 0, 1, 0, 0]);
        assert_eq!(bincount(&[0, 1, 1, 3], 2), vec![1, 2, 0, 1]);
        assert_eq!(bincount(&[], 3), vec![0, 0, 0]);
        assert!(bincount(&[], 0).is_empty());
    }

    #[test]
    fn argrelmin_ignores_boundaries_and_plateaus() {
        assert_eq!(argrelmin(&[3.0, 1.0, 2.0, 0.5, 4.0]), vec![1, 3]);
        assert_eq!(argrelmin(&[1.0, 1.0, 1.0]), Vec::<usize>::new());
        assert_eq!(argrelmin(&[0.0, 1.0, 2.0]), Vec::<usize>::new());
        assert_eq!(argrelmin(&[1.0, 0.0]), Vec::<usize>::new());
        assert_eq!(argrelmin(&[]), Vec::<usize>::new());
    }

    #[test]
    fn variance_ddof_conventions() {
        let v = [1.0, 2.0, 3.0, 4.0];
        assert!((variance(&v, 0) - 1.25).abs() < EPS);
        assert!((variance(&v, 1) - (5.0 / 3.0)).abs() < EPS);
        assert!(variance(&[1.0], 1).is_nan());
    }

    #[test]
    fn pearson_r_of_identical_series_is_one() {
        let a = [1.0, 2.0, 3.0, 4.0];
        assert!((pearson_r(&a, &a).unwrap() - 1.0).abs() < EPS);
        assert_eq!(pearson_r(&a, &[1.0, 1.0, 1.0, 1.0]), None);
        assert_eq!(pearson_r(&a, &[1.0]), None);
    }

    /// Reference values dumped from the live Python oracle by
    /// `port_reports/umi_tools/data_exports/export_distance_numeric.py`.
    #[test]
    fn matches_numpy_reference_values() {
        let raw = include_str!("../tests/fixtures/reference_distance_numeric.json");
        let doc: serde_json::Value = serde_json::from_str(raw).expect("fixture parses");

        let as_vec = |v: &serde_json::Value| -> Vec<f64> {
            v.as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_f64().unwrap())
                .collect()
        };

        for case in doc["numeric"].as_array().unwrap() {
            let input = as_vec(&case["input"]);
            assert!(
                (median(&input) - case["median"].as_f64().unwrap()).abs() < EPS,
                "median {input:?}"
            );
            let expected_cumsum = as_vec(&case["cumsum"]);
            for (got, want) in cumsum(&input).iter().zip(expected_cumsum.iter()) {
                assert!((got - want).abs() < EPS, "cumsum {input:?}");
            }
            assert_eq!(
                argmax(&input),
                Some(usize::try_from(case["argmax"].as_u64().unwrap()).unwrap()),
                "argmax {input:?}"
            );
            assert!((sum(&input) - case["sum"].as_f64().unwrap()).abs() < 1e-6);
            assert_eq!(max(&input), Some(case["max"].as_f64().unwrap()));
            assert_eq!(min(&input), Some(case["min"].as_f64().unwrap()));
            if let Some(want) = case["log10"].as_array() {
                for (got, w) in log10(&input).iter().zip(want.iter()) {
                    assert!((got - w.as_f64().unwrap()).abs() < EPS, "log10 {input:?}");
                }
            }
        }

        for case in doc["linspace"].as_array().unwrap() {
            let lo = case["lo"].as_f64().unwrap();
            let hi = case["hi"].as_f64().unwrap();
            let n = usize::try_from(case["n"].as_u64().unwrap()).unwrap();
            let xs = linspace(lo, hi, n);
            assert_eq!(xs[0], case["first"].as_f64().unwrap());
            assert_eq!(xs[xs.len() - 1], case["last"].as_f64().unwrap());
            for (got, want) in xs.iter().zip(as_vec(&case["head"]).iter()) {
                assert!((got - want).abs() < EPS, "linspace({lo}, {hi}, {n})");
            }
        }

        for case in doc["digitize_bincount"].as_array().unwrap() {
            let values = as_vec(&case["values"]);
            let bins = as_vec(&case["bins"]);
            let binned = digitize_right(&values, &bins);
            let want: Vec<usize> = case["digitize_right"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| usize::try_from(v.as_u64().unwrap()).unwrap())
                .collect();
            assert_eq!(binned, want, "digitize {values:?}");

            let want_min6: Vec<usize> = case["bincount_minlength_6"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| usize::try_from(v.as_u64().unwrap()).unwrap())
                .collect();
            assert_eq!(bincount(&binned, 6), want_min6);

            let want_default: Vec<usize> = case["bincount_default"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| usize::try_from(v.as_u64().unwrap()).unwrap())
                .collect();
            assert_eq!(bincount(&binned, 0), want_default);
        }
    }
}
