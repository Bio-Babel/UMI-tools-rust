//! Per-UMI deduplication statistics
//! (Python source: `umi_tools/dedup.py`'s `--output-stats` path).
//!
//! `dedup --output-stats` writes three TSVs that let a user judge whether
//! deduplication worked: the distribution of reads per UMI per position, a
//! per-UMI summary, and the average edit distance between UMIs at a position
//! compared against a null drawn from the observed UMI frequencies.
//!
//! The Python original builds these with pandas — a `DataFrame`, a
//! `groupby().agg()`, a `merge(how='left')` and a `fillna(0).astype(int)`.
//! That is the only pandas use in the whole package, over a two-column frame,
//! so it is reimplemented here rather than taking a dataframe dependency. The
//! column order, the sort order and the integer cast are all mirrored, because
//! they are visible in the files the test suite compares.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

use crate::error::{PortError, Result};
use crate::numeric;
use crate::utilities;

/// One UMI's summary across every position it was seen at
/// (Python source: the columns of `aggregateStatsDF`'s returned frame).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UmiStats {
    /// The median of the per-position counts, excluding zeros.
    ///
    /// Python's `groupby().agg('median')` produces a float here, and the
    /// column is only cast to `int` later, in `main`, after the pre/post merge
    /// — so a median of 1.5 becomes 1 in the output file, truncated toward
    /// zero rather than rounded. [`write_per_umi`] performs that cast.
    pub median_counts: f64,
    /// How many positions the UMI was observed at.
    pub times_observed: u64,
    /// The total reads across all those positions.
    pub total_counts: u64,
}

/// UMI to its summary, ordered by UMI — pandas sorts a `groupby` key by
/// default, and `merge(sort=True)` sorts again, so the output is
/// UMI-ordered.
pub type UmiStatsMap = BTreeMap<Vec<u8>, UmiStats>;

/// Aggregates per-position observations into per-UMI summaries
/// (Python source: `dedup.aggregateStatsDF`).
///
/// The input is the flat `(UMI, count)` list `dedup` accumulates, one entry
/// per UMI per position. Python does:
///
/// ```text
/// grouped = stats_df.groupby("UMI")
/// agg_df = grouped.agg({'counts': ['median', len, 'sum']})
/// agg_df.columns = ['median_counts', 'times_observed', 'total_counts']
/// ```
///
/// # Examples
///
/// ```
/// use umi_tools_core::stats::aggregate_stats;
///
/// let observations = vec![
///     (b"AAA".to_vec(), 1), (b"AAA".to_vec(), 3),
///     (b"CCC".to_vec(), 5),
/// ];
/// let stats = aggregate_stats(&observations);
/// assert_eq!(stats[&b"AAA".to_vec()].median_counts, 2.0);
/// assert_eq!(stats[&b"AAA".to_vec()].times_observed, 2);
/// assert_eq!(stats[&b"AAA".to_vec()].total_counts, 4);
/// ```
#[must_use]
pub fn aggregate_stats(observations: &[(Vec<u8>, u64)]) -> UmiStatsMap {
    let mut grouped: BTreeMap<Vec<u8>, Vec<u64>> = BTreeMap::new();
    for (umi, count) in observations {
        grouped.entry(umi.clone()).or_default().push(*count);
    }
    grouped
        .into_iter()
        .map(|(umi, counts)| {
            let as_f64: Vec<f64> = counts.iter().map(|&c| c as f64).collect();
            let stats = UmiStats {
                median_counts: numeric::median(&as_f64),
                times_observed: counts.len() as u64,
                total_counts: counts.iter().sum(),
            };
            (umi, stats)
        })
        .collect()
}

/// Writes the per-UMI-per-position count distribution
/// (Python source: the `_per_umi_per_position.tsv` block of `dedup.main`).
///
/// Each row is a read count and how many (UMI, position) pairs had that count,
/// before and after deduplication. The row set is the union of the counts
/// seen on either side, sorted ascending — so a count present only after
/// deduplication still gets a row, with `instances_pre` at zero.
///
/// # Errors
///
/// Returns [`PortError::Io`] if the file cannot be written.
pub fn write_per_umi_per_position(
    pre: &[(Vec<u8>, u64)],
    post: &[(Vec<u8>, u64)],
    path: &Path,
) -> Result<()> {
    let tally = |rows: &[(Vec<u8>, u64)]| -> BTreeMap<u64, u64> {
        let mut out = BTreeMap::new();
        for (_, count) in rows {
            *out.entry(*count).or_insert(0) += 1;
        }
        out
    };
    let pre_counts = tally(pre);
    let post_counts = tally(post);

    let mut index: Vec<u64> = pre_counts
        .keys()
        .chain(post_counts.keys())
        .copied()
        .collect();
    index.sort_unstable();
    index.dedup();

    let mut out = utilities::create_file(path, false)?;
    writeln!(out, "counts\tinstances_pre\tinstances_post").map_err(|e| PortError::io(path, e))?;
    for count in index {
        writeln!(
            out,
            "{count}\t{}\t{}",
            pre_counts.get(&count).copied().unwrap_or(0),
            post_counts.get(&count).copied().unwrap_or(0)
        )
        .map_err(|e| PortError::io(path, e))?;
    }
    out.finish()
}

/// Writes the per-UMI summary table
/// (Python source: the `_per_umi.tsv` block of `dedup.main`).
///
/// The pre- and post-deduplication summaries are joined on the UMI with a
/// **left** join from the pre side, so a UMI that survived deduplication but
/// was never seen before it cannot appear — and a UMI that disappeared gets
/// zeros on the post side. Python spells that as
/// `merge(..., how='left', sort=True).fillna(0).astype(int)`; the `astype(int)`
/// truncates the float medians toward zero, which is reproduced here.
///
/// # Errors
///
/// - [`PortError::Io`] if the file cannot be written.
/// - [`PortError::Value`] if a UMI is not valid UTF-8, where Python's
///   `x.decode()` raises `UnicodeDecodeError`.
pub fn write_per_umi(pre: &UmiStatsMap, post: &UmiStatsMap, path: &Path) -> Result<()> {
    let mut out = utilities::create_file(path, false)?;
    writeln!(
        out,
        "UMI\tmedian_counts_pre\ttimes_observed_pre\ttotal_counts_pre\
         \tmedian_counts_post\ttimes_observed_post\ttotal_counts_post"
    )
    .map_err(|e| PortError::io(path, e))?;

    for (umi, pre_stats) in pre {
        let name = utilities::decode_utf8(umi)?;
        let post_stats = post.get(umi).copied().unwrap_or(UmiStats {
            median_counts: 0.0,
            times_observed: 0,
            total_counts: 0,
        });
        writeln!(
            out,
            "{name}\t{}\t{}\t{}\t{}\t{}\t{}",
            pre_stats.median_counts as i64,
            pre_stats.times_observed,
            pre_stats.total_counts,
            post_stats.median_counts as i64,
            post_stats.times_observed,
            post_stats.total_counts
        )
        .map_err(|e| PortError::io(path, e))?;
    }
    out.finish()
}

/// The binned average-edit-distance table
/// (Python source: the `edit_distance_df` built at the end of `dedup.main`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EditDistanceTable {
    /// The bin labels, with the lowest relabelled `Single_UMI`.
    pub edit_distance: Vec<String>,
    /// Positions with one UMI, pre-deduplication.
    pub unique: Vec<usize>,
    /// The null expectation for `unique`.
    pub unique_null: Vec<usize>,
    /// Positions after deduplication by the chosen method.
    pub method: Vec<usize>,
    /// The null expectation for `method`.
    pub method_null: Vec<usize>,
}

/// Bins the average edit distances and tallies each bin
/// (Python source: the `bin_clusters`/`tallyCounts` helpers and the
/// `edit_distance_df` assembly in `dedup.main`).
///
/// Bins run from `-1` to `max + 1` inclusive. The `-1` bin catches positions
/// with a single UMI, whose average distance is the `-1.0` sentinel
/// [`crate::umi_methods::get_average_umi_distance`] returns; that bin is
/// relabelled `Single_UMI` in the output, which is what Python's
/// `edit_distance_df.loc[0, 'edit_distance'] = "Single_UMI"` does.
///
/// # Errors
///
/// Returns [`PortError::Value`] when every input series is empty, where
/// Python's `max(map(max, ...))` raises `ValueError`.
pub fn bin_edit_distances(
    pre: &[f64],
    post: &[f64],
    pre_null: &[f64],
    post_null: &[f64],
) -> Result<EditDistanceTable> {
    let max_of = |series: &[f64]| numeric::max(series);
    let max_ed = [pre, post, pre_null, post_null]
        .iter()
        .filter_map(|s| max_of(s))
        .fold(f64::NEG_INFINITY, f64::max);
    if !max_ed.is_finite() {
        return Err(PortError::value("max() arg is an empty sequence"));
    }
    let max_ed = max_ed as i64;

    // Python: `cluster_bins = range(-1, int(max_ed) + 2)`.
    let bins: Vec<f64> = (-1..=max_ed + 1).map(|b| b as f64).collect();
    let min_length = (max_ed + 3).max(0) as usize;

    let tally = |series: &[f64]| -> Vec<usize> {
        numeric::bincount(&numeric::digitize_right(series, &bins), min_length)
    };

    let mut edit_distance: Vec<String> = bins.iter().map(|b| format!("{}", *b as i64)).collect();
    // The lowest bin holds the single-UMI sentinel.
    if let Some(first) = edit_distance.first_mut() {
        *first = "Single_UMI".to_string();
    }

    Ok(EditDistanceTable {
        edit_distance,
        unique: tally(pre),
        unique_null: tally(pre_null),
        method: tally(post),
        method_null: tally(post_null),
    })
}

/// Writes the binned edit-distance table
/// (Python source: `edit_distance_df.to_csv(..., index=False, sep="\t")`).
///
/// `method_name` becomes two column headers, `<method>` and `<method>_null`,
/// exactly as Python's `options.method` does — so the file's header depends on
/// which clustering method was run.
///
/// # Errors
///
/// Returns [`PortError::Io`] if the file cannot be written.
pub fn write_edit_distance(
    table: &EditDistanceTable,
    method_name: &str,
    path: &Path,
) -> Result<()> {
    let mut out = utilities::create_file(path, false)?;
    writeln!(
        out,
        "unique\tunique_null\t{method_name}\t{method_name}_null\tedit_distance"
    )
    .map_err(|e| PortError::io(path, e))?;

    let rows = table.edit_distance.len();
    let at = |v: &[usize], i: usize| v.get(i).copied().unwrap_or(0);
    for i in 0..rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            at(&table.unique, i),
            at(&table.unique_null, i),
            at(&table.method, i),
            at(&table.method_null, i),
            table.edit_distance[i]
        )
        .map_err(|e| PortError::io(path, e))?;
    }
    out.finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn obs(pairs: &[(&str, u64)]) -> Vec<(Vec<u8>, u64)> {
        pairs
            .iter()
            .map(|(u, c)| (u.as_bytes().to_vec(), *c))
            .collect()
    }

    #[test]
    fn aggregation_matches_pandas() {
        let rows = obs(&[
            ("AAA", 1),
            ("AAA", 3),
            ("CCC", 5),
            ("GGG", 2),
            ("GGG", 4),
            ("GGG", 6),
        ]);
        let stats = aggregate_stats(&rows);
        assert_eq!(stats.len(), 3);
        assert_eq!(stats[&b"AAA".to_vec()].median_counts, 2.0);
        assert_eq!(stats[&b"AAA".to_vec()].times_observed, 2);
        assert_eq!(stats[&b"AAA".to_vec()].total_counts, 4);
        assert_eq!(stats[&b"CCC".to_vec()].median_counts, 5.0);
        assert_eq!(stats[&b"GGG".to_vec()].median_counts, 4.0);
        assert_eq!(stats[&b"GGG".to_vec()].total_counts, 12);
    }

    #[test]
    fn even_group_medians_are_averaged_like_pandas() {
        let stats = aggregate_stats(&obs(&[("AAA", 1), ("AAA", 2)]));
        assert_eq!(stats[&b"AAA".to_vec()].median_counts, 1.5);
    }

    #[test]
    fn aggregation_of_nothing_is_empty() {
        assert!(aggregate_stats(&[]).is_empty());
    }

    #[test]
    fn per_umi_per_position_unions_both_count_sets() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pos.tsv");
        // pre: two UMIs at count 1, one at 2; post: one at 1, two at 3.
        let pre = obs(&[("A", 1), ("B", 1), ("C", 2)]);
        let post = obs(&[("A", 1), ("B", 3), ("C", 3)]);
        write_per_umi_per_position(&pre, &post, &path).unwrap();

        let body = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = body.trim_end().split('\n').collect();
        assert_eq!(lines[0], "counts\tinstances_pre\tinstances_post");
        assert_eq!(lines[1], "1\t2\t1");
        assert_eq!(lines[2], "2\t1\t0");
        assert_eq!(lines[3], "3\t0\t2");
    }

    #[test]
    fn per_umi_left_joins_and_truncates_medians() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("umi.tsv");
        // AAA has a 1.5 median pre and nothing post; CCC survives.
        let pre = aggregate_stats(&obs(&[("AAA", 1), ("AAA", 2), ("CCC", 4)]));
        let post = aggregate_stats(&obs(&[("CCC", 4)]));
        write_per_umi(&pre, &post, &path).unwrap();

        let body = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = body.trim_end().split('\n').collect();
        assert!(lines[0].starts_with("UMI\tmedian_counts_pre"));
        // 1.5 truncates to 1, and the missing post side is zero-filled.
        assert_eq!(lines[1], "AAA\t1\t2\t3\t0\t0\t0");
        assert_eq!(lines[2], "CCC\t4\t1\t4\t4\t1\t4");
    }

    #[test]
    fn per_umi_rejects_non_utf8_umis() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("umi.tsv");
        let pre: UmiStatsMap = [(
            vec![0xff, 0xfe],
            UmiStats {
                median_counts: 1.0,
                times_observed: 1,
                total_counts: 1,
            },
        )]
        .into_iter()
        .collect();
        assert!(write_per_umi(&pre, &UmiStatsMap::new(), &path).is_err());
    }

    #[test]
    fn edit_distance_bins_run_from_the_single_umi_sentinel() {
        // -1.0 is the single-UMI sentinel; 0.0 and 2.0 are real distances.
        let table = bin_edit_distances(&[-1.0, 0.0, 2.0], &[-1.0], &[1.0], &[1.5]).unwrap();
        assert_eq!(table.edit_distance[0], "Single_UMI");
        assert_eq!(table.edit_distance[1], "0");
        assert_eq!(table.edit_distance.last().unwrap(), "3");
        // pre = [-1, 0, 2] bins to [0, 1, 3]; NumPy gives [1, 1, 0, 1, 0].
        assert_eq!(table.unique, vec![1, 1, 0, 1, 0]);
        // post = [-1] is the single-UMI bin alone.
        assert_eq!(table.method, vec![1, 0, 0, 0, 0]);
        // pre_null = [1.0] bins to index 2.
        assert_eq!(table.unique_null, vec![0, 0, 1, 0, 0]);
        // post_null = [1.5] bins to index 3.
        assert_eq!(table.method_null, vec![0, 0, 0, 1, 0]);
        assert_eq!(table.unique.len(), table.edit_distance.len());
    }

    #[test]
    fn edit_distance_of_nothing_is_a_value_error() {
        let err = bin_edit_distances(&[], &[], &[], &[]).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m == "max() arg is an empty sequence"));
    }

    #[test]
    fn edit_distance_file_names_columns_after_the_method() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("ed.tsv");
        let table = bin_edit_distances(&[-1.0, 1.0], &[-1.0], &[0.0], &[0.0]).unwrap();
        write_edit_distance(&table, "directional", &path).unwrap();

        let body = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = body.trim_end().split('\n').collect();
        assert_eq!(
            lines[0],
            "unique\tunique_null\tdirectional\tdirectional_null\tedit_distance"
        );
        assert!(lines[1].ends_with("\tSingle_UMI"));
        assert_eq!(lines.len(), table.edit_distance.len() + 1);
    }
}
