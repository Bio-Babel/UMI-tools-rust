//! Fixture-driven parity tests against reference values dumped from the live
//! Python package.
//!
//! The JSON files under `tests/fixtures/` are produced by the exporters in
//! `port_reports/umi_tools/data_exports/`, which run the *original* Python
//! implementation in a live interpreter — they are not hand-written
//! expectations. Regenerate them by running those scripts against the Python
//! oracle environment.
//!
//! These complement, rather than replace, the same-process differential
//! harnesses under `validation/`: these run in CI with no Python installed,
//! while the harnesses compare against a live oracle.

use std::collections::{BTreeMap, BTreeSet};

use umi_tools_core::fuzzy;
use umi_tools_core::network::{
    breadth_first_search, build_substr_idx, get_substr_slices, iter_nearest_neighbours,
    remove_umis, AdjMap, CellClusterer, ClusterMethod, CountMap, Umi, UmiClusterer,
};
use umi_tools_core::whitelist::{
    check_error, error_detect_above_threshold, get_error_correct_mapping, knee_distance,
    knee_estimate_distance, most_common, Barcode, CellCounts, ResolutionMethod, TrueToFalseMap,
};

fn load(name: &str) -> serde_json::Value {
    let raw = match name {
        "network" => include_str!("fixtures/reference_network.json"),
        "whitelist" => include_str!("fixtures/reference_whitelist.json"),
        other => panic!("unknown fixture {other}"),
    };
    serde_json::from_str(raw).expect("fixture parses")
}

fn umi(s: &str) -> Umi {
    s.as_bytes().to_vec()
}

fn strings(v: &serde_json::Value) -> Vec<String> {
    v.as_array()
        .expect("array")
        .iter()
        .map(|x| x.as_str().expect("string").to_string())
        .collect()
}

fn umis(v: &serde_json::Value) -> Vec<Umi> {
    strings(v).iter().map(|s| umi(s)).collect()
}

fn counts_from(v: &serde_json::Value) -> CountMap {
    v.as_object()
        .expect("object")
        .iter()
        .map(|(k, c)| (umi(k), c.as_u64().expect("count")))
        .collect()
}

#[test]
fn substr_slices_match_python() {
    let doc = load("network");
    let cases = doc["helpers"]["get_substr_slices"].as_array().unwrap();
    assert!(!cases.is_empty());
    for case in cases {
        let umi_length = case["umi_length"].as_u64().unwrap() as usize;
        let idx_size = case["idx_size"].as_u64().unwrap() as usize;
        let expected: Vec<(usize, usize)> = case["slices"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| {
                let pair = s.as_array().unwrap();
                (
                    pair[0].as_u64().unwrap() as usize,
                    pair[1].as_u64().unwrap() as usize,
                )
            })
            .collect();
        assert_eq!(
            get_substr_slices(umi_length, idx_size),
            expected,
            "get_substr_slices({umi_length}, {idx_size})"
        );
    }
}

#[test]
fn breadth_first_search_matches_python() {
    let doc = load("network");
    for case in doc["helpers"]["breadth_first_search"].as_array().unwrap() {
        let adj: AdjMap = case["adj"]
            .as_object()
            .unwrap()
            .iter()
            .map(|(k, v)| (umi(k), umis(v)))
            .collect();
        let node = umi(case["node"].as_str().unwrap());
        let got: Vec<Umi> = breadth_first_search(&node, &adj).into_iter().collect();
        assert_eq!(got, umis(&case["component"]), "bfs from {:?}", case["node"]);
    }
}

#[test]
fn remove_umis_matches_python() {
    let doc = load("network");
    let adj_case = &doc["helpers"]["breadth_first_search"].as_array().unwrap()[0];
    let adj: AdjMap = adj_case["adj"]
        .as_object()
        .unwrap()
        .iter()
        .map(|(k, v)| (umi(k), umis(v)))
        .collect();

    for case in doc["helpers"]["remove_umis"].as_array().unwrap() {
        let cluster = umis(&case["cluster"]);
        let nodes = umis(&case["nodes"]);
        let got: Vec<Umi> = remove_umis(&adj, &cluster, &nodes).into_iter().collect();
        assert_eq!(
            got,
            umis(&case["remaining"]),
            "remove_umis {:?}",
            case["nodes"]
        );
    }
}

#[test]
fn nearest_neighbour_pairs_match_python() {
    let doc = load("network");
    for case in doc["helpers"]["iter_nearest_neighbours"]
        .as_array()
        .unwrap()
    {
        let items = umis(&case["umis"]);
        let min_edit = case["min_edit"].as_u64().unwrap() as u32;
        let idx = build_substr_idx(&items, items[0].len(), min_edit);

        // Python yields unordered pairs; the fixture stores them sorted, so
        // both sides are normalised to a sorted pair set.
        let mut got: Vec<[Umi; 2]> = iter_nearest_neighbours(&items, &idx)
            .into_iter()
            .map(|(a, b)| if a <= b { [a, b] } else { [b, a] })
            .collect();
        got.sort();

        let mut expected: Vec<[Umi; 2]> = case["pairs"]
            .as_array()
            .unwrap()
            .iter()
            .map(|p| {
                let pair = p.as_array().unwrap();
                let (a, b) = (
                    umi(pair[0].as_str().unwrap()),
                    umi(pair[1].as_str().unwrap()),
                );
                if a <= b {
                    [a, b]
                } else {
                    [b, a]
                }
            })
            .collect();
        expected.sort();

        assert_eq!(got, expected, "iter_nearest_neighbours min_edit={min_edit}");
    }
}

/// The full clustering chain over 70 cases — 10 hand-chosen (including the
/// upstream repository's own worked example and both sides of the 25-UMI
/// substring-index boundary) and 60 seeded-random ones.
#[test]
fn clustering_matches_python_for_every_method() {
    let doc = load("network");
    let cases = doc["clustering"].as_array().unwrap();
    assert!(
        cases.len() >= 70,
        "expected the full case set, got {}",
        cases.len()
    );

    for case in cases {
        let name = case["name"].as_str().unwrap();
        let counts = counts_from(&case["counts"]);
        let threshold = case["threshold"].as_u64().unwrap() as u32;
        let items: Vec<Umi> = counts.keys().cloned().collect();

        for method in ClusterMethod::all() {
            let expected = &case["methods"][method.as_str()];
            let clusterer = UmiClusterer::new(method);

            // Stage 1: the adjacency list, compared with each neighbour list
            // sorted, since Python appends in pair-iteration order.
            let adj = clusterer.get_adj_list(&items, &counts, threshold);
            match (adj.graph(), expected["adj_list"].as_object()) {
                (None, None) => {}
                (Some(graph), Some(want)) => {
                    assert_eq!(graph.len(), want.len(), "{name}/{method}: adj_list size");
                    for (node, neighbours) in graph {
                        let key = String::from_utf8_lossy(node).into_owned();
                        let mut got: Vec<String> = neighbours
                            .iter()
                            .map(|n| String::from_utf8_lossy(n).into_owned())
                            .collect();
                        got.sort();
                        assert_eq!(
                            got,
                            strings(&want[&key]),
                            "{name}/{method}: neighbours of {key}"
                        );
                    }
                }
                (got, want) => panic!(
                    "{name}/{method}: adjacency presence differs (rust={:?}, python={:?})",
                    got.is_some(),
                    want.is_some()
                ),
            }

            // Stage 2: the final groups, exactly — including order, because
            // the exporter feeds Python the counts in sorted-key order to
            // match this crate's BTreeMap iteration. Only the trailing
            // members of an adjacency group are order-tolerant.
            let mut c = UmiClusterer::new(method);
            let got = c.cluster(&counts, threshold).expect("clustering succeeds");
            let want: Vec<Vec<Umi>> = expected["final"]
                .as_array()
                .unwrap()
                .iter()
                .map(umis)
                .collect();

            if method == ClusterMethod::Adjacency {
                assert_eq!(got.len(), want.len(), "{name}/{method}: group count");
                for (g, w) in got.iter().zip(want.iter()) {
                    assert_eq!(g[0], w[0], "{name}/{method}: group leader");
                    let gs: BTreeSet<&Umi> = g.iter().collect();
                    let ws: BTreeSet<&Umi> = w.iter().collect();
                    assert_eq!(gs, ws, "{name}/{method}: group membership");
                }
            } else {
                assert_eq!(got, want, "{name}/{method}: groups");
            }
        }
    }
}

#[test]
fn cell_clusterer_matches_python() {
    let doc = load("network");
    for case in doc["cell_clusterer"]["cases"].as_array().unwrap() {
        let fuzzy_match = case["fuzzy"].as_bool().unwrap();
        let counts = counts_from(&case["counts"]);
        let items: Vec<Umi> = counts.keys().cloned().collect();
        let threshold = case["dir_threshold"].as_f64().unwrap_or(10.0);

        let clusterer = CellClusterer::new("directional", threshold, fuzzy_match).unwrap();
        let adj = clusterer.adj_list_directional(&items, &counts);
        let want = case["adj_list"].as_object().unwrap();

        assert_eq!(
            adj.len(),
            want.len(),
            "adjacency size (fuzzy={fuzzy_match})"
        );
        for (node, neighbours) in &adj {
            let key = String::from_utf8_lossy(node).into_owned();
            let mut got: Vec<String> = neighbours
                .iter()
                .map(|n| String::from_utf8_lossy(n).into_owned())
                .collect();
            got.sort();
            assert_eq!(
                got,
                strings(&want[&key]),
                "neighbours of {key} (fuzzy={fuzzy_match})"
            );
        }
    }

    // The upstream `__call__` defect is recorded in the fixture too.
    let recorded = doc["cell_clusterer"]["call_error"].as_str().unwrap();
    assert!(
        recorded.starts_with("TypeError: "),
        "the Python original is expected to raise TypeError, got {recorded}"
    );
    let c = CellClusterer::new("directional", 10.0, false).unwrap();
    assert!(c.call(&[], &CountMap::new()).is_err());
}

/// Builds a counter from the fixture's JSON object.
///
/// The exporter dumps `dict(Counter)`, whose JSON object preserves insertion
/// order, and `serde_json` reads it in that order (its `preserve_order`
/// feature is off, so this relies on the map being read sequentially) — which
/// is the order `most_common` breaks ties in.
fn cell_counts_from(v: &serde_json::Value) -> CellCounts {
    let mut counts = CellCounts::new();
    for (barcode, count) in v.as_object().expect("object") {
        counts.set(barcode.as_bytes(), count.as_u64().expect("count"));
    }
    counts
}

#[test]
fn distance_knee_matches_python() {
    let doc = load("whitelist");
    let mut checked = 0;
    for case in doc["knee"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let counts = cell_counts_from(&case["counts"]);

        // The intermediate geometry first, so a divergence is localised.
        let cumulative: Vec<i64> = case["cumulative"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_i64().unwrap())
            .collect();
        let ranked = most_common(&counts);
        let raw: Vec<i64> = ranked.iter().map(|(_, c)| *c as i64).collect();
        let got_cumulative = umi_tools_core::numeric::cumsum_i64(&raw);
        assert_eq!(got_cumulative, cumulative, "{name}: cumulative counts");

        let (dist, idx) = knee_distance(&cumulative).expect("knee distance");
        let want_dist: Vec<f64> = case["dist_to_line"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap())
            .collect();
        assert_eq!(dist.len(), want_dist.len(), "{name}: distance curve length");
        for (g, w) in dist.iter().zip(want_dist.iter()) {
            // Tier 1: same algorithm, same operation order, no BLAS.
            assert!(
                (g - w).abs() <= 1e-9 * w.abs().max(1.0),
                "{name}: distance to line ({g} vs {w})"
            );
        }
        assert_eq!(
            idx,
            case["argmax"].as_u64().unwrap() as usize,
            "{name}: argmax of the distance curve"
        );

        if let Some(want) = case["distance_whitelist"].as_array() {
            let result = knee_estimate_distance(&counts, None).expect("distance knee");
            let got: Vec<String> = result
                .barcodes
                .clone()
                .unwrap()
                .iter()
                .map(|b| String::from_utf8_lossy(b).into_owned())
                .collect();
            let want: Vec<String> = want
                .iter()
                .map(|v| v.as_str().unwrap().to_string())
                .collect();

            // The knee takes the top `idx + 1` ranked barcodes. When the count
            // at the cut equals the count just past it, which barcodes land
            // inside is decided by the ranking's tie order — Python's
            // `Counter.most_common` breaks ties on the counter's insertion
            // order, this crate breaks them on the barcode. Both are
            // deterministic; neither is more correct, and the *number* of
            // barcodes and the threshold count are identical either way. See
            // 05_design.md § Determinism decisions.
            let idx = result.knee_index.expect("a knee index");
            let boundary = ranked[idx].1;
            let tied_at_cut = ranked.get(idx + 1).is_some_and(|(_, c)| *c == boundary);

            if tied_at_cut {
                assert_eq!(
                    got.len(),
                    want.len(),
                    "{name}: selected barcode count (tie at cut)"
                );
                for barcode in &got {
                    assert!(
                        counts.get(barcode.as_bytes()).unwrap_or(0) >= boundary,
                        "{name}: {barcode} selected below the cut count {boundary}"
                    );
                }
            } else {
                assert_eq!(got, want, "{name}: selected barcodes");
            }
            checked += 1;
        }
    }
    assert!(
        checked >= 5,
        "expected several distributions, checked {checked}"
    );
}

#[test]
fn density_kde_matches_scipy() {
    let doc = load("whitelist");
    let mut checked = 0;
    for case in doc["knee"].as_array().unwrap() {
        let Some(want_cov) = case["density_covariance"].as_f64() else {
            continue;
        };
        let name = case["name"].as_str().unwrap();
        let log_counts: Vec<f64> = case["density_log_counts"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap())
            .collect();

        let kde = umi_tools_core::whitelist::GaussianKde::new(&log_counts, 0.1).expect("kde fits");
        assert!(
            (kde.covariance() - want_cov).abs() <= 1e-12 * want_cov.abs().max(1.0),
            "{name}: kernel covariance ({} vs {want_cov})",
            kde.covariance()
        );

        let lo = log_counts.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = log_counts.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let xx = umi_tools_core::numeric::linspace(lo, hi, 10_000);
        let density = kde.eval_many(&xx);

        let want_sample: Vec<f64> = case["density_values_sample"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap())
            .collect();
        for (k, w) in want_sample.iter().enumerate() {
            let g = density[k * 500];
            assert!(
                (g - w).abs() <= 1e-9 * w.abs().max(1.0),
                "{name}: density at index {} ({g} vs {w})",
                k * 500
            );
        }

        let want_mins: Vec<usize> = case["local_mins"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap() as usize)
            .collect();
        assert_eq!(
            umi_tools_core::numeric::argrelmin(&density),
            want_mins,
            "{name}: local minima"
        );
        checked += 1;
    }
    assert!(
        checked >= 3,
        "expected several KDE cases, checked {checked}"
    );
}

#[test]
fn density_knee_selection_matches_python() {
    let doc = load("whitelist");
    let mut checked = 0;
    for case in doc["knee"].as_array().unwrap() {
        if case["density_covariance"].as_f64().is_none() {
            continue;
        }
        let name = case["name"].as_str().unwrap();
        let counts = cell_counts_from(&case["counts"]);
        let result = umi_tools_core::whitelist::knee_estimate_density(&counts, None, None);

        match (case["density_whitelist"].as_array(), result) {
            (Some(want), Ok(got)) => {
                let got: Vec<String> = got
                    .barcodes
                    .expect("a whitelist was selected")
                    .iter()
                    .map(|b| String::from_utf8_lossy(b).into_owned())
                    .collect();
                let want: Vec<String> = want
                    .iter()
                    .map(|v| v.as_str().unwrap().to_string())
                    .collect();
                assert_eq!(got, want, "{name}: density-selected barcodes");
                checked += 1;
            }
            (None, Ok(got)) => {
                assert!(
                    got.barcodes.is_none(),
                    "{name}: Python selected nothing but Rust selected {} barcodes",
                    got.len()
                );
                checked += 1;
            }
            (_, Err(e)) => panic!("{name}: density knee failed: {e}"),
        }
    }
    assert!(
        checked >= 3,
        "expected several density cases, checked {checked}"
    );
}

#[test]
fn error_correct_mapping_matches_python() {
    let doc = load("whitelist");
    for case in doc["error_correct_mapping"].as_array().unwrap() {
        let barcodes: Vec<Barcode> = strings(&case["barcodes"])
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect();
        let whitelist: BTreeSet<Barcode> = strings(&case["whitelist"])
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect();
        let threshold = case["threshold"].as_u64().unwrap() as u32;

        let got = get_error_correct_mapping(&barcodes, &whitelist, threshold);
        let want = case["mapping"].as_object().unwrap();
        assert_eq!(
            got.len(),
            want.len(),
            "mapping size (threshold={threshold})"
        );
        for (parent, children) in &got {
            let key = String::from_utf8_lossy(parent).into_owned();
            let got_children: Vec<String> = children
                .iter()
                .map(|c| String::from_utf8_lossy(c).into_owned())
                .collect();
            assert_eq!(got_children, strings(&want[&key]), "children of {key}");
        }
    }
}

#[test]
fn check_error_matches_python() {
    let doc = load("whitelist");
    for case in doc["check_error"].as_array().unwrap() {
        let barcode = case["barcode"].as_str().unwrap().as_bytes().to_vec();
        let whitelist: Vec<Barcode> = strings(&case["whitelist"])
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect();
        let errors = case["errors"].as_u64().unwrap() as u32;

        let got: Vec<String> = check_error(&barcode, &whitelist, errors)
            .iter()
            .map(|b| String::from_utf8_lossy(b).into_owned())
            .collect();
        assert_eq!(
            got,
            strings(&case["near_matches"]),
            "check_error({:?}, errors={errors})",
            case["barcode"]
        );
    }
}

#[test]
fn error_detect_above_threshold_matches_python() {
    let doc = load("whitelist");
    for case in doc["error_detect_above_threshold"].as_array().unwrap() {
        let counts = cell_counts_from(&case["counts"]);
        // Most-abundant-first, as the estimator returns them.
        let whitelist: Vec<Barcode> = {
            let selected: BTreeSet<Barcode> = strings(&case["whitelist_in"])
                .iter()
                .map(|s| s.as_bytes().to_vec())
                .collect();
            most_common(&counts)
                .into_iter()
                .filter(|(b, _)| selected.contains(b))
                .map(|(b, _)| b)
                .collect()
        };
        let seed: TrueToFalseMap = case["true_to_false_in"]
            .as_object()
            .unwrap()
            .iter()
            .map(|(k, v)| {
                (
                    k.as_bytes().to_vec(),
                    strings(v).iter().map(|s| s.as_bytes().to_vec()).collect(),
                )
            })
            .collect();
        let method: ResolutionMethod = case["resolution_method"].as_str().unwrap().parse().unwrap();

        let (out, map) =
            error_detect_above_threshold(&counts, &whitelist, Some(&seed), 1, method).unwrap();

        let got: Vec<String> = out
            .iter()
            .map(|b| String::from_utf8_lossy(b).into_owned())
            .collect();
        assert_eq!(
            got,
            strings(&case["whitelist_out"]),
            "{method:?}: whitelist"
        );

        let want_map = case["true_to_false_out"].as_object().unwrap();
        assert_eq!(map.len(), want_map.len(), "{method:?}: correction map size");
        for (parent, children) in &map {
            let key = String::from_utf8_lossy(parent).into_owned();
            let got_children: Vec<String> = children
                .iter()
                .map(|c| String::from_utf8_lossy(c).into_owned())
                .collect();
            assert_eq!(got_children, strings(&want_map[&key]), "{method:?}: {key}");
        }
    }
}

/// The three fuzzy patterns `umi_tools` builds, over every pattern/text pair
/// the exporter enumerated against the live PyPI `regex` module.
#[test]
fn fuzzy_matching_matches_python_regex() {
    let doc = load("whitelist");
    let cases = doc["fuzzy"].as_array().unwrap();
    assert!(
        cases.len() >= 100,
        "expected the full grid, got {}",
        cases.len()
    );

    for case in cases {
        let pattern = case["pattern"].as_str().unwrap().as_bytes();
        let text = case["text"].as_str().unwrap().as_bytes();
        let label = format!(
            "pattern={} text={}",
            case["pattern"].as_str().unwrap(),
            case["text"].as_str().unwrap()
        );
        assert_eq!(
            fuzzy::matches_within_errors(pattern, text, 1),
            case["e1"].as_bool().unwrap(),
            "{{e<=1}} {label}"
        );
        assert_eq!(
            fuzzy::matches_within_errors(pattern, text, 2),
            case["e2"].as_bool().unwrap(),
            "{{e<=2}} {label}"
        );
        assert_eq!(
            fuzzy::matches_within_substitutions(pattern, text, 1),
            case["s1"].as_bool().unwrap(),
            "{{s<=1}} {label}"
        );
        assert_eq!(
            fuzzy::matches_within_insertions(pattern, text, 1),
            case["i1"].as_bool().unwrap(),
            "{{i<=1}} {label}"
        );
    }
}

/// `BTreeMap`-keyed counts feed the clusterer in sorted order; this guards the
/// assumption the fixture comparison rests on.
#[test]
fn count_map_iterates_in_sorted_key_order() {
    let counts: BTreeMap<Umi, u64> = [(umi("TTTT"), 1), (umi("AAAA"), 2)].into_iter().collect();
    let keys: Vec<Umi> = counts.keys().cloned().collect();
    assert_eq!(keys, vec![umi("AAAA"), umi("TTTT")]);
}

/// Every `--bc-pattern` the package's test suite and documentation use, run
/// against 245 real reads from the checked-in FASTQ fixtures plus synthetic
/// edge cases — 1715 comparisons in all. Both whether the pattern matched and
/// the per-group spans are compared, because `core::extract` slices the read
/// by those spans.
#[test]
fn extraction_patterns_match_python_regex() {
    let doc: serde_json::Value =
        serde_json::from_str(include_str!("fixtures/reference_patterns.json"))
            .expect("fixture parses");
    let cases = doc["cases"].as_array().unwrap();
    assert!(
        cases.len() >= 1000,
        "expected the full grid, got {}",
        cases.len()
    );

    let mut matched = 0usize;
    for case in cases {
        let source = case["pattern"].as_str().unwrap();
        let subject = case["subject"].as_str().unwrap();
        let pattern = umi_tools_core::pattern::Pattern::parse(source)
            .unwrap_or_else(|e| panic!("pattern {source} failed to parse: {e}"));

        let got = pattern.match_at_start(subject);
        let want_matched = case["matched"].as_bool().unwrap();
        assert_eq!(
            got.is_some(),
            want_matched,
            "match/no-match for {source} on {subject:?}"
        );

        let Some(got) = got else { continue };
        matched += 1;

        let want_span = case["span"].as_array().unwrap();
        assert_eq!(
            got.span,
            (
                want_span[0].as_u64().unwrap() as usize,
                want_span[1].as_u64().unwrap() as usize
            ),
            "overall span for {source} on {subject:?}"
        );

        let want_groups = case["groups"].as_object().unwrap();
        assert_eq!(
            got.groups.len(),
            want_groups.len(),
            "captured group count for {source} on {subject:?}"
        );
        for (name, span) in want_groups {
            let span = span.as_array().unwrap();
            assert_eq!(
                got.span(name),
                Some((
                    span[0].as_u64().unwrap() as usize,
                    span[1].as_u64().unwrap() as usize
                )),
                "span of group {name} for {source} on {subject:?}"
            );
        }
    }
    assert!(matched > 100, "expected many matches, saw {matched}");
}

/// Builds the Rust config matching one exporter label. Kept beside the
/// exporter's `CONFIGS` table so the two stay in step.
fn extract_config_for(label: &str) -> umi_tools_core::extract::ExtractConfig {
    use umi_tools_core::extract::{ExtractConfig, ExtractMethod};
    use umi_tools_core::umi_methods::QualityEncoding;

    let base = ExtractConfig::new();
    match label {
        "string_simple" => ExtractConfig {
            pattern: Some("NNNXXXXNN".into()),
            ..base
        },
        "string_3prime" => ExtractConfig {
            pattern: Some("NNNNN".into()),
            prime3: true,
            ..base
        },
        "string_cell" => ExtractConfig {
            pattern: Some("CCCCCCNNNNNNNNNN".into()),
            extract_cell: true,
            ..base
        },
        "regex_5prime" => ExtractConfig {
            method: ExtractMethod::Regex,
            pattern: Some("^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})".into()),
            ..base
        },
        "regex_3prime" => ExtractConfig {
            method: ExtractMethod::Regex,
            pattern: Some(".*(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})$".into()),
            ..base
        },
        "regex_cell" => ExtractConfig {
            method: ExtractMethod::Regex,
            pattern: Some("^(?P<cell_1>.{6})(?P<umi_1>.{10})".into()),
            extract_cell: true,
            ..base
        },
        "regex_retain_umi" => ExtractConfig {
            method: ExtractMethod::Regex,
            pattern: Some("^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})".into()),
            retain_umi: true,
            ..base
        },
        "regex_quality_filter" => ExtractConfig {
            method: ExtractMethod::Regex,
            pattern: Some("^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})".into()),
            quality_encoding: Some(QualityEncoding::Phred33),
            quality_filter_threshold: Some(30),
            ..base
        },
        "regex_quality_mask" => ExtractConfig {
            method: ExtractMethod::Regex,
            pattern: Some("^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})".into()),
            quality_encoding: Some(QualityEncoding::Phred33),
            quality_filter_mask: Some(30),
            ..base
        },
        "regex_indrop_fuzzy" => ExtractConfig {
            method: ExtractMethod::Regex,
            pattern: Some(
                "(?P<cell_1>.{8,12})(?P<discard_2>GAGTGATTGCTTGTGACGCCTT{s<=2})\
                 (?P<cell_3>.{8})(?P<umi_1>.{6})T{3}.*"
                    .into(),
            ),
            extract_cell: true,
            ..base
        },
        "string_paired" => ExtractConfig {
            pattern: Some("NNN".into()),
            pattern2: Some("NNN".into()),
            ..base
        },
        "regex_scrb_paired" => ExtractConfig {
            method: ExtractMethod::Regex,
            pattern: Some("^(?P<cell_1>.{6})(?P<umi_1>.{10})".into()),
            extract_cell: true,
            ..base
        },
        other => panic!("no Rust config defined for exporter label {other}"),
    }
}

fn record_from(v: &serde_json::Value) -> umi_tools_core::Record {
    umi_tools_core::Record::new(
        v["identifier"].as_str().unwrap(),
        v["seq"].as_str().unwrap(),
        v["quals"].as_str().unwrap(),
    )
}

/// End-to-end extraction over 1360 real reads from the checked-in FASTQ
/// fixtures, across every `--bc-pattern`/`--extract-method` combination the
/// package's own test suite exercises — including the inDrop fuzzy pattern and
/// both quality modes.
///
/// The rewritten identifier, sequence and quality string are compared byte for
/// byte, as is the final `read_counts` tally that `umi_tools extract` prints.
#[test]
fn extraction_matches_python_on_real_reads() {
    let doc: serde_json::Value =
        serde_json::from_str(include_str!("fixtures/reference_extract.json"))
            .expect("fixture parses");

    let mut compared = 0usize;

    for group in doc["single"].as_array().unwrap() {
        let label = group["label"].as_str().unwrap();
        let mut extractor = umi_tools_core::ExtractFilterAndUpdate::new(extract_config_for(label))
            .unwrap_or_else(|e| panic!("{label}: config rejected: {e}"));

        for (i, case) in group["cases"].as_array().unwrap().iter().enumerate() {
            let mut read = record_from(&case["input"]);
            let got = extractor.process(&mut read, None);

            match (case["error"].as_str(), case["output"].as_object(), got) {
                (Some(py_err), _, Ok(_)) => {
                    panic!("{label}[{i}]: Python raised {py_err} but Rust succeeded")
                }
                (Some(_), _, Err(_)) => {}
                (None, Some(_), Ok(true)) => {
                    let want = &case["output"];
                    assert_eq!(
                        read.identifier,
                        want["identifier"].as_str().unwrap(),
                        "{label}[{i}]: identifier"
                    );
                    assert_eq!(read.seq, want["seq"].as_str().unwrap(), "{label}[{i}]: seq");
                    assert_eq!(
                        read.quals,
                        want["quals"].as_str().unwrap(),
                        "{label}[{i}]: quals"
                    );
                    compared += 1;
                }
                (None, None, Ok(false)) => compared += 1,
                (None, want, got) => panic!(
                    "{label}[{i}]: Python {} but Rust {:?}",
                    if want.is_some() {
                        "kept the read"
                    } else {
                        "discarded it"
                    },
                    got
                ),
            }
        }

        let want_counts = group["read_counts"].as_object().unwrap();
        for (event, count) in want_counts {
            assert_eq!(
                extractor
                    .read_counts()
                    .get(event.as_str())
                    .copied()
                    .unwrap_or(0),
                count.as_u64().unwrap(),
                "{label}: read_counts[{event}]"
            );
        }
        assert_eq!(
            extractor.read_counts().len(),
            want_counts.len(),
            "{label}: read_counts key set ({:?} vs {:?})",
            extractor.read_counts().keys().collect::<Vec<_>>(),
            want_counts.keys().collect::<Vec<_>>()
        );
    }

    for group in doc["paired"].as_array().unwrap() {
        let label = group["label"].as_str().unwrap();
        let mut extractor = umi_tools_core::ExtractFilterAndUpdate::new(extract_config_for(label))
            .unwrap_or_else(|e| panic!("{label}: config rejected: {e}"));

        for (i, case) in group["cases"].as_array().unwrap().iter().enumerate() {
            let mut r1 = record_from(&case["input1"]);
            let mut r2 = record_from(&case["input2"]);
            let got = extractor.process(&mut r1, Some(&mut r2));

            match (case["error"].as_str(), case["output1"].as_object(), got) {
                (Some(py_err), _, Ok(_)) => {
                    panic!("{label}[{i}]: Python raised {py_err} but Rust succeeded")
                }
                (Some(_), _, Err(_)) => {}
                (None, Some(_), Ok(true)) => {
                    for (got, want, side) in
                        [(&r1, &case["output1"], 1), (&r2, &case["output2"], 2)]
                    {
                        assert_eq!(
                            got.identifier,
                            want["identifier"].as_str().unwrap(),
                            "{label}[{i}]: read{side} identifier"
                        );
                        assert_eq!(
                            got.seq,
                            want["seq"].as_str().unwrap(),
                            "{label}[{i}]: read{side} seq"
                        );
                        assert_eq!(
                            got.quals,
                            want["quals"].as_str().unwrap(),
                            "{label}[{i}]: read{side} quals"
                        );
                    }
                    compared += 1;
                }
                (None, None, Ok(false)) => compared += 1,
                (None, want, got) => panic!(
                    "{label}[{i}]: Python {} but Rust {:?}",
                    if want.is_some() {
                        "kept the pair"
                    } else {
                        "discarded it"
                    },
                    got
                ),
            }
        }

        let want_counts = group["read_counts"].as_object().unwrap();
        for (event, count) in want_counts {
            assert_eq!(
                extractor
                    .read_counts()
                    .get(event.as_str())
                    .copied()
                    .unwrap_or(0),
                count.as_u64().unwrap(),
                "{label}: read_counts[{event}]"
            );
        }
    }

    assert!(
        compared > 1000,
        "expected the full read set, compared {compared}"
    );
}
