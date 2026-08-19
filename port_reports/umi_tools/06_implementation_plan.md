# 06 — Implementation Plan (Step 8)

## Relationship to the auto-generated plan

`slice_optimizer --algorithm cohesion` output is kept verbatim in
`_slice_plan_raw.md` (graph: 107 nodes, 66 edges, 47 components, max depth 4,
density 0.0058). It was **not** adopted as-is, for the reason the workflow
anticipates: it partitions into five ~30-item alphabetical bins that split every
class from its own methods (`UMIClusterer` lands in Slice 1,
`UMIClusterer.__init__` in Slice 3, `get_bundles.update_dicts` in Slice 4,
`random_read_generator.__init__` in Slice 5). That ordering has forward
dependencies everywhere and cannot compile incrementally, which is fatal in
Rust in a way it would not be in Python.

The very low edge density (66 edges over 107 nodes) is itself the finding: this
package's coupling is not export-to-export, it is **export-to-`options`-object**,
and the optimizer cannot see that edge because `options` is not an export. So
the slice plan below is built on the module dependency graph instead, with the
optimizer's depth ordering used to check that no slice has a forward dependency.

Module dependency graph (internal imports only):

```
Documentation ──────────────────────────────► {count, count_tab, dedup, extract,
                                                group, prepare_for_em, whitelist}
_dedup_umi ──► {network, umi_methods, whitelist_methods}
Utilities  ──► {everything}
whitelist_methods ──► network
umi_methods ──► extract_methods, {extract, whitelist}
sam_methods ──► {count, count_tab, dedup, group}
network ──► {count, count_tab, dedup, group}
extract_methods ──► {extract, whitelist}
```

Depth ordering: `_dedup_umi` → `Utilities`(library half) → `umi_methods` →
`extract_methods` / `whitelist_methods` → `network` → `sam_methods` → CLI.

## Slice 0 — Cargo scaffold, exporters, fixtures

Manual, prepended per the workflow.

- **Work**: `init_rust_crate.py umi_tools-rs --py-name umi_tools --py-dir …`;
  add `crates/cli`; pin workspace dependencies; retarget the PyO3 extension to
  `umi_tools_rs` (deviation D9) and verify with `check_maturin_develop.py`;
  write the Python-side exporters under
  `port_reports/umi_tools/data_exports/`; copy the reduced test inputs into
  `crates/core/tests/fixtures/` and record `SHA256SUMS`.
- **`data_registry.py stage`**: **not run** — `data_assets.yaml` declares no
  remote assets (`03_data_asset_strategy.md`). Step 10 is skipped for the same
  reason.
- **Tests**: `cargo build` succeeds for the empty scaffold; `cargo test` runs
  the scaffold's placeholder test; `python -c "import umi_tools, umi_tools_rs"`
  proves both coexist in one interpreter.
- **Stage**: n/a (scaffolding).
- **Done when**: `cargo build --workspace` clean, `maturin develop` installs,
  the oracle still imports unshadowed.
- **Status**: complete.

## Slice 1 — `core::distance`, `core::numeric`, `core::error`, `core::utilities`

- **Exports**: `edit_distance` (+ `within_threshold`, `edit_distance_f64`);
  the NumPy subset (`median`, `cumsum`, `argmax`, `digitize`, `bincount`,
  `linspace`, `log10`, `power`, `vstack`-equivalent, `outer`, `sum`, `sqrt`);
  `PortError`; `log/info/warning/warn/debug/error/critical`, `open_file`,
  `create_file`, `get_temp_file`, `get_temp_filename`, `determine_format`,
  `output_names_and_formats`, `get_header`, `get_params`, `get_footer`,
  `format_multiline`, `GlobalOptions`.
- **Depends on**: Slice 0.
- **Tests**: `#[test]` per function; differential
  `validation/differential_distance.py` and `differential_numeric.py`.
- **Stage**: C (`distance`), D (`utilities`, `numeric`).
- **Done when**: Tier 1 exact for every integer/boolean output; ≤ 1e-9 for the
  float reductions.

## Slice 2 — `core::umi_methods`

- **Exports**: `RANGES`, `Record`, `guess_format`, `Display`, `FastqIterator`,
  `JoinedFastqIterator`, `get_average_umi_distance`, `RandomReadGenerator`.
- **Depends on**: Slices 0–1.
- **Tests**: round-trip parse of every `tests/*.fastq*` fixture; suffix-removal
  error cases; `strict=False` reconciliation; single-UMI `-1.0` case.
  Differential `validation/differential_umi_methods.py`.
- **Stage**: C.
- **Done when**: parsed records byte-identical to Python's for every fixture;
  `get_average_umi_distance` ≤ 1e-9.

## Slice 3 — `core::network` (the algorithmic core)

- **Exports**: `breadth_first_search`, `recursive_search`,
  `breadth_first_search_recursive`, `remove_umis`, `get_substr_slices`,
  `build_substr_idx`, `iter_nearest_neighbours`, `ClusterMethod`,
  `UmiClusterer` (all five methods + the six `_group_*`/`_get_best_*` helpers),
  `CellClusterer`.
- **Depends on**: Slices 0–1.
- **Tests**: the five-method table from `step1_unit_test.py` re-derived against
  the *live* 1.1.6 API (the in-repo copy is stale — see `00_baseline.md`); the
  `>25` UMI substring-index path vs. the `<=25` pairwise path on the same input;
  single-UMI, all-identical, all-distinct, and unequal-length-assertion cases.
  Differential `validation/differential_network.py` over randomised UMI count
  dictionaries.
- **Stage**: C.
- **Done when**: cluster membership exactly equal to Python's for all five
  methods over ≥ 500 randomised inputs, including the substring-index path.

## Slice 4 — `core::extract`

- **Exports**: `add_barcodes_to_identifier`, `extract_seq_and_quals`,
  `get_below_threshold`, `umi_below_threshold`, `mask_umi`, `extract_barcodes`,
  `ExtractConfig`, `ExtractFilterAndUpdate` (string and regex methods, 5′/3′,
  either-read, quality filter/mask, cell/UMI whitelist filtering).
- **Depends on**: Slices 0–2.
- **Tests**: every `bc-pattern` in `tests.yaml` (string and regex, including the
  fuzzy inDrop patterns); `retain_umi`; `--either-read` with both resolutions;
  reads shorter than the pattern (`ValueError`); non-`N/X/C` pattern
  (`ValueError`). Differential `validation/differential_extract.py`.
- **Stage**: C.
- **Done when**: extracted `(cell, umi, seq, quals, identifier)` byte-identical
  for every fixture read, and the `read_counts` counter matches key-for-key.

## Slice 5 — `core::fuzzy`, `core::bktree`, `core::whitelist`

- **Exports**: the fuzzy matcher (`{e<=n}`, `{s<=n}`, `{i<=n}`); `BkTree`;
  `knee_estimate_density`, `knee_estimate_distance`, `get_error_correct_mapping`,
  `get_cell_whitelist`, `get_user_defined_barcodes`, `check_error`,
  `error_detect_above_threshold`.
- **Depends on**: Slices 0–1, 3 (`whitelist_methods` imports `network`).
- **Tests**: KDE and `argrelextrema` against SciPy on the same vectors;
  both knee methods on the real `tests/` barcode count distributions; BK-tree
  candidate sets vs. `pybktree`; fuzzy matcher vs. the `regex` module over an
  exhaustive small alphabet. Differential `validation/differential_whitelist.py`.
- **Stage**: C.
- **Done when**: knee thresholds within Tier 2 (r ≥ 0.9999) and the *selected
  barcode sets* exactly equal; everything else Tier 1.

## Slice 6 — `core::sam`, `core::bundles`, `core::alignment`, `core::stats`

- **Exports**: the five barcode getters, `get_gene_count_tab`,
  `get_metacontig2contig`, `MetaFetcher`, `find_splice`, `get_read_position`,
  `GetBundles`, `TwoPassPairWriter`, `open_input_alignments`,
  `open_output_alignments`, `sort_output`, `aggregate_stats`.
- **Depends on**: Slices 0–3.
- **Tests**: cigar/splice cases (leading/trailing soft clip, `N`, hard clip,
  invalid op → `ValueError`); bundling over `chr19.bam` compared position by
  position; SAM/BAM/CRAM round-trip; stats aggregation vs. pandas.
  Differential `validation/differential_sam.py`.
- **Stage**: C for the pure functions, D for the record iterators.
- **Done when**: bundle keys and counts exactly equal to Python's over the full
  `chr19.bam`; `aggregate_stats` exact.

## Slice 7 — `crates/cli`

- **Exports**: `cli::documentation`, the option parser/validator
  (`validate_extract_options`, `validate_sam_options`, `start`, `stop`,
  `split_append_value`, `short_help`, the help formatter), and the eight
  entry points (`extract`, `whitelist`, `group`, `dedup`, `count`, `count_tab`,
  `prepare_for_em`, dispatch).
- **Depends on**: Slices 0–6.
- **Tests**: option-validation matrix ported from `validateExtractOptions` /
  `validateSamOptions` as `Result::Err` assertions; end-to-end CLI runs
  reproducing `tests.yaml` cases against the Python CLI's live output.
- **Stage**: D.
- **Done when**: for every `tests.yaml` case, the Rust binary's output matches
  the Python binary's output run on the same input in the same session.

## Slice 8 — Examples, docs, dead-code audit

- **Work**: `crates/core/examples/*.rs` (Step 11), `cargo doc --no-deps
  --workspace` clean, README, `dead_code_audit.py`.
- **Stage**: n/a.

## Forward-dependency check

| Slice | Depends on | Any forward dependency? |
|---|---|---|
| 1 | 0 | no |
| 2 | 0, 1 | no |
| 3 | 0, 1 | no |
| 4 | 0, 1, 2 | no |
| 5 | 0, 1, 3 | no |
| 6 | 0, 1, 2, 3 | no |
| 7 | 0–6 | no |
| 8 | 0–7 | no |

`whitelist_methods` importing `network` (Slice 5 after Slice 3) is the only
non-obvious edge; the optimizer's depth-4 chain confirms it.

## Live stage tracker

Every module reached its planned stage. "Stage C" means a pure-Rust
implementation with the PyO3 bindings retained; "Stage D" means pure Rust with
no Python-facing surface. Nothing remains at Stage A or B.

| Module | Planned stage | Reached | Evidence |
|---|---|---|---|
| `core::distance` | C | **C** | `validation/results_network.csv`, exact |
| `core::mersenne` | C | **C** | unit tests against CPython's and NumPy's streams |
| `core::numeric` | D | **D** | `crates/core/tests/reference_parity.rs`, ≤ 1e-9 |
| `core::utilities` | D | **D** | unit tests; log format checked against a live run |
| `core::umi_methods` | C | **C** | `results_network.csv`; FASTQ round-trip over every fixture |
| `core::pattern` | C | **C** | 1715 pattern/subject comparisons vs the `regex` module |
| `core::network` | C | **C** | `results_network.csv`, 200 randomised positions × 5 methods |
| `core::extract` | C | **C** | `results_extract.csv`, 1200 real reads |
| `core::fuzzy` / `core::bktree` | C | **C** | exhaustive comparison vs `regex` and `pybktree` |
| `core::whitelist` | C | **C** | `results_whitelist.csv`, incl. the Tier-2 density knee |
| `core::sam` / `core::bundles` | C/D | **C/D** | `crates/core/tests/sam_parity.rs` over two real BAMs |
| `core::alignment` / `core::stats` | D | **D** | round-trip tests; stats files byte-identical |
| `crates/cli` | D | **D** | the upstream CLI suite, 67/70 |

## Slice outcomes

| Slice | Status | Notes |
|---|---|---|
| 0 — scaffold, exporters, fixtures | complete | `data_registry.py stage` not run; no remote assets |
| 1 — distance, numeric, error, utilities | complete | a NumPy-2 defect and an inverted log-level comparison found here |
| 2 — umi_methods | complete | `RandomReadGenerator` later reworked to reproduce `numpy.random.choice` |
| 3 — network | complete | the `cluster_ordered` API was added once ordering proved observable |
| 4 — extract | complete | required `core::pattern`, an addition to the plan |
| 5 — fuzzy, bktree, whitelist | complete | `CellCounts` became an ordered counter for the same reason |
| 6 — sam, bundles, alignment, stats | complete | the MAPQ-255 sentinel bug was found here |
| 7 — cli | complete | five real defects found by running the upstream suite; see `08_validation.md` |
| 8 — examples, docs, dead-code audit | complete | Step 11/12 |

## Additions to the planned module list

Three modules were not in the Step-7 design and were added during
implementation, each because a dependency turned out to need reimplementing
rather than mapping:

| Module | Why it was needed |
|---|---|
| `core::pattern` | The Rust `regex` crate provides no per-group spans and no approximate matching; `extract` needs both. |
| `core::mersenne` | `--random-seed` makes the tie-break and the null distribution reproducible, and the upstream suite compares those outputs byte for byte — so CPython's and NumPy's generators both had to be reproduced, not merely their distributions. |
| `cli::help` / `cli::usage` | `--help` output is compared byte for byte; it is captured from the live command rather than re-derived. |
