# 08 — Validation (Steps 9 and 12)

## Headline

| Check | Result |
|---|---|
| Upstream CLI acceptance suite, run **unchanged** against the Rust binary | **67 / 70 passing** |
| Same-process differential harnesses (`validation/`) | **40 / 40 comparisons pass** |
| Python-vs-Rust CLI comparison, case by case, same session | **66 / 66 comparable cases byte-identical**; 3 blocked by documented gaps |
| `cargo test --workspace` | **395 passed, 0 failed** |
| `cargo clippy --workspace --all-targets -- -D warnings` | clean |
| `cargo fmt --all --check` | clean |
| `cargo doc --no-deps --workspace` | clean, no missing-doc or broken-link warnings |
| `maturin build --release` | produces an installable wheel |
| Dead-code audit | 0 findings in every category |

The 3 non-passing upstream cases are the two capability gaps in
§ Known gaps below, not behaviour differences.

## Which validation method applies, and why

`references/differential-testing-and-validation.md` allows either a
differential harness or `cargo test` coverage parity plus documented
deviations. This port uses **both**, plus a third method the package's own
structure makes available:

1. **Same-process differential harnesses** (`validation/differential_*.py`).
   The original Python package and the `maturin develop`-installed Rust
   extension are imported into one interpreter and called on the same inputs.
   No serialisation round-trip, so a difference cannot be an artefact of the
   comparison.
2. **Fixture-driven parity tests** (`crates/core/tests/*.rs`). Reference values
   dumped from the live Python package by the exporters in
   `port_reports/umi_tools/data_exports/`, compared by `cargo test` with no
   Python present — so CI needs neither a Python environment nor a network.
3. **The upstream acceptance suite, run against the Rust binary.** The
   package's own `tests/test_umi_tools.py` resolves `umi_tools` from `PATH`,
   so putting the Rust binary first runs its 70 end-to-end cases against the
   port with no modification to the harness. This is the strongest available
   evidence and the reason the port reaches byte-level parity rather than
   approximate parity.

Method 3 exists because of what `01_audit.md` found: the Python package has
**no library-level unit tests at all**. Its suite is entirely CLI-level. The
scaffolding script reached the same conclusion independently, classifying
`test_umi_tools.py` as "CLI/subprocess-driven — run unchanged against the
ported CLI rather than translated".

## Differential results

Recorded via `ResultRecorder` in `validation/results_*.csv`.

### `results_network.csv` — 10/10

| Function | Metric | Value | Threshold | Tier |
|---|---|---:|---:|:--:|
| `edit_distance` | mismatches | 0 | 0 | 1 |
| `edit_distance` unequal lengths | returns `None` | yes | — | 1 |
| `get_adj_list[adjacency]` | mismatched positions | 0 | 0 | 1 |
| `get_adj_list[directional]` | mismatched positions | 0 | 0 | 1 |
| `UMIClusterer[unique]` | mismatched positions | 0 / 200 | 0 | 1 |
| `UMIClusterer[percentile]` | mismatched positions | 0 / 200 | 0 | 1 |
| `UMIClusterer[cluster]` | mismatched positions | 0 / 200 | 0 | 1 |
| `UMIClusterer[adjacency]` | mismatched positions | 0 / 200 | 0 | 1 |
| `UMIClusterer[directional]` | mismatched positions | 0 / 200 | 0 | 1 |
| `get_average_umi_distance` | max abs diff | 0.0 | 1e-9 | 1 |

200 seeded-random UMI count dictionaries, 2–60 UMIs each, spanning both the
pairwise and the substring-index code paths (which diverge above 25 UMIs).

### `results_extract.csv` — 20/20

Every `--bc-pattern` in the package's test suite and documentation, over 150
real reads per configuration — 1200 reads in total — comparing the rewritten
identifier, sequence and quality string, and the `read_counts` tally.

| Configuration | Mismatched reads | Differing event counts |
|---|---:|---:|
| `string_simple` (`NNNXXXXNN`) | 0 / 150 | 0 |
| `string_3prime` | 0 / 150 | 0 |
| `string_cell` (`CCCCCCNNNNNNNNNN`) | 0 / 150 | 0 |
| `regex_5prime` | 0 / 150 | 0 |
| `regex_3prime` | 0 / 150 | 0 |
| `regex_retain_umi` | 0 / 150 | 0 |
| `regex_quality_mask` | 0 / 150 | 0 |
| `regex_indrop_fuzzy` (`{s<=2}`) | 0 / 150 | 0 |

Plus `mask_umi`, `umi_below_threshold`, `extractSeqAndQuals` and
`match_pattern` (the regex group spans, over 372 pattern/subject pairs) — all
exact.

### `results_whitelist.csv` — 6/6

| Function | Metric | Value | Threshold | Tier |
|---|---|---:|---:|:--:|
| `getKneeEstimateDistance` | mismatched distributions | 0 / 5 | 0 | 1 |
| `getKneeEstimateDensity` | mismatched distributions | 0 / 5 | 0 | **2** |
| `getErrorCorrectMapping` | mismatched distributions | 0 / 5 | 0 | 1 |
| `checkError` | mismatches | 0 | 0 | 1 |
| `getUserDefinedBarcodes` | mismatches | 0 | 0 | 1 |
| `errorDetectAboveThreshold` | mismatched distributions | 0 / 10 | 0 | 1 |

The density knee is the one genuinely Tier-2 comparison — a different numeric
backend for one step, since Python's kernel density estimate is SciPy's and
this crate's is its own. The *selected barcode sets* nonetheless match
exactly, which is the property that matters: a knee threshold either includes
a barcode or it does not.

The underlying curve is checked more tightly than a Tier-2 threshold would
require. `crates/core/tests/reference_parity.rs` compares the kernel
covariance to **1e-12 relative** and the density at 20 points across the
evaluation grid to **1e-9 relative** against SciPy 1.17.1, and requires the
set of local minima to be identical. So the Tier-2 allowance is not being
leaned on.

### `results_cli.csv` — end-to-end

Runs every case in `tests/tests.yaml` twice in the same session, once per
implementation, and compares the outputs each case declares. This compares the
two implementations against **each other**, not against checked-in references,
so a difference cannot be blamed on a stale fixture.

| Metric | Value |
|---|---:|
| Comparable cases producing identical output | **66 / 66** |
| Cases differing in output | **0** |
| Unexpected exit-code differences | **0** |
| Cases blocked by a documented gap | 3 |

Normalisation matches the upstream harness: `#`-prefixed and timestamped log
lines are dropped (both implementations stamp their own runtime, host and pid),
`.gz` outputs are decompressed before comparison (a gzip header embeds a
modification time), and cases declaring `sort: True` are sorted first — the
upstream suite does all three for the same reasons.

## Fixture-driven parity tests

Run by `cargo test` with no Python present. Reference values come from the
exporters in `port_reports/umi_tools/data_exports/`, each of which runs the
live Python package.

| Suite | What it compares | Volume |
|---|---|---|
| `reference_parity.rs` — clustering | adjacency list **and** final groups, all five methods | 70 cases × 5 methods |
| `reference_parity.rs` — helpers | `get_substr_slices`, BFS, `remove_umis`, nearest-neighbour pairs | all recorded cases |
| `reference_parity.rs` — knee | cumulative curve, distance-to-line curve, argmax, selected set | 6 distributions |
| `reference_parity.rs` — KDE | covariance, density at 20 grid points, local minima | 4 distributions |
| `reference_parity.rs` — whitelist | error-correct mapping, `checkError`, `errorDetectAboveThreshold`, user-defined barcodes | all recorded cases |
| `reference_parity.rs` — fuzzy | `{e<=n}`, `{s<=n}`, `{i<=n}` vs the `regex` module | 180 comparisons |
| `reference_parity.rs` — patterns | match/no-match **and** per-group spans | 1715 comparisons |
| `reference_parity.rs` — extraction | rewritten records and event tallies | 1360 reads |
| `sam_parity.rs` — coordinates | grouping position, flush coordinate, splice offset | 800 reads, 2 real BAMs |
| `sam_parity.rs` — bundling | bundle contents and read-event tallies | 3 option scenarios |

## Coverage parity

`scripts/test_coverage_diff.py` reports 3 Python tests against 355 Rust tests.
That comparison needs the mapping in `test_coverage.md` § Interpretation: the
Python side is three `pytest` functions, one parameterised into 70 cases, and
the script counts functions. Each is covered:

| Python test | Rust counterpart |
|---|---|
| `test_script[...]` (70 cases) | the same suite, run unchanged against the Rust binary — 67/70 |
| `test_script_has_main[script]` | `cli::run`'s dispatch plus `cli::tests::usage_lists_every_tool` |
| `test_style` (`pycodestyle`) | `cargo fmt --check` and `cargo clippy -- -D warnings`, both clean |

Every one of the 355 Rust tests is **new** coverage: the Python package has no
library-level unit tests, so nothing was translated — the unit tests were
written from the Python source and validated against live values.

## Defects found by validation

Six defects were found by comparing against a running oracle that reading the
source did not reveal. They are documented in full in `07_translation_log.md`;
in summary:

| # | Defect | Found by |
|---|---|---|
| 1 | `edit_distance` raises under NumPy 2 (stale generated `.c` references the removed `np.Inf`) | calling it |
| 2 | `CellClusterer.__call__` is unreachable — three independent upstream bugs | calling it |
| 3 | `noodles` reports MAPQ 255 as `None` while `pysam` returns 255 — `--mapping-quality=20` discarded an entire file | `sam_parity.rs` |
| 4 | The port's log-level comparison was inverted, silently swallowing every error and warning | running the upstream suite and getting no diagnosis |
| 5 | `Utilities.Start` seeds the stdlib RNG — the Step-5 analysis concluded it did not, and recorded a deviation that was wrong | comparing three Python runs |
| 6 | Python sorts its output even when writing to a pipe | sorted-vs-unsorted diff |

Defects 3, 4 and 6 were bugs **in the port**. Defects 1 and 2 are upstream.
Defect 5 was an error in the port's own analysis, now corrected — deviation
D11 is withdrawn in `05_design.md` rather than deleted, because the mistaken
reading shaped the first implementation.

## Intentional deviations

Enumerated with reasoning in `05_design.md` § Intentional deviations. Summary:

| # | Deviation | Observable? |
|---|---|---|
| D1 | `edit_distance` unequal lengths returns `None` (the `.pyx` intent) rather than raising | only on input the shipped Python build rejects |
| D2 | An unknown cluster method is rejected at construction, not first use | no; the error type is the same |
| D3 | `CellClusterer::call` returns `Err(UpstreamDefect)` | no; the Python path always raises |
| D4 | `Record`'s discarded `entry_format` parameter is absent | no; it provably had no effect |
| D5 | Whitelist PNGs are not drawn; their data and the `_cell_thresholds.tsv` sidecar are | yes, if a user wanted the plots |
| D6 | FASTQ parse errors return `Err` instead of `sys.exit(1)` from a library | no; the CLI reproduces the message and exit code |
| D7 | `getErrorCorrectMapping` drops a `str()` coercion that would corrupt `bytes` keys | no; in-package callers pass `str` |
| D8 | Provably-dead parameters of `TwoPassPairWriter` are absent | no |
| D9 | The extension installs as `umi_tools_rs`, not `umi_tools._native` | yes, to a Python caller — and deliberately, so both packages import |
| D10 | *(resolved)* CLI help layout | **no longer a deviation** — the help text is captured from the live command and is byte-identical |
| D11 | *(withdrawn)* the tie-break RNG | **no longer a deviation** — CPython's generator is reproduced exactly |
| D12 | Remote reference URLs are not fetched | yes; see Known gaps |
| D13 | CRAM output of sequence-less records is refused rather than crashing | yes; see Known gaps |

## Known gaps

Two, both detected and reported with an actionable message. Together they
account for the 3 upstream cases that do not pass.

### D12 — remote reference sequences

htslib fetches a `--reference-filename` given as an `http(s)`/`ftp` URL; this
port reads local files only. Adding it would mean an HTTP client and a TLS
stack in a crate that otherwise touches no network.
`core::alignment::is_remote_reference` detects the URL and returns an error
naming the workaround rather than surfacing "No such file or directory".

All three CRAM cases in the upstream suite pass a GitHub URL for a reference
that is also checked in locally. **With the local path substituted,
`dedup` from CRAM produces output byte-identical to Python's** — verified
directly, so the gap is in fetching the reference, not in reading CRAM.

### D13 — CRAM output of records with no sequence

A SAM record may carry `*` for its sequence, and the package's own
`chr19_1mb.cram` fixture is entirely such records. `noodles-cram`'s writer
assumes a sequence is present when converting the CIGAR into CRAM features and
**panics** — an upstream panic on ordinary user input. Since a panic on user
input is what this port's error model forbids, the case is detected before the
writer sees it and returns `PortError::Alignment`. CRAM *input* of the same
records works, as does SAM and BAM output.

## Reproducibility

### Data version

No external reference database. Fixtures are committed; their SHA-256 digests
are in `crates/core/tests/fixtures/SHA256SUMS`.

### Code version

- `Cargo.lock` is committed.
- The Python oracle is pinned by `data_exports/pip_freeze.txt`:
  Python 3.11.9, numpy 2.4.6, scipy 1.17.1, pysam 0.24.0, pandas 3.0.5,
  regex 2026.7.19, pybktree 1.1.
- Rust: 1.97.1 stable, `aarch64-apple-darwin`.

### Numeric version

The reference values were produced by NumPy 2.4.6 and SciPy 1.17.1 on
`aarch64-apple-darwin`. The port uses **no BLAS and no parallel reduction** in
any float-sensitive path, so summation order is fixed by construction rather
than by a linked library — which is why the density-knee comparison lands
within 1e-9 rather than needing a Tier-2 correlation threshold.

Both random-number streams are reproduced exactly rather than approximated:
CPython's `random.random()` (`init_by_array` seeding) for the read tie-break,
and NumPy's legacy `RandomState` (`init_genrand` seeding) plus
`numpy.random.choice`'s inverse-transform sampling for the null distribution.
Both are verified against the live interpreters in `core::mersenne`'s and
`core::umi_methods`' unit tests.

## Behavioural notes that are not deviations

- **The output sort happens as records are written.** Python writes an
  unsorted temporary file and then runs `samtools sort` on it;
  `alignment::CoordinateSorter` sorts the records as they are produced, with
  samtools' order (reference, position, strand, then input order), spilling
  sorted runs of 256 MiB to `--temp-dir` and merging them. The output is the
  same; the unsorted temporary file is never written.
- **Log timestamps are UTC.** Python renders them in local time. The format is
  identical and the test harness strips these lines; taking a timezone
  dependency for a cosmetic log field was not worth it.
- **Output headers are the input's text.** Python writes with the input as
  `template`, which keeps each header line and its field order; SAM and BAM
  output here is written from the input's header text too (sorted output with
  `@HD` updated as `samtools sort` updates it), not re-serialised by `noodles`.
  CRAM input and output still go through `noodles`' header.
- **`--gene-transcript-map` transcript order.** Python keeps each gene's
  transcripts in a `set`, so the order it reads them in — and with it which
  reads `dedup` keeps and how `group` numbers its groups — changes with the
  interpreter's hash seed from run to run. The port reads them in the map
  file's order, one of Python's orders; outputs were compared with Python
  run under a hash seed that gives that order. Like Python's `fetch(contig)`,
  this needs a `.bai` or `.csi` index; `--chrom` uses one when present.
- **The input is streamed, never held in memory.** Like Python, `dedup`
  re-opens the input for the passes beside the main one — the multimapping
  tag check, the `--output-stats` UMI frequency table, and
  `TwoPassPairWriter`'s mate searches. The mate search per contig is a forward
  scan of a second reader rather than an index lookup, which visits the same
  records in the same order for a coordinate-sorted input.

## Release readiness

- `Cargo.toml` metadata is complete for `cargo publish`: `description`,
  `license` (MIT, as the original), `repository`, `keywords`, `categories`,
  `readme`.
- `crates/core` sets `exclude = ["tests/fixtures/**"]`, keeping the published
  `.crate` well inside crates.io's 10 MiB cap while the fixtures stay in the
  repository for `cargo test`.
- `maturin build --release` produces an installable wheel.
- `cargo doc --no-deps --workspace` is clean; every `pub` item has a doc
  comment, and fallible ones have an `# Errors` section.
