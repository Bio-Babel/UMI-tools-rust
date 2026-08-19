# 03 — Data Asset Strategy (Step 4)

## Summary

`umi_tools` ships **no bundled runtime data and no reference databases**, and
downloads nothing at runtime. All data in the source tree is test fixtures.

Consequence: `data_assets.yaml` stays `assets: []`, and Step 10 (data staging /
Zenodo) is **skipped**. The remaining work in this step is classifying the test
fixtures and deciding which the Rust crate carries.

## Evidence for "no runtime data"

| Check | Result |
|---|---|
| `MANIFEST.in` | ships only `umi_tools/*.pyx`, `umi_tools/*.c`, `requirements.txt` |
| `[tool.setuptools]` package-data | not declared |
| `importlib.resources` / `pkg_resources` use | none |
| `__file__`-relative data loads | none |
| `requests` / `urllib` / any HTTP client | none imported anywhere in the package |
| `data/` style directory in the package | none |
| Pickle (`.pkl`) fixtures | none |
| `.npy` / `.npz` fixtures | none |

Every file the package reads at runtime is supplied by the user on the command
line: a BAM/SAM/CRAM alignment, a FASTQ, a whitelist TSV, a gene–transcript map
TSV, or a CRAM reference FASTA.

## Classification of the source tree's data

| Asset group | Files | Size | Role | Distribution | Rust target |
|---|---|---:|---|---|---|
| Alignment inputs (`chr19*.bam/.cram/.sam` + indices) | ~20 | ~2.5 MB | `test_fixture` | bundled in-repo | copied into `umi_tools-rs/tests/fixtures/`, read by `noodles` directly |
| FASTQ inputs (`slim.fastq.gz`, `scrb_seq_fastq.*.gz`, `indrop.fastq.*.gz`, …) | ~14 | ~1.6 MB | `test_fixture` | bundled in-repo | `tests/fixtures/`, read by `flate2` + the ported FASTQ iterator |
| Whitelist / barcode TSVs (`scrb_seq_barcodes`, `indrop_fuzzy_whitelist.tsv`, `indrop_blacklist.tsv`, `chr19_gene_assigned*.tsv`) | ~8 | ~1 MB | `test_fixture` | bundled in-repo | `tests/fixtures/` |
| CRAM reference (`chr19_1mb.fa` + `.fai`) | 2 | ~1 MB | `test_fixture` | bundled in-repo | `tests/fixtures/` |
| Expected-output references (`group_*_py3.sam/.tsv`, `processed_*.fastq`, `count_*.tsv`, `*_help`) | ~75 | ~85 MB | `test_fixture` (oracle output) | **not copied** | superseded by live differential comparison against the Python oracle |
| Python-side reference dumps produced during this port | — | small | `test_fixture` | generated | `port_reports/umi_tools/data_exports/*.json` → `crates/core/tests/fixtures/` |

### Why the ~85 MB of expected-output references are not copied

Those files are md5-compared snapshots of the Python CLI's output. Carrying
them would (a) blow past crates.io's 10 MiB package cap by an order of
magnitude, and (b) freeze a *snapshot* as the oracle when a *live* oracle is
available in the same process via PyO3. The port compares against the running
Python package instead (Step 9), which is strictly stronger. The reduced inputs
are kept because `cargo test` must not require the Python environment or network
access.

Total fixture footprint committed to the crate: **~6 MB**, all inputs, no
references. Still above the 10 MiB cap? No — but to keep `cargo publish` lean,
`Cargo.toml` sets `exclude = ["tests/fixtures/**"]`, so fixtures live in the
repository for `cargo test` and are omitted from the published `.crate`.

## Exporter / builder / loader contract

Because there is no runtime data asset, the three-part contract collapses to its
test-fixture form:

1. **Python-side exporter** — `port_reports/umi_tools/data_exports/export_*.py`
   (written in Slice 0). Each runs a specific Python function on a fixed set of
   inputs and dumps `{args, result}` records as JSON. These serve double duty as
   the Step-8 Tier-1 reference values and as `cargo test` fixtures, which is
   exactly the contract's intent.
2. **Rust-side builder** — none required; the JSON is consumed directly by
   `serde_json` in `#[cfg(test)]` code. No build-time codegen, so no `build.rs`.
3. **Loader API** — for fixtures, `core::testing::load_fixture(name)`
   (`#[cfg(test)]`-gated). For user-supplied data at runtime, the loader *is* the
   public API: `open_input_alignments`, `open_file`, `get_user_defined_barcodes`.

No `include_bytes!`, no `rust-embed`, no `registry.rs`, no `download.rs`, no
`umi_tools_data/` staging directory — none of the machinery those exist for is
present in this package. Creating them would be scaffolding with nothing behind
it.

## Reproducibility axes

1. **Data version** — no versioned external database. Fixtures are pinned by
   being committed; their SHA-256 digests are recorded in
   `crates/core/tests/fixtures/SHA256SUMS` so drift is detectable.
2. **Code version** — `Cargo.lock` committed; Python oracle pinned via
   `data_exports/pip_freeze.txt` (see `00_baseline.md`).
3. **Numeric version** — the only floating-point-sensitive paths are the two
   knee estimators and the percentile clustering method. Reference numbers were
   produced by NumPy 2.4.6 / SciPy 1.17.1 on `aarch64-apple-darwin`; the Rust
   side uses no BLAS and no parallel reduction in those functions, so summation
   order is fixed by construction. Recorded in `08_validation.md`.

## `data_assets.yaml`

Left as `assets: []`. Step 10 is skipped, per the decision rule
"skip if no remote assets".
