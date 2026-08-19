# 00 — Executable Baseline (Step 1)

Status: **PASS** (both environments).

## Python Environment (the oracle)

User-provided, reused as instructed — no new environment created.

| Item | Value |
|---|---|
| Interpreter | `$PY_ENV/bin/python` |
| Python version | 3.11.9 (venv) |
| `umi_tools.__file__` | `$UMI_TOOLS_SRC/umi_tools/__init__.py` |
| Install mode | editable (`pip install -e .`) — points at the live source tree, so edits are visible to differential tests |
| `umi_tools.__version__` | `1.1.6` |
| CLI entry point | `$PY_ENV/bin/umi_tools` |

### Key dependency versions

| Package | Version | Note |
|---|---|---|
| numpy | 2.4.6 | NumPy 2.x — see the `np.Inf` finding below |
| pandas | 3.0.5 | only used by `dedup --output-stats` |
| pysam | 0.24.0 | BAM/SAM/CRAM I/O + `samtools sort` shim |
| scipy | 1.17.1 | `gaussian_kde`, `argrelextrema` (density knee) |
| matplotlib | 3.11.1 | whitelist diagnostic plots only |
| regex | 2026.7.19 | fuzzy-match regex (`{e<=1}`, `{s<=1}`, `{i<=1}`) |
| pybktree | 1.1 | BK-tree used by `getErrorCorrectMapping` |
| maturin | 1.14.1 | present in the same env — PyO3 bindings buildable |

Full `pip freeze`: `data_exports/pip_freeze.txt`.

### Smoke test

```python
>>> from umi_tools.network import UMIClusterer
>>> UMIClusterer(cluster_method='directional')({b'ACGT':10, b'ACGA':2, b'TTTT':5}, threshold=1)
[[b'ACGT', b'ACGA'], [b'TTTT']]
```

Passes. `step1_unit_test.py` in the repo root is *stale*: it calls
`clusterer(input_data.keys(), input_data, threshold=1)` (3 args) while the
1.1.6 `UMIClusterer.__call__` signature is `(umis, threshold)` (2 args). This
is the first concrete instance of the "read source, not names" rule — the
in-repo test is not a valid oracle and was not used as one.

### Existing Python test suite

`tests/test_umi_tools.py` + `tests/tests.yaml` — a data-driven **end-to-end CLI**
harness (~62 parameterised cases). It shells out to the `umi_tools` binary,
compares stdout/output files against checked-in references by md5 (with
per-format decoding for CRAM). It contains **zero** unit tests of the library
API. Consequence for Step 8: `cargo test` coverage parity cannot be achieved by
mechanical `pytest`→`#[test]` translation; unit coverage must be authored from
the Python source and validated differentially. Recorded again in
`08_validation.md`.

Must be run from the repo root (`pytest tests/test_umi_tools.py`), not from
`tests/` — the harness resolves `tests/tests.yaml` relative to `os.getcwd()`.

Baseline run result: see `08_validation.md` § Python baseline test run.

### Runtime finding — `_dedup_umi.edit_distance` is broken for unequal lengths in this env

`umi_tools/_dedup_umi.pyx` returns `np.inf` when `len(a) != len(b)`. The
**checked-in, pre-generated** `_dedup_umi.c` (the file `pyproject.toml`
actually compiles) emits `np.Inf`, which NumPy 2.0 removed. So:

```python
>>> edit_distance(b"ACGT", b"ACG")
AttributeError: `np.Inf` was removed in the NumPy 2.0 release. Use `np.inf` instead.
```

Also, the function takes `char *`, so `str` arguments raise
`TypeError: expected bytes, str found`.

The *intended* semantics (from the `.pyx` gold standard) are: unequal lengths →
infinite distance, i.e. never `<= threshold`. The Rust port implements the
intent, not the NumPy-2 accident. Recorded as an intentional deviation in
`05_design.md`.

## Rust Toolchain

User environment; nothing installed by this port.

| Item | Value |
|---|---|
| `rustc` | 1.97.1 (8bab26f4f 2026-07-14) |
| `cargo` | 1.97.1 (c980f4866 2026-06-30) |
| host triple | `aarch64-apple-darwin` |
| channel | stable (no nightly features required) |
| `maturin` | 1.14.1 (in the Python env above) |

### Gate: does the toolchain actually link?

```
$ cargo new --lib _probe && cd _probe && cargo build
   Compiling _probe v0.1.0
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 0.41s
```

Passes. Probe crate deleted.

No native system libraries are needed: the port has no BLAS/LAPACK/HDF5
dependency (see `01_audit.md` § Dependency trace — the only linear algebra is
elementwise vector arithmetic in the distance-knee estimator).

## Decision

Both gates pass. Proceed to Step 2.
