# 01 — Package Audit & Dependency Trace (Steps 2 & 6)

## Package audit

### Identity

| Item | Value |
|---|---|
| Python import name | `umi_tools` |
| Distribution name | `umi_tools` (PyPI: `umi_tools`) |
| Version | 1.1.6 |
| Source | `$UMI_TOOLS_SRC` |
| Upstream | https://github.com/CGATOxford/UMI-tools |
| License | MIT |
| Purpose | Deduplicating, grouping and counting sequencing reads by Unique Molecular Identifier (UMI), and extracting/whitelisting cell barcodes |

Machine-generated audit dump: `_audit_raw.md` (from `python_package_audit.py`).

### `__all__` and the public surface

**`__all__` is not declared** in any module. `__init__.py` re-exports exactly one
name:

```python
try:
    from umi_tools.network import UMIClusterer
except ImportError:
    pass
```

The `try/except ImportError` is deliberate (it lets `setup.py` import the
package before dependencies exist), so even `UMIClusterer` is not guaranteed to
be bound at import time. This is itself a translation risk: the package's own
authors have not declared what is public. Per the workflow's fallback rule the
public surface is taken as *every top-level non-underscore name in every module*
— 69 functions and 14 classes — and the whole of it is in scope. `_dedup_umi`
(underscore-prefixed, but the algorithmic core) is also in scope because three
modules import from it.

The real user-facing surface is the **CLI**: `umi_tools <tool>` dispatching to
seven subcommands via `importlib.import_module("umi_tools." + command)`.

### Module inventory

| Module | LOC | Role |
|---|---:|---|
| `Utilities.py` | 1775 | `optparse` machinery, logging, temp files, pysam open/sort helpers, option validation |
| `Documentation.py` | 532 | pure docstring constants concatenated into `__doc__` |
| `sam_methods.py` | 727 | BAM read → (umi, cell) extraction; `get_bundles` position/gene bundler; paired-end writer; cigar/splice logic |
| `whitelist_methods.py` | 718 | knee estimation (distance & density), BK-tree error correction, whitelist file parsing, above-knee error detection |
| `extract_methods.py` | 615 | barcode extraction (string pattern & regex), quality filtering/masking, read rewriting |
| `network.py` | 544 | **the algorithmic core** — `UMIClusterer` (5 clustering methods), `ReadDeduplicator`, `CellClusterer`, substring-index nearest-neighbour search, BFS |
| `dedup.py` | 465 | `dedup` CLI |
| `whitelist.py` | 531 | `whitelist` CLI |
| `extract.py` | 535 | `extract` CLI |
| `group.py` | 326 | `group` CLI |
| `prepare_for_em.py` | 285 | `prepare-for-em` CLI (RSEM/Salmon mate fixing) |
| `count.py` | 233 | `count` CLI |
| `umi_methods.py` | 239 | FASTQ `Record` + iterators, average UMI distance, random UMI generator |
| `count_tab.py` | 172 | `count_tab` CLI |
| `umi_tools.py` | 70 | top-level CLI dispatcher |
| `version.py` | 1 | `__version__ = "1.1.6"` |
| `_dedup_umi.pyx` | 19 | Cython `edit_distance` (Hamming) |

### Native code

- `umi_tools/_dedup_umi.pyx` — 19 lines of Cython. Declared in `pyproject.toml`
  as `ext-modules = [{name = "umi_tools._dedup_umi", sources = ["umi_tools/_dedup_umi.c"]}]`
  — note it compiles the **checked-in generated `.c`**, not the `.pyx`.
- Compiled artefact present: `_dedup_umi.cpython-311-darwin.so`.
- Imported by `network.py`, `umi_methods.py`, `whitelist_methods.py` (all
  importing exactly one symbol, `edit_distance`).
- No `cffi`, `ctypes`, or `pybind11`.

This is the entire native surface: a Hamming distance over two equal-length byte
strings. It is trivially portable and, per the workflow, is the *easiest* piece —
already statically typed. See `04_compiled_code_strategy.md`.

### Dynamic-behaviour risk

The audit script reports 0 markers (no metaclasses, no `__getattr__`, no
`exec`/`eval`, no signature-rewriting decorators, no `*args`/`**kwargs` in public
functions). That number understates the real risk, which is concentrated in four
patterns that only surfaced by running the code:

1. **Method-pointer dispatch as configuration.** `UMIClusterer.__init__`,
   `CellClusterer.__init__` and `ExtractFilterAndUpdate.__init__` assign *bound
   methods* to instance attributes (`self.get_adj_list = self._get_adj_list_directional`).
   The set of assigned attributes differs per method string, and an unrecognised
   string assigns nothing — so `UMIClusterer("bogus")` constructs fine and fails
   only at call time:
   `AttributeError: 'UMIClusterer' object has no attribute 'get_adj_list'`
   (verified at runtime). In Rust this becomes an enum resolved in the
   constructor — a deliberate deviation, recorded in `05_design.md`.
2. **`optparse.Values` as a duck-typed god-object.** `get_bundles`,
   `ReadDeduplicator` and every `main()` take an `options` object and read ~40
   attributes off it. There is no class, no schema, no type. Resolving this into
   concrete option structs is the single largest translation task.
3. **`bytes` vs `str` inconsistency.** `sam_methods` barcode getters return
   `bytes` (`.encode('utf-8')`); `extract_methods` works in `str`;
   `whitelist_methods.getErrorCorrectMapping` does `set([str(x) for x in whitelist])`
   — which stringifies bytes to `"b'ACGT'"` if handed bytes. `ReadDeduplicator`
   calls `cluster[0].decode()`. `count_tab.main` does `try: gene.decode() except AttributeError: pass`.
   The Rust port normalises on `Vec<u8>`/`&[u8]` for barcodes and documents each
   boundary.
4. **Upstream defects reachable from the public surface** (all verified by
   running, not reading):
   - `CellClusterer.__call__` is **dead code**: it calls
     `self.get_connected_components(umis, adj_list, counts)` (3 args) against a
     method whose signature is `(self, graph, counts)`, so it always raises
     `TypeError: CellClusterer._get_connected_components_adjacency() takes 3 positional arguments but 4 were given`.
     It would then call `self.get_groups`, which `CellClusterer.__init__` never
     assigns, and `_get_connected_components_adjacency` calls `.sort()` on the
     `set` returned by `breadth_first_search`. Three independent bugs. Nothing in
     the package calls it.
   - `edit_distance` on unequal-length inputs raises `AttributeError` under
     NumPy ≥ 2 (`np.Inf` in the generated `.c`). See `00_baseline.md`.
   - `errorDetectAboveThreshold` contains `error_counter["error_discarded_mt_1"]`
     — a bare subscript with no `+= 1`. The counter is only *read* (defaults to
     0), so the "discarded due to errors from multiple other CBs" log line always
     prints 0. Faithfully mirrored, with a note.
   - `ExtractFilterAndUpdate.__call__` computes `new_identifier2` in the
     `either_read` branch and then assigns `read2.identifier = new_identifier`
     (not `…2`). Both are built from the same `umi`/`cell` here so the observable
     result is identical; mirrored as written.

### Data assets

- **No bundled runtime data.** `MANIFEST.in` ships only `*.pyx`, `*.c` and
  `requirements.txt`. There is no `package-data`, no `importlib.resources` use,
  no `__file__`-relative data loading, and no `data/` directory.
- **No remote downloads.** No `requests`, no `urllib` anywhere in the package.
- Test fixtures live in `tests/` (~121 files: BAM/CRAM/FASTQ/TSV + references)
  and are *not* installed with the package.

Full classification in `03_data_asset_strategy.md`.

### Existing tests

`tests/test_umi_tools.py` drives `tests/tests.yaml` (65 named cases; 70
parameterised tests after CRAM/sort variants). Every case is an end-to-end CLI
invocation compared against a checked-in reference file. Distribution by
subcommand: dedup 19, group 16, extract 15, whitelist 10, count 4, count_tab 3.
`tests/test_style.py` is a pycodestyle check.

**Baseline run: 70 passed in 129.55 s** (from repo root, oracle env).

There are **no library-level unit tests**. Implication for Step 8 recorded in
`00_baseline.md` and `08_validation.md`.

---

## Dependency trace (Step 6)

Only paths that actually affect `umi_tools` semantics are traced. Upstream
dependencies are not bulk-ported.

### Symbol-level trace

| Python dependency | Symbols actually used | Where | Rust decision | Rationale |
|---|---|---|---|---|
| `numpy` | `np.median`, `np.log10`, `np.linspace`, `np.power`, `np.vstack`, `np.sum`, `np.sqrt`, `np.outer`, `np.argmax`, `np.cumsum`, `np.digitize`, `np.bincount`, `np.random.choice`, `np.random.seed`, `np.inf` | `network`, `whitelist_methods`, `dedup`, `umi_methods` | **reimplement locally**, in `crates/core/src/numeric.rs` | Every use is a 1-D reduction or elementwise op over `Vec<f64>`/`Vec<i64>`. Pulling `ndarray` in for `median`/`cumsum`/`argmax` buys nothing and adds a semantic-divergence surface. Median tie-handling, `digitize(right=True)` and `bincount(minlength=)` are reimplemented to NumPy's exact definition and differentially tested. |
| `numpy.matlib` | `npm.repmat` | `whitelist_methods.getKneeEstimateDistance` | reimplement (it is a broadcast of a 2-vector across `nPoints` rows) | `numpy.matlib` is deprecated upstream; the operation is a dot product in disguise. |
| `scipy.stats` | `gaussian_kde(x, bw_method=0.1)` and evaluation at 10 000 points | `getKneeEstimateDensity` | **reimplement locally** — Gaussian KDE with Scott/`bw_method` scalar factor | 1-D KDE with a scalar bandwidth factor is ~30 lines. The only candidate crate (`criterion`-adjacent stats crates) does not reproduce SciPy's covariance-factor convention. Fidelity requires matching `covariance = cov(dataset, bias=False) * factor**2`; that is specified exactly and Tier-2 validated. |
| `scipy.signal` | `argrelextrema(y, np.less)` | `getKneeEstimateDensity` | reimplement (strict-less-than comparison against both neighbours, `order=1`, `mode='clip'`) | 10 lines; SciPy's exact boundary semantics are mirrored. |
| `pandas` | `DataFrame`, `groupby().agg({'counts': ['median', len, 'sum']})`, `merge(how='left')`, `fillna(0).astype(int)`, `to_csv(sep='\t')` | `dedup.py` `--output-stats` only | **reimplement locally** with `BTreeMap` aggregation | Only one code path uses pandas, over a two-column frame. `polars` would be a heavyweight dependency for a group-by-median-len-sum. The exact column order, sort order and int cast are mirrored. |
| `pysam` | `AlignmentFile`, `.fetch`, `.references`, `AlignedSegment` accessors (`is_read2`, `is_unmapped`, `mapq`, `cigar`, `cigarstring`, `pos`, `aend`, `tlen`, `get_tag`, `set_tag`, `to_dict`/`from_dict`), `pysam.sort` | `sam_methods`, `dedup`, `group`, `count`, `prepare_for_em`, `Utilities` | **adopt `noodles`** (`noodles-sam`, `noodles-bam`, `noodles-cram`, `noodles-bgzf`) | `noodles` is the maintained pure-Rust htslib-equivalent stack and covers SAM/BAM/CRAM read+write, which is exactly the surface used. `pysam.sort` (a samtools shim) is replaced by an in-crate coordinate sort. |
| `regex` (PyPI, fuzzy) | `regex.compile("(%s){e<=%i}")`, `{s<=%i}`, `{i<=1}`, `.match`, `.groupdict`, `.span`, `.groupindex` | `whitelist_methods.checkError`, `errorDetectAboveThreshold`, `CellClusterer` fuzzy path, `extract_methods` regex path | **hybrid**: `regex` crate for the ordinary named-group patterns; a purpose-built fuzzy matcher in `crates/core/src/fuzzy.rs` for `{e<=n}`/`{s<=n}`/`{i<=n}` | The Rust `regex` crate has no approximate-matching extension, and no crates.io crate implements Python-`regex`'s fuzzy syntax. The three fuzzy uses are all "does this whole barcode match this whole barcode within *n* errors", i.e. bounded edit / substitution / insertion distance — implementable exactly and much easier to validate than a general fuzzy engine. Documented as an input-form narrowing in `05_design.md`. |
| `pybktree` | `BKTree(hamming_distance, whitelist)`, `.find(item, threshold)` | `getErrorCorrectMapping` | **reimplement locally** (`crates/core/src/bktree.rs`) | ~60 lines; a BK-tree over a caller-supplied metric. `bk-tree` on crates.io exists but its `Metric` trait and tie-ordering differ; local control makes candidate *ordering* reproducible, which matters because the caller branches on `len(candidates)`. |
| `matplotlib` | `pyplot`, `lines.Line2D`, `figure`, `savefig` | `whitelist_methods` plotting, `whitelist --plot-prefix` | **not ported as plotting**; the *data* behind every plot is emitted instead | Plot rendering is not library semantics and has no Rust equivalent with matching output bytes. The `_cell_thresholds.tsv` sidecar file (which *is* semantic and is compared by the test suite) is fully ported. Recorded as an intentional deviation. |
| `future` / `builtins` / `six` | `iteritems`, `dict` | several | **drop** | Python-2 compatibility shims with no meaning in Rust. |
| `optparse` | `OptionParser`, `Option`, `OptionGroup`, `IndentedHelpFormatter`, `SUPPRESS_HELP` | `Utilities`, every `main()` | **adopt `clap`** (builder API) + a bespoke help formatter | `optparse` is deprecated in Python itself. `clap`'s builder API reproduces the option set; the `BetterFormatter`/`AppendCommaOption` behaviours (comma-splitting `append` options, `""`/`"None"` treated as defaults, indented wrapped help) are reimplemented on top. Deviation: help *text layout* is not byte-identical. |
| `logging` | `Formatter`, handlers, levels | `Utilities` | reimplement a minimal logger writing the same line format | `MultiLineFormatter` prefixes continuation lines with `#` + indent; that exact format is part of the log files the test suite compares. |
| `gzip` | `gzip.open` (rt/wt) | `Utilities.openFile` | **adopt `flate2`** | Transparent `.gz`/`.z` handling on read and write, with `compresslevel` honoured. |
| `tempfile`, `os`, `time`, `uuid`, `inspect`, `textwrap`, `copy`, `itertools`, `collections`, `random`, `re` | assorted | everywhere | std / `tempfile` / `uuid` / `rand` crates | Direct equivalents. `itertools.combinations`/`product` → `itertools` crate. `collections.Counter`/`defaultdict` → `HashMap`/`BTreeMap` (ordering decisions in `05_design.md`). |

### Algorithm backend selection

For each nontrivial algorithm: Python entry point → actual backend → Rust choice.

| Algorithm | Python entry point | Actual backend | Rust implementation | Risk |
|---|---|---|---|---|
| Hamming distance | `_dedup_umi.edit_distance` | hand-written C loop over `char*` | `core::distance::edit_distance` — byte loop, `Option<u32>` (`None` = unequal length) | none; Tier 1 exact |
| Directional/adjacency UMI networks | `UMIClusterer._get_adj_list_*` | pure Python + `itertools.combinations` | `core::network` | tie-ordering (see below) |
| Approximate nearest neighbours | `build_substr_idx` / `iter_nearest_neighbours` | pure Python `defaultdict(set)` | `core::network::SubstrIdx` | **set-iteration order** — see below |
| Connected components | `breadth_first_search` | pure Python `set`-based BFS | `core::network::breadth_first_search` → `BTreeSet` | component contents are sorted before use, so order-independent |
| Density knee | `getKneeEstimateDensity` | `scipy.stats.gaussian_kde` + `scipy.signal.argrelextrema` | `core::knee::density` | Tier 2 — float divergence in the KDE sum |
| Distance knee | `getKneeEstimateDistance` | NumPy vector ops | `core::knee::distance` | Tier 1/2 — pure float arithmetic, same op order |
| BK-tree error correction | `getErrorCorrectMapping` | `pybktree` | `core::bktree` | candidate order (branch is on `len`, so order-insensitive) |
| Fuzzy barcode match | `checkError`, `errorDetectAboveThreshold` | `regex` module fuzzy | `core::fuzzy` | narrowed input form; documented |
| Read bundling | `sam_methods.get_bundles` | pysam iteration + `random.random()` tie-break | `core::bundles` + `noodles` | **RNG stream** — see below |
| Per-UMI stats aggregation | `dedup.aggregateStatsDF` | pandas groupby | `core::stats` | median convention mirrored |

### Known numeric / determinism risks

1. **Python `set` iteration order.** `iter_nearest_neighbours` yields pairs from
   a `set` union; `_group_adjacency` does `list(connected_nodes - observed)`.
   Where the downstream consumer sorts (BFS components, `_group_cluster`,
   `_group_directional`) the order washes out. Where it does not
   (`_group_adjacency`'s trailing members, `remove_umis` cardinality) the Python
   output itself is not reproducible across runs with `PYTHONHASHSEED` unset for
   `str` keys — though UMI keys are `bytes`, whose hash *is* stable within a
   build. The Rust port uses `BTreeSet`/sorted iteration to make the order
   *deterministic*, and differential tests compare as **sets** where Python's own
   order is an artefact. Each such site is listed in `07_translation_log.md`.
2. **`random.random()` in `get_bundles.update_dicts`.** Reservoir sampling picks
   the retained read among equal-MAPQ, equal-multimapping candidates. Python
   seeds this via `np.random.seed(options.random_seed)` — but `update_dicts` uses
   the **stdlib `random`** module, which `np.random.seed` does *not* affect. So
   the read chosen among ties is genuinely non-reproducible in Python. Tier 3;
   validated on counts and on the set of retained UMIs, not on which specific
   read object survives.
3. **`np.random.choice` in `random_read_generator`.** Only feeds
   `dedup --output-stats` null distributions. Tier 3 (summary statistics).
4. **Float summation order** in the KDE and knee-distance paths. Same order is
   preserved deliberately (no rayon parallel reduction in these functions) to
   keep the divergence at Tier 1/2 rather than Tier 3.
