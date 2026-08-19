# 04 — Native Code & Dynamic-Behaviour Strategy (Step 5)

Two independent inputs feed the staged-architecture decision. For this package
they point in opposite directions: the native code is trivial, the dynamic
Python is where the work is.

## Part 1 — Existing native code

### Inventory

| File | Lines | Content |
|---|---:|---|
| `umi_tools/_dedup_umi.pyx` | 19 | one `cpdef int edit_distance(a, b)` |
| `umi_tools/_dedup_umi.c` | 251 407 | Cython-generated from the above (checked in) |
| `umi_tools/_dedup_umi.cpython-311-darwin.so` | 85 KB | compiled artefact |

No `cffi`, no `ctypes`, no `pybind11`, no `-sys` crate candidate, no vendored
third-party C library. `pyproject.toml` compiles the checked-in `.c`, not the
`.pyx`, which is why the generated file's stale `np.Inf` reference survives.

### The entire native surface

```cython
cpdef int edit_distance(a, b):
    cdef char * aa = a
    cdef char * bb = b
    cdef int k, l, c
    c = 0
    la = len(a)
    lb = len(b)
    if la != lb:
        return np.inf
    for k from 0 <= k < la:
        if aa[k] != bb[k]: c += 1
    return c
```

Hamming distance over two equal-length byte strings, with a sentinel for
unequal lengths.

### Coupling assessment

- **Coupling to CPython internals**: only the `char *` coercion from `bytes`.
- **Coupling to NumPy**: only the `np.inf` sentinel — and the `cpdef int`
  return declaration contradicts it (an `int`-typed function returning a float
  infinity), which is why the return type in Rust has to be a decision, not a
  transcription.
- **Extractability**: total. Three callers, each importing only this symbol.

### Strategy: **reimplement in safe Rust** (not re-host, not wrap)

Re-hosting the `.so` through `bindgen`/a `-sys` crate would be strictly worse:
it would drag in a CPython dependency for 12 lines of arithmetic and would
preserve the NumPy-2 defect. There is no performance argument either — the Rust
byte loop is at least as fast as the C one, and the port can additionally give
the compiler a slice pair to autovectorise.

Signature resolution:

| Python | Rust |
|---|---|
| `edit_distance(a: bytes, b: bytes) -> int \| float('inf')` | `pub fn edit_distance(a: &[u8], b: &[u8]) -> Option<u32>` |

`None` means "unequal lengths / infinite distance". Every call site in the
Python source is of the form `edit_distance(x, y) <= threshold`, so `None` maps
to `false` — captured once in `core::distance::within_threshold(a, b, t) -> bool`
so no caller re-derives it. A `edit_distance_f64` convenience returning
`f64::INFINITY` is provided for the `get_average_umi_distance` path, which
averages the raw values.

Deviation recorded: Python (as built here, NumPy ≥ 2) *raises* on unequal
lengths; the `.pyx` gold standard *returns infinity*. The port implements the
`.pyx` intent. `05_design.md` § Intentional deviations.

## Part 2 — Pure-Python dynamic-behaviour risk

This is the hard part, and the audit script's "0 markers" score does not see it.
Four patterns, each with a required Step-7 resolution.

### 2.1 Method-pointer dispatch as configuration (**high**)

`UMIClusterer.__init__`, `CellClusterer.__init__` and
`ExtractFilterAndUpdate.__init__` do not store a mode — they rebind instance
attributes to different bound methods:

```python
if cluster_method == "adjacency":
    self.get_adj_list = self._get_adj_list_adjacency
    self.get_connected_components = self._get_connected_components_adjacency
    self.get_groups = self._group_adjacency
elif ...
```

Three consequences that only appear at runtime:

- An unrecognised method string binds **nothing**, so construction succeeds and
  the first call fails with `AttributeError: 'UMIClusterer' object has no
  attribute 'get_adj_list'`. Verified.
- The `null` variants (`_get_adj_list_null`, `_get_connected_components_null`)
  return `None` / the input unchanged, so the "adjacency list" parameter is
  polymorphic across `dict[bytes, list[bytes]]` and `None`, and the "clusters"
  parameter across `list[list[bytes]]` and `list[bytes]`. `_group_percentile`
  exploits this: it receives the flat UMI list and calls
  `_get_best_percentile(clusters, counts)` on it.
- `ExtractFilterAndUpdate` additionally rebinds `self.extract` / `self.joiner`
  for 5′ vs 3′, and leaves `self.pattern_length` etc. **unset** when no pattern
  is given.

**Resolution**: `enum ClusterMethod { Unique, Percentile, Cluster, Adjacency, Directional }`
with `FromStr`; `UmiClusterer::new` returns `Result<_, PortError::UnknownClusterMethod>`,
moving the failure from call time to construction time. The polymorphic
intermediates become `enum AdjList { Adjacency(BTreeMap<..>), Null }` and
`enum Clusters { Components(Vec<Vec<Umi>>), Flat(Vec<Umi>) }`. This is the
single most important type decision in the port; it is what makes the five
clustering methods expressible without `dyn Any`.

### 2.2 `optparse.Values` as an untyped god-object (**high**)

`get_bundles`, `ReadDeduplicator` and all seven `main()` functions take an
`options` object and read attributes off it that were defined across four
different `add_*_options` blocks in `Utilities.Start`. There is no class and no
schema; the set of valid attributes is only discoverable by reading `Start`'s
1000-line body and every `parser.set_defaults` call in every subcommand.

Attributes actually read (extracted from source, not guessed): `method`,
`threshold`, `per_cell`, `per_gene`, `per_contig`, `gene_tag`, `assigned_tag`,
`skip_regex`, `gene_transcript_map`, `get_umi_method`, `umi_sep`, `umi_tag`,
`cell_tag`, `umi_tag_split`, `umi_tag_delim`, `cell_tag_split`,
`cell_tag_delim`, `paired`, `spliced`, `ignore_tlen`, `ignore_umi`,
`read_length`, `soft_clip_threshold`, `subset`, `mapping_quality`,
`detection_method`, `whole_contig`, `chrom`, `unmapped_reads`, `unpaired_reads`,
`chimeric_pairs`, `output_unmapped`, `filter_umi`, `umi_whitelist`,
`umi_whitelist_paired`, `random_seed`, `stdin`, `stdout`, `stdlog`, `tmpdir`,
`in_sam`, `out_sam`, `in_format`, `out_format`, `no_sort_output`,
`input_options`, `output_options`, `reference_filename`, `compresslevel`,
`loglevel`, `timeit_file`, `timeit_name`, `timeit_header`.

**Resolution**: split into three concrete structs with `Default` impls —
`BundleOptions` (what `get_bundles` reads), `DedupOptions` (what
`ReadDeduplicator` reads), `IoOptions` (streams/format/tmpdir) — plus per-
subcommand `clap` structs that build them. The library crate never sees a
`clap` type; the CLI crate does the translation. This keeps `crates/core`
usable as a library, which the Python package's design does not allow.

### 2.3 `bytes`/`str` inconsistency across module boundaries (**medium**)

- `sam_methods.*` barcode getters return `bytes` (they `.encode('utf-8')`).
- `extract_methods.*` operates entirely on `str`.
- `whitelist_methods.getErrorCorrectMapping` does `set([str(x) for x in whitelist])`
  — passing `bytes` there produces the literal string `"b'ACGT'"`.
- `ReadDeduplicator.__call__` calls `cluster[0].decode()` for whitelist lookup.
- `count_tab.main` does `try: gene = gene.decode() except AttributeError: pass`.

**Resolution**: barcodes are `Vec<u8>`/`&[u8]` throughout `core`, because that
is what the BAM path produces and what `edit_distance` needs. Conversion to
`str` happens only at output formatting, via a single
`core::utilities::decode_utf8_lossy_strict` helper that mirrors Python's
`.decode()` (error on invalid UTF-8, matching Python's default `strict`).
`getErrorCorrectMapping` is ported to take `&[Vec<u8>]` — its `str()` call is a
Python idiom for normalising an unknown input type, not intended semantics, and
carrying it forward would mean reproducing `"b'ACGT'"`. Recorded as an
intentional deviation.

### 2.4 Upstream defects reachable from the public surface (**medium**)

All four verified by execution, not reading. Full detail in `01_audit.md`
§ Dynamic-behaviour risk. Per-item disposition:

| Defect | Disposition |
|---|---|
| `CellClusterer.__call__` — 3 independent bugs, always raises `TypeError` | ported as `CellClusterer::call() -> Result<_, PortError::UpstreamDefect>`; the working pieces (`adj_list_directional`, `connected_components_adjacency`, `get_best_min_account`) are ported and usable |
| `edit_distance` unequal-length `AttributeError` under NumPy ≥ 2 | `.pyx` intent implemented (see Part 1) |
| `errorDetectAboveThreshold` — `error_counter["error_discarded_mt_1"]` with no increment | mirrored exactly (the counter stays 0, the log line prints 0) — this one is *observable* in the log output the test suite compares |
| `ExtractFilterAndUpdate.__call__` either-read branch assigns `new_identifier` where `new_identifier2` was computed | mirrored; both values are equal here so behaviour is unchanged |

The rule applied: mirror a defect when it is observable in output the oracle
produces; correct it only when the Python code is *unreachable* (`CellClusterer.__call__`)
or when the defect is an artefact of the build environment rather than the
source (`np.Inf`). Every case is listed in `05_design.md` § Intentional
deviations.

### 2.5 Non-dynamic but worth flagging: recursion and iteration order

- `network.py` calls `sys.setrecursionlimit(10000)` at import and
  `breadth_first_search_recursive` catches `RecursionError` to fall back to the
  iterative BFS. Rust has no catchable stack overflow; the port implements the
  recursive variant with an explicit depth counter and the same fallback, and
  documents the mechanism difference. (`breadth_first_search_recursive` is dead
  code upstream — the caller is commented out — but it is public, so it is in
  scope.)
- `recursive_search` stores state on the **function object**
  (`recursive_search.component = set(...)`). Rust has no function attributes;
  the port threads the accumulator as a `&mut BTreeSet` parameter, with the
  public `recursive_search(node, adj_list)` wrapper preserving the Python
  signature. This is an API-shape adaptation, not a behaviour change.

## Staged-architecture consequence

The native code is a Stage-D candidate on day one (pure Rust, no Python
dependency, trivially validated). The dynamic-Python core is where the Stage-B
differential-testing window earns its keep: `UMIClusterer` and
`ExtractFilterAndUpdate` must be exercised against the live oracle through PyO3
before being declared parity-correct, because their behaviour is defined by
attribute rebinding that no amount of static reading fully pins down.

Per-module stage assignments: `05_design.md` § Staged-architecture module table.
