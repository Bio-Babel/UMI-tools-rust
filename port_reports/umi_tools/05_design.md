# 05 — Design: Translation Map & Staged Architecture (Step 7)

## Crate layout

```
umi_tools-rs/                      workspace root
├── Cargo.toml                     [workspace] members = core, cli, pyo3-bindings
├── Cargo.lock                     committed
├── rust-toolchain.toml            stable + clippy + rustfmt
├── crates/core/                   package `umi_tools-core`, lib `umi_tools_core`
│   └── src/
│       ├── lib.rs                 crate root, `pub use` re-exports = public API
│       ├── error.rs               PortError hierarchy (Python exception mapping)
│       ├── distance.rs            _dedup_umi.pyx
│       ├── network.rs             network.py
│       ├── umi_methods.rs         umi_methods.py
│       ├── extract.rs             extract_methods.py
│       ├── whitelist.rs           whitelist_methods.py
│       ├── sam.rs                 sam_methods.py (pure functions)
│       ├── bundles.rs             sam_methods.get_bundles (split for size)
│       ├── alignment.rs           Utilities.py alignment open/sort helpers
│       ├── utilities.rs           Utilities.py logging/IO/format helpers
│       ├── options.rs             the optparse god-object, made concrete
│       ├── numeric.rs             the NumPy/SciPy subset actually used
│       ├── stats.rs               dedup.aggregateStatsDF (the pandas path)
│       ├── bktree.rs              pybktree
│       ├── fuzzy.rs               the `regex` module's fuzzy subset
│       └── documentation.rs       Documentation.py
├── crates/cli/                    package `umi_tools-cli`, bin `umi_tools`
│   └── src/                       one module per subcommand + dispatch/parser
├── crates/pyo3-bindings/          package `umi_tools_rs`, module umi_tools_rs._native
└── validation/                    differential harnesses (not published)
```

`crates/cli` is an addition to the scaffold: the seven `main()` functions and
the `optparse` machinery are a genuine third concern, and keeping them out of
`core` is what makes `core` usable as a library — something the Python package's
structure does not permit, since every module imports `Utilities`, which owns
global mutable state (`global_options`, `global_benchmark`).

## Ownership model

The Python data flow is dominated by dictionaries of UMIs to counts and to
reads. Decisions, per the "who owns what" rule:

| Python object | Lifetime/sharing in Python | Rust ownership |
|---|---|---|
| UMI / cell barcode (`bytes`) | immutable, copied freely, used as dict key | `Vec<u8>` owned in collections, `&[u8]` in parameters. Never `String` — the BAM path produces arbitrary bytes and `edit_distance` is byte-wise. |
| `counts: dict[bytes, int]` | read-only inside clustering | `&BTreeMap<Vec<u8>, u64>` borrowed. `BTreeMap` not `HashMap`: several call sites iterate and the Python order is an artefact (see § Determinism). |
| `adj_list: dict[bytes, list[bytes]]` | built then read | owned `BTreeMap<Vec<u8>, Vec<Vec<u8>>>` returned by value; borrowed by consumers. |
| `bundle: dict[umi, {"read":…, "count":…}]` | built by `get_bundles`, consumed by `ReadDeduplicator` | `BTreeMap<Vec<u8>, BundleEntry>` moved into the consumer. `BundleEntry` owns its read(s). |
| pysam `AlignedSegment` | GC'd, aliased freely (stored in dicts *and* written out) | `noodles::sam::alignment::RecordBuf` — owned, cloned where Python aliased. Cloning a record buf is a memcpy of a few hundred bytes; `Rc<RefCell<_>>` would buy nothing since reads are effectively immutable after bundling. |
| `options` (`optparse.Values`) | read-only, passed everywhere | `&BundleOptions` / `&DedupOptions` borrowed. Never cloned into structs. |
| `ExtractFilterAndUpdate.read_counts` (Counter) | mutated on every read | owned `BTreeMap<String, u64>` field, mutated through `&mut self`. This is why `process` takes `&mut self` where Python's `__call__` looks pure. |
| `random_read_generator` internal buffers | mutated | owned fields, `&mut self` methods. |
| `TwoPassPairWriter` (holds two open files) | mutated, must be closed | owns both handles; `close(self)` consumes `self`, so use-after-close is a compile error rather than a runtime surprise. |

No `Rc`, no `Arc`, no `RefCell`, no `Mutex` anywhere in the port. The Python
code has no genuine shared-mutable-state graph: it has one god-object read
immutably, plus per-object mutable accumulators. Where `.clone()` appears it is
on `RecordBuf` at a point where Python aliased a read into two places, and it is
noted at the call site — never used as a borrow-checker escape hatch.

## Error model

Python raises `ValueError` almost exclusively, plus a few `AssertionError`s and
one `NotImplementedError`. `U.error()` is *not* an exception — it logs at
CRITICAL and calls `sys.exit(1)`.

```rust
pub enum PortError {
    Value(String),              // ValueError
    Assertion(String),          // AssertionError
    NotImplemented(String),     // NotImplementedError
    Io { path, source },        // IOError/OSError
    Format(String),             // malformed FASTQ/SAM/TSV input
    UnknownClusterMethod(String),
    UpstreamDefect { symbol, detail },
    Alignment(String),          // noodles/samtools-equivalent failures
}
```

Mapping table (`crates/pyo3-bindings` reconstructs the Python type, never a
generic `RuntimeError`):

| Rust variant | Python exception | Example site |
|---|---|---|
| `Value` | `ValueError` | `validate_sam_options`, `get_barcode_read_id`, `knee_estimate_distance`'s "Something's gone wrong here!!" |
| `Assertion` | `AssertionError` | `UmiClusterer::cluster`'s equal-length check, `get_gene_count_tab`'s 2-column check |
| `NotImplemented` | `NotImplementedError` | `open_file` with a mode other than `r`/`w` on a `.gz` path |
| `Io` | `OSError` | any file open/read/write failure |
| `Format` | `ValueError` | FASTQ parse errors (Python calls `U.error()` → `sys.exit`; see deviation D6) |
| `UnknownClusterMethod` | `AttributeError` | `UmiClusterer::from_name("bogus")` (Python defers this to call time; deviation D2) |
| `UpstreamDefect` | `TypeError` | `CellClusterer::call` (deviation D3) |
| `Alignment` | `RuntimeError` | sort/write failures |

`panic!`/`unwrap`/`expect` appear only where a violated invariant is a
programmer error in *this crate* (e.g. an index derived from a length just
computed). No user input reaches a panic; every parse, every file open, every
option combination returns `Result`.

## API mapping table

Surface-type column: `fn`→function, `method`→method, `prop`→property-as-method,
`ctor`→constructor. Parity column: **A** = adaptation (shape changed, semantics
preserved), **D** = intentional deviation (behaviour changed, documented), **=**
= direct.

### `_dedup_umi.pyx`

| Python | Rust | Ownership | Errors | Surface | Parity |
|---|---|---|---|---|---|
| `edit_distance(a, b) -> int\|inf` | `edit_distance(a: &[u8], b: &[u8]) -> Option<u32>` | borrow | none | fn | **D1** |
| — | `edit_distance_f64(a, b) -> f64` | borrow | none | fn | A (helper for the averaging path) |
| — | `within_threshold(a, b, t: u32) -> bool` | borrow | none | fn | A (captures the `<= threshold` idiom once) |

### `network.py`

| Python | Rust | Ownership | Errors | Surface | Parity |
|---|---|---|---|---|---|
| `breadth_first_search(node, adj_list)` | `breadth_first_search(node: &[u8], adj: &AdjMap) -> BTreeSet<Vec<u8>>` | borrow → owned | none | fn | A (`set` → `BTreeSet`) |
| `recursive_search(node, adj_list)` | `recursive_search(node: &[u8], adj: &AdjMap) -> BTreeSet<Vec<u8>>` | borrow → owned | none | fn | A (function-attribute state → threaded accumulator) |
| `breadth_first_search_recursive(node, adj_list)` | same signature as above | borrow → owned | none | fn | A (`RecursionError` → explicit depth guard) |
| `remove_umis(adj_list, cluster, nodes)` | `remove_umis(adj, cluster: &[Vec<u8>], nodes: &[Vec<u8>]) -> BTreeSet<Vec<u8>>` | borrow | none | fn | = |
| `get_substr_slices(umi_length, idx_size)` | `get_substr_slices(umi_length: usize, idx_size: usize) -> Vec<(usize, usize)>` | owned | none | fn | = |
| `build_substr_idx(umis, umi_length, min_edit)` | `build_substr_idx(umis: &[Vec<u8>], umi_length, min_edit) -> SubstrIdx` | borrow → owned | none | fn | A (nested defaultdict → newtype) |
| `iter_nearest_neighbours(umis, substr_idx)` | `iter_nearest_neighbours(umis, idx) -> Vec<(Vec<u8>, Vec<u8>)>` | borrow → owned | none | fn | A (generator → eager `Vec`; see § Determinism) |
| `UMIClusterer(cluster_method="directional")` | `UmiClusterer::new(m: ClusterMethod)`, `UmiClusterer::from_name(&str) -> Result<Self>` | owned | `UnknownClusterMethod` | ctor | **D2** |
| `UMIClusterer.__call__(umis, threshold)` | `UmiClusterer::cluster(&mut self, counts: &CountMap, threshold: u32) -> Result<Vec<Vec<Vec<u8>>>>` | borrow, `&mut self` | `Assertion` | method | A (`&mut` because it updates position stats) |
| `.positions`, `.total_umis_per_position`, `.max_umis_per_position` | same-named `pub fn` accessors | borrow | none | prop | A |
| `ReadDeduplicator(options)` | `ReadDeduplicator::new(o: &DedupOptions) -> Result<Self>` | borrow → owned | `Io`, `Value` | ctor | A |
| `ReadDeduplicator.__call__(bundle, threshold)` | `ReadDeduplicator::dedup(&mut self, bundle: &Bundle, threshold) -> Result<DedupResult>` | borrow | `Assertion`, `Value` | method | A (3-tuple → named struct) |
| `CellClusterer(...)` | `CellClusterer::new(method, dir_threshold: f64, fuzzy_match: bool) -> Result<Self>` | owned | `Value` | ctor | = |
| `CellClusterer.__call__` | `CellClusterer::call(...)` → always `Err(UpstreamDefect)` | — | `UpstreamDefect` | method | **D3** |

### `umi_methods.py`

| Python | Rust | Ownership | Errors | Surface | Parity |
|---|---|---|---|---|---|
| `RANGES` | `pub const RANGES: [(&str, (u8, u8)); 3]` + `range_for(&str) -> Option<(u8,u8)>` | static | none | const | A |
| `Record(identifier, seq, quals, entry_format=None)` | `Record { identifier: String, seq: String, quals: String }` | owned | none | ctor | **D4** (`entry_format` is discarded by Python's own `__init__`) |
| `Record.guessFormat()` | `Record::guess_format(&self) -> Vec<&'static str>` | borrow | none | method | A (order fixed; see § Determinism) |
| `Record.__str__` | `impl Display for Record` | borrow | none | trait | A |
| `fastqIterate(infile, remove_suffix=False)` | `FastqIterator::new(r: impl BufRead, remove_suffix: bool)`, `Iterator<Item = Result<Record>>` | owned reader | `Format`, `Io` | fn→struct | **D6** |
| `joinedFastqIterate(it1, it2, strict=True)` | `JoinedFastqIterator::new(a, b, strict)`, `Iterator<Item = Result<(Record, Record)>>` | owned | `Value`, `Format` | fn→struct | A |
| `get_average_umi_distance(umis)` | `get_average_umi_distance(umis: &[Vec<u8>]) -> f64` | borrow | none | fn | = (returns `-1.0` for a single UMI, as Python does) |
| `random_read_generator(bamfile, chrom, barcode_getter)` | `RandomReadGenerator::new(path, chrom: Option<&str>, getter: &BarcodeGetter, seed) -> Result<Self>` | owned | `Io`, `Alignment` | ctor | A |
| `.getUmis(n)` | `RandomReadGenerator::get_umis(&mut self, n: usize) -> Vec<Vec<u8>>` | `&mut self` | none | method | = |

### `extract_methods.py`

| Python | Rust | Ownership | Errors | Surface | Parity |
|---|---|---|---|---|---|
| `addBarcodesToIdentifier(read, UMI, cell, umi_separator)` | `add_barcodes_to_identifier(read: &Record, umi: &str, cell: &str, sep: &str) -> String` | borrow | none | fn | = |
| `extractSeqAndQuals(seq, quals, umi_bases, cell_bases, discard_bases, retain_umi=False)` | `extract_seq_and_quals(...) -> SeqQuals` | borrow | none | fn | A (4-tuple → named struct) |
| `get_below_threshold(umi_quals, quality_encoding, quality_filter_threshold)` | `get_below_threshold(quals: &str, enc: QualityEncoding, thr: i32) -> Vec<bool>` | borrow | none | fn | A (`str` → enum) |
| `umi_below_threshold(...)` | `umi_below_threshold(...) -> bool` | borrow | none | fn | = |
| `mask_umi(umi, umi_quals, encoding, threshold)` | `mask_umi(...) -> String` | borrow | none | fn | = |
| `ExtractBarcodes(read, match, extract_umi, extract_cell, discard, retain_umi)` | `extract_barcodes(read: &Record, m: &RegexMatch, opts: ExtractBarcodeOpts) -> Result<Barcodes>` | borrow | `Value` | fn | A (4 bool params → options struct; 6-tuple → named struct) |
| `ExtractFilterAndUpdate(method, pattern, …, umi_separator="_")` (14 params) | `ExtractFilterAndUpdate::new(cfg: ExtractConfig) -> Result<Self>` with `ExtractConfig::default()` | owned | `Value` | ctor | A (14 positional params → config struct) |
| `.__call__(read1, read2=None)` | `process(&mut self, r1: &mut Record, r2: Option<&mut Record>) -> Result<bool>` | `&mut` | `Value` | method | A (mutates in place, returns "kept"; Python returns the mutated read(s) or `None`) |
| `.getBarcodes(read1, read2)` | `get_barcodes(&mut self, …) -> Result<Option<BarcodeValues>>` | `&mut` | `Value` | method | A (dispatch attribute → `match self.method`) |
| `.getCellBarcode(read1, read2)` | `get_cell_barcode(&mut self, …) -> Result<Option<String>>` | `&mut` | `Value` | method | A |
| `.filterQuality/.maskQuality/.filterCellBarcode/.filterUMIBarcode` | snake_case methods, same semantics | `&mut self` | none | method | = |
| `.getReadCounts()` | `read_counts(&self) -> &BTreeMap<String, u64>` | borrow | none | prop | A |

### `whitelist_methods.py`

| Python | Rust | Ownership | Errors | Surface | Parity |
|---|---|---|---|---|---|
| `getKneeEstimateDensity(counts, expect_cells=False, cell_number=False, plotfile_prefix=None)` | `knee_estimate_density(counts: &CellCounts, expect_cells: Option<u64>, cell_number: Option<usize>, plot_prefix: Option<&Path>) -> Result<KneeResult>` | borrow | `Io`, `Value` | fn | **D5** (`False`-as-absent → `Option`; plots → data) |
| `getKneeEstimateDistance(counts, cell_number=False, plotfile_prefix=None)` | `knee_estimate_distance(...) -> Result<KneeResult>` | borrow | `Value`, `Io` | fn | **D5** |
| `getErrorCorrectMapping(cell_barcodes, whitelist, threshold=1)` | `get_error_correct_mapping(barcodes: &[Vec<u8>], whitelist: &BTreeSet<Vec<u8>>, threshold: u32) -> BTreeMap<Vec<u8>, BTreeSet<Vec<u8>>>` | borrow | none | fn | **D7** (drops the `str(x)` coercion) |
| `getCellWhitelist(...)` | `get_cell_whitelist(...) -> Result<(Option<BTreeSet<Vec<u8>>>, Option<TrueToFalseMap>)>` | borrow | `Value`, `Io` | fn | A |
| `getUserDefinedBarcodes(tsv, tsv2=None, getErrorCorrection=False, deriveErrorCorrection=False, threshold=1)` | `get_user_defined_barcodes(tsv: &Path, tsv2: Option<&Path>, mode: ErrorCorrectionMode, threshold: u32) -> Result<(BTreeSet<Vec<u8>>, Option<FalseToTrueMap>)>` | borrow | `Io`, `Assertion`, `Format` | fn | A (two mutually-exclusive bools → enum) |
| `checkError(barcode, whitelist, errors=1)` | `check_error(barcode: &[u8], whitelist: &[Vec<u8>], errors: u32) -> Vec<Vec<u8>>` | borrow | none | fn | = (early return after 2 matches preserved) |
| `errorDetectAboveThreshold(counts, whitelist, t2f_map, errors=1, resolution_method="discard")` | `error_detect_above_threshold(…, method: ResolutionMethod) -> Result<(BTreeSet<Vec<u8>>, TrueToFalseMap)>` | borrow | `Assertion` | fn | = (including the never-incremented counter) |

`FalseToTrueMap` = `BTreeMap<Vec<u8>, Option<Vec<u8>>>` — the `None` value is
load-bearing: `getUserDefinedBarcodes` sets it for barcodes reachable from more
than one whitelist entry, and `filterUMIBarcode` branches on it.

### `sam_methods.py`

| Python | Rust | Ownership | Errors | Surface | Parity |
|---|---|---|---|---|---|
| `get_barcode_read_id(read, cell_barcode=False, sep="_")` | `get_barcode_read_id(name: &str, cell: bool, sep: &str) -> Result<Barcode>` | borrow | `Value` | fn | = |
| `get_barcode_tag(read, …8 params)` | `get_barcode_tag(read: &RecordBuf, cfg: &TagConfig) -> Result<Barcode>` | borrow | `Value` | fn | A |
| `get_umi_read_string(read_id, sep)` / `get_cell_umi_read_string(read_id, sep)` | same names, `(&str, &str) -> Result<Barcode>` | borrow | `Value` | fn | = |
| `get_barcode_umis(read, cell_barcode=False)` | `get_barcode_umis(name: &str, cell: bool) -> Result<Barcode>` | borrow | `Value` | fn | = |
| `get_gene_count_tab(infile, bc_getter)` | `GeneCountTab::new(r: impl BufRead, getter: StringBarcodeGetter)`, `Iterator<Item = Result<(String, CellUmiCounts)>>` | owned | `Assertion`, `Value` | fn→struct | A (callable param → enum, not `dyn Fn`) |
| `getMetaContig2contig(bamfile, map_path)` | `get_metacontig2contig(refs: &BTreeSet<String>, path: &Path) -> Result<MetaContigMap>` | borrow | `Io`, `Format` | fn | A (takes reference names, not the file handle) |
| `metafetcher(bamfile, m2c, metatag)` | `MetaFetcher` iterator | owned | `Alignment` | fn→struct | A |
| `find_splice(cigar)` | `find_splice(cigar: &[(CigarOp, usize)]) -> Result<Option<usize>>` | borrow | `Value` | fn | A (`int \| False` → `Option<usize>`) |
| `get_read_position(read: &RecordBuf, soft_clip_threshold)` | `get_read_position(...) -> Result<ReadPosition>` | borrow | `Value` | fn | A (3-tuple → named struct) |
| `get_bundles(options, …5 flags)` | `GetBundles::new(o: &BundleOptions, flags: BundleFlags) -> Result<Self>` | borrow → owned | `Value` | ctor | A |
| `.__call__(inreads)` | `GetBundles::bundles(&mut self, reads) -> impl Iterator<Item = Result<BundleItem>>` | `&mut self` | `Alignment`, `Value` | method | A (3-tuple yields → `enum BundleItem { Single(RecordBuf), Bundle{..} }`) |
| `.update_dicts(read, pos, key, umi)` / `.check_output()` | same names, snake_case | `&mut self` / `&self` | none | method | = |
| `TwoPassPairWriter(infile, outfile, tags=False)` | `TwoPassPairWriter::new(in_path: &Path, out: Writer) -> Result<Self>` | owned | `Io` | ctor | **D8** (`tags` param unused by Python's own `__init__`) |
| `.write(read, unique_id=None, umi=None, unmapped=False)` | `write(&mut self, read: &RecordBuf, unmapped: bool) -> Result<()>` | `&mut self` | `Io` | method | **D8** (`unique_id`/`umi` unused by Python's body) |
| `.write_mates()` / `.close()` | `write_mates(&mut self)` / `close(self)` | `&mut self` / consuming | `Io` | method | A |

### `Utilities.py` and the CLI

Row-by-row mapping in `07_translation_log.md` (mostly mechanical:
`log/info/warn/debug/error/critical` → same-named free functions over a
process-global logger; `openFile` → `open_file`/`create_file`, split by mode
because Rust's `File` is not mode-polymorphic; `determine_format`,
`output_names_and_formats`, `get_temp_filename` are direct). The `optparse`
classes map to `clap` plumbing in `crates/cli`.

`DefaultOptions` (a class used as a namespace of mutable globals) becomes
`GlobalOptions`, held in a `OnceLock<RwLock<GlobalOptions>>` in
`core::utilities` — the one piece of genuinely global mutable state in the
Python design, kept because `U.info()` is called from library code and must
respect `--log`/`--verbose` without threading a logger through 40 functions.

## Container decision checkpoint

1. **Does the package define a primary data container class?** Yes, two: the
   "bundle" (`dict[umi -> {"read", "count"}]`) and the FASTQ `Record`.
2. **Does the Rust ecosystem have a standard container?** For the alignment
   records themselves, **yes** — `noodles::sam::alignment::RecordBuf` is the
   domain-standard owned SAM/BAM record, and it is what `noodles`'s readers and
   writers both speak. For the bundle and the FASTQ record, **no** — these are
   package-specific aggregations.
3. **Consequence**: `RecordBuf` is adopted wholesale; no wrapper type, no
   re-derived field accessors. The `pysam` attribute → `noodles` accessor
   mapping is tabulated in `07_translation_log.md` (`read.pos` →
   `alignment_start()`, `read.aend` → `alignment_end()`, `read.tlen` →
   `template_length()`, `read.mapq` → `mapping_quality()`,
   `read.opt(tag)`/`get_tag` → `data().get()`, …).
4. For the bundle and `Record`, bespoke structs per `rust-porting-patterns.md`:

```rust
pub struct BundleEntry { pub reads: BundleReads, pub count: u64 }
pub enum BundleReads { One(Box<RecordBuf>), Many(Vec<RecordBuf>), CountOnly }
pub type Bundle = BTreeMap<Vec<u8>, BundleEntry>;
```

`BundleReads` is an enum rather than a `Vec` in every case because
`get_bundles` has three genuinely different retention modes (`all_reads`,
`only_count_reads`, single-best) that Python distinguishes only by which dict
keys it happens to populate — reading `bundle[umi]["read"]` is a `KeyError` in
`only_count_reads` mode. The enum makes that a compile-time distinction.

## Determinism decisions

Where Python's output order is an artefact of `set`/`dict` iteration, the port
picks a *defined* order and the differential test compares as sets.

| Site | Python order | Rust order | Test comparison |
|---|---|---|---|
| `breadth_first_search` return | `set` | `BTreeSet` | set equality |
| `_get_connected_components_adjacency` | `sorted(graph, key=counts, reverse=True)` then `component.sort()` | same | exact (Python already sorts) |
| `iter_nearest_neighbours` | `set` union, `difference_update` | `BTreeSet`, same difference | pair-set equality |
| `_group_adjacency` trailing members | `list(connected_nodes - observed)` — genuinely unordered | sorted | first element exact, remainder as sets |
| `Record.guessFormat` | `dict` iteration over `RANGES` | declaration order | set equality |
| `get_user_defined_barcodes` return | `set(whitelist)` | `BTreeSet` | set equality |

### Ordering that turned out to be observable, and is therefore reproduced

The Step-7 analysis treated Python's dictionary and set ordering as an artefact
to be replaced with a canonical order. Implementation showed that in four
places the order is **observable in the output the package's own test suite
compares**, so it is carried explicitly instead. Each was found by a
byte-for-byte comparison against a live run, not by re-reading the source.

| Site | What order Python uses | How the port carries it |
|---|---|---|
| A bundle's UMIs, feeding the clusterer | `dict` insertion order — first sighting at that position | [`bundles::BundleEntry::first_seen`] plus [`bundles::insertion_order`]; `UmiClusterer::cluster_ordered` takes the order explicitly |
| Component seeding in `_get_connected_components_adjacency` | `sorted(graph, key=counts.get, reverse=True)` — a *stable* sort over the adjacency dict's insertion order, which is the UMI order | the same UMI order is threaded into `connected_components_adjacency` |
| `Counter.most_common` for cell barcodes | stable sort over the counter's insertion order | [`whitelist::CellCounts`] is an ordered counter, not a `BTreeMap` |
| The whitelist handed to `errorDetectAboveThreshold` | a **list** in `most_common` order, not a set | [`whitelist::KneeResult::ordered_barcodes`], and the function takes a slice |

Where the order genuinely is an artefact — the trailing members of an
adjacency group, `breadth_first_search`'s return, `Record::guess_format` — the
port still picks a defined order and the differential tests compare as sets.

`sorted(x, key=counts.get, reverse=True)` in Python is a **stable** sort, so
equal-count UMIs keep insertion order. Rust's `sort_by` is also stable — the
port uses `sort_by` (never `sort_unstable_by`) at every count-keyed sort so tie
order matches given the same input order. This is load-bearing for
`_group_directional`, where the first element of each group is the reported
parent UMI.

## Intentional deviations

| # | Deviation | Python behaviour | Rust behaviour | Why |
|---|---|---|---|---|
| **D1** | `edit_distance` unequal lengths | `AttributeError` (NumPy ≥ 2 removed `np.Inf`, which the checked-in generated `.c` still references) | `None`, i.e. infinite distance, never `<= threshold` | Implements the `.pyx` gold standard's intent. The `AttributeError` is a build artefact of a stale generated file, not designed behaviour. |
| **D2** | Unknown cluster method | constructs fine; `AttributeError` on first call | `UmiClusterer::from_name` returns `Err(UnknownClusterMethod)` | Rust cannot leave methods unbound. Moving the error earlier is strictly more informative and cannot change a working program's behaviour. |
| **D3** | `CellClusterer.__call__` | always `TypeError` (three independent upstream bugs) | `Err(UpstreamDefect)` carrying the Python error text | The Python path is unreachable dead code. Guessing the intended algorithm would invent behaviour; erroring preserves "this does not work" without fabricating a result. The three working private methods *are* ported and public. |
| **D4** | `Record.__init__`'s `entry_format` | accepted, then discarded by a rebinding bug in the tuple assignment | not a constructor parameter | The parameter provably has no effect (`self.format` is never set; `guessFormat` recomputes). |
| **D5** | Whitelist plots | writes 3 PNGs via matplotlib | no PNGs; `_cell_thresholds.tsv` written identically, and `KneeResult` exposes every quantity the plots draw | Plot pixels are not library semantics and no Rust stack reproduces matplotlib's output bytes. The *data* is fully preserved and is what the test suite compares. |
| **D6** | FASTQ parse errors | `U.error()` → log CRITICAL + `sys.exit(1)` from inside a generator | `Err(PortError::Format)` | A library must not terminate the process. The CLI maps `Format` back to the same message and exit code, so CLI behaviour is unchanged. |
| **D7** | `getErrorCorrectMapping`'s `set([str(x) for x in whitelist])` | stringifies; given `bytes` produces `"b'ACGT'"` | takes `&BTreeSet<Vec<u8>>`, no coercion | The `str()` is defensive normalisation of an unknown input type, not intended semantics; reproducing it would mean reproducing a corrupted key. All in-package callers pass `str`, where the coercion is a no-op. |
| **D8** | `TwoPassPairWriter`'s `tags` ctor param and `write`'s `unique_id`/`umi` params | accepted and never read | not present | Provably dead parameters. |
| **D9** | Extension module name | n/a | `umi_tools_rs`, not `umi_tools._native` | The scaffold's mixed layout shadows the editable-installed oracle in the shared env (a hazard `differential-testing-and-validation.md` documents as confirmed). A distinct top-level name keeps both importable in one interpreter, which is what same-process differential testing actually requires. |
| **D10** | CLI help text layout | `optparse` + custom `BetterFormatter` | `clap` + a formatter reproducing the indent/wrap rules | Byte-identical `optparse` help is not achievable through `clap` and is not semantics. Option *names, arities, defaults and validation* are identical. |
| **D11** | *(withdrawn)* the `random.random()` tie-break was first recorded here as non-reproducible in Python | — | — | **This was wrong.** `Utilities.Start` calls `random.seed(options.random_seed)`, which the first reading of `Start` missed; the tie-break is fully reproducible given `--random-seed`, and the package's own test suite depends on it. The port therefore reproduces CPython's Mersenne Twister exactly (`core::mersenne`) rather than substituting a different generator. No deviation remains. Recorded rather than deleted because the mistaken analysis shaped the first implementation. |
| **D12** | Remote reference sequences | htslib fetches a `--reference-filename` given as an `http(s)`/`ftp` URL | [`core::alignment::is_remote_reference`] detects the URL and returns `PortError::Alignment` naming the workaround | Reading references over the network would mean an HTTP client and a TLS stack in a crate that otherwise touches no network. The gap is detected and reported with an actionable message rather than surfacing as "No such file or directory". Affects all three CRAM cases in the upstream suite, which pass a GitHub URL; with the same reference downloaded locally, `dedup` from CRAM is byte-identical to Python's output. |
| **D13** | CRAM **output** of records with no sequence | htslib writes them | `PortError::Alignment` naming the record | `noodles-cram`'s writer indexes past the end of an empty sequence while converting the CIGAR to CRAM features, and *panics*. A panic on ordinary user input is exactly what this port's error model forbids, so the case is detected before the writer sees it. CRAM *input* of the same records works, as does SAM/BAM output. Affects `dedup_to_cram` and `dedup_bam_to_cram`, whose fixture is entirely sequence-less records. |

## Staged-architecture module table

Assignments as of Step 7; revisited at the end of each slice (current values
tracked in `06_implementation_plan.md`).

| Python module | Rust module | Target stage | Rationale |
|---|---|---|---|
| `_dedup_umi.pyx` | `core::distance` | **C** | Trivially validated; a Python caller may want the fast path. |
| `network.py` | `core::network` | **C** | The algorithmic core and the most valuable thing to expose to Python users. |
| `umi_methods.py` | `core::umi_methods` | **C** | `Record`/FASTQ iteration is directly useful from Python. |
| `extract_methods.py` | `core::extract` | **C** | Same. |
| `whitelist_methods.py` | `core::whitelist` | **C** | Same. |
| `sam_methods.py` | `core::sam`, `core::bundles` | **C** | Bindings for the pure functions; the `RecordBuf`-typed iterators stay Rust-only, because marshalling `noodles` records into `pysam` objects would cost more than it buys. |
| `Utilities.py` (library half) | `core::utilities`, `core::alignment`, `core::options` | **D** | Logging/temp-file/format helpers have no value as a Python API — Python has better ones. |
| `Utilities.py` (optparse half) | `cli::parser`, `cli::help_fmt` | **D** | CLI-only. |
| `Documentation.py` | `cli::documentation` | **D** | CLI-only string constants. |
| the seven `main()`s + `umi_tools.py` | `cli::*` | **D** | Entry points; the deliverable is the `umi_tools` binary. |

No module is left at Stage A: every one has a real Rust implementation by the
end of the port. Stage B is the *transient* state each module passes through
during its slice — implemented in Rust, exercised against the live oracle
through `umi_tools_rs`, before its stage is advanced.

## Public API surface (`crates/core/src/lib.rs` re-exports)

```rust
pub mod bktree;  pub mod bundles; pub mod distance; pub mod error;
pub mod extract; pub mod fuzzy;   pub mod network;  pub mod numeric;
pub mod options; pub mod sam;     pub mod stats;    pub mod umi_methods;
pub mod utilities; pub mod whitelist;

pub use error::{PortError, Result};
pub use distance::{edit_distance, within_threshold};
pub use network::{CellClusterer, ClusterMethod, ReadDeduplicator, UmiClusterer};
pub use umi_methods::{FastqIterator, JoinedFastqIterator, Record};
pub use extract::{ExtractConfig, ExtractFilterAndUpdate};
pub use whitelist::{knee_estimate_density, knee_estimate_distance, KneeResult};
pub use bundles::{Bundle, BundleEntry, BundleItem, GetBundles};
pub const VERSION: &str = "1.1.6";
```

`pub(crate)`/private is the default for anything not listed. `numeric`, `bktree`
and `fuzzy` are `pub` because they are faithful reimplementations of
third-party behaviour a user may reasonably want to inspect or reuse, and
because their differential tests are part of the crate's evidence.
