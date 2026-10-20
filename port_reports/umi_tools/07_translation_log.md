# 07 — Translation Log

A running record of decisions taken while implementing, and of what running the
Python oracle revealed that reading it did not. The API mapping tables live in
`05_design.md`; this file is the narrative and the surprises.

## Session — 2026-08-18

- Resume point: fresh port, Step 0.
- Build status at close: `cargo build --workspace` clean, `cargo clippy
  --workspace --all-targets -- -D warnings` clean, `cargo fmt --check` clean,
  `cargo test --workspace` 395 passed / 0 failed.
- Upstream CLI acceptance suite, run unchanged against the Rust binary:
  **67 of 70 passing**.

---

## Defects found by *running* the Python oracle, not by reading it

Each of these changed the implementation. None was visible from the source
alone; every one surfaced from a differential comparison against a live run.

### 1. `_dedup_umi.edit_distance` raises under NumPy 2

`_dedup_umi.pyx` returns `np.inf` for unequal-length inputs, but the
**checked-in generated `_dedup_umi.c`** — which is what `pyproject.toml`
compiles — emits `np.Inf`, removed in NumPy 2.0. So the shipped extension
raises `AttributeError` rather than returning infinity. The port implements the
`.pyx` intent (deviation D1). Found by calling `edit_distance(b"ACGT", b"ACG")`
in the oracle.

### 2. `CellClusterer.__call__` cannot run at all

Three independent defects: it passes four arguments to a three-parameter
method, then calls `self.get_groups`, which `__init__` never assigns, and
`_get_connected_components_adjacency` calls `.sort()` on the `set` that
`breadth_first_search` returns. Confirmed by calling it:
`TypeError: CellClusterer._get_connected_components_adjacency() takes 3
positional arguments but 4 were given`. Ported as
`Err(PortError::UpstreamDefect)` (deviation D3); the three working private
methods are ported and public.

### 3. `noodles` reports MAPQ 255 as `None`; `pysam` returns 255

Every read in the package's own `chr19.bam` carries MAPQ 255, the SAM
specification's "not available" sentinel. Treating `noodles`' `None` as 0 made
`--mapping-quality=20` discard the **entire file**. Fixed by
[`sam::raw_mapping_quality`], which maps `None` back to 255. Found by the
`bundles_mapq20` scenario in `crates/core/tests/sam_parity.rs` producing zero
bundles where Python produced 41.

### 4. The log-level comparison was inverted

`LogLevel` is ordered by verbosity (`Error < Info < Debug`), and `emit` tested
`level < loglevel` — so at the default verbosity **every error and warning was
silently dropped**. The symptom was that every CLI failure printed only "check
the log file" and the log contained nothing. Found by running the upstream
suite against the Rust binary and having no diagnosis available.

### 5. `Utilities.Start` seeds the stdlib RNG — the first reading missed it

`Start` calls `random.seed(options.random_seed)` at line 1110. The Step-5
analysis concluded the opposite (that only `np.random.seed` was called, which
does not affect the stdlib `random` module) and recorded it as deviation D11 —
"Python is non-reproducible here". It is not: three runs of `umi_tools dedup
--random-seed=123456789` produce byte-identical output, and the checked-in
reference matches. D11 is withdrawn in `05_design.md`, and the port reproduces
CPython's Mersenne Twister exactly (`core::mersenne`).

### 6. Python sorts its output even when writing to a pipe

`output_names_and_formats` routes records through a temporary file whenever
`--no-sort-output` is absent, *including* when the eventual output is standard
output — then sorts that file to the pipe. The port initially skipped the sort
for stdout, producing the right 728 reads in the wrong order. Found by
`diff <(sort rust) <(sort python)` being empty while the unsorted diff was not.

### 7. Four ordering sites are observable, not artefacts

See `05_design.md` § Ordering that turned out to be observable. Each was found
the same way: identical *content*, different *order*, in a suite that compares
by md5.

---

## Dependency-mapping decisions taken during implementation

### `core::pattern` — a new module, not in the Step-7 design

`--extract-method=regex` needs two things the Rust `regex` crate does not
provide: **per-group spans** (which bases each named group covered, so the read
can be sliced) and **approximate matching** (`{s<=n}` and friends). No
crates.io crate offers the second. A focused backtracking matcher over the
pattern subset the tool uses was the answer.

The critical semantic finding, again from running rather than reading: in

```text
(?P<discard_2>GAGTGATTGCTTGTGACGCCTT{s<=2})
```

the `{s<=2}` binds to the **final `T` alone**, not to the 22-base linker.
Substituting the linker's last base still matches; substituting its first does
not; substituting the last *two* does not either. Reading it as "the linker
within 2 substitutions" — the obvious reading — would have made the inDrop
patterns match far more reads than they should.

Validated by 1715 pattern/subject comparisons against the live `regex` module,
covering every `--bc-pattern` in the test suite and documentation over 245 real
reads. All match, spans included.

### `core::mersenne` — a new module, for the same reason

Reproducing the *distribution* of `random.random()` is not enough when the
output is compared byte for byte; the *stream* has to match. CPython's
`random.seed(n)` runs `init_by_array` over the seed's 32-bit words, while
NumPy's `RandomState.seed(n)` runs `init_genrand` on the low 32 bits — the same
algorithm, seeded two different ways, and `dedup --output-stats` uses both.
Both are implemented and verified against the live interpreters.

### `pysam` → `noodles`

`dedup`, `group` and `count` hold reads as `record::BamRecord` — the record's
BAM encoding, which is also what a `pysam.AlignedSegment` wraps — and read
fields through `noodles_bam::RecordRef`. Other commands use `noodles`'
decoded `RecordBuf`.

| `pysam` | Rust | Note |
|---|---|---|
| `read.qname` | `record.name()` | `None` for `*` |
| `read.pos` | `record.position()` | 0-based, as pysam |
| `read.aend` | `pos` + the CIGAR's reference span | computed in `sam::get_read_position` |
| `read.mapq` | `record.mapping_quality()` | the raw byte, 255 included; see defect 3 above |
| `read.cigar` | `record.view().cigar()` | `(Kind, usize)` pairs rather than `(int, int)` |
| `read.tlen` | `record.template_length()` | |
| `read.is_reverse` | `record.flags().is_reverse_complemented()` | |
| `read.is_read2` | `flags.is_segmented() && flags.is_last_segment()` | pysam's `is_read2` implies paired |
| `read.get_tag(t)` / `read.opt(t)` | `sam::tag_text(record, t)` | coerces any scalar type to text, as pysam does |
| `read.set_tag(t, v)` | `record.set_tag(t, v)` | moves an existing tag to the end and picks the smallest integer type, as pysam does |
| `AlignmentFile.references` | `header.reference_sequences().keys()` | |
| `AlignmentFile.fetch()` | `OpenedAlignment::records(until_eof)` | streamed; `fetch()` skips the unplaced block, `fetch(until_eof=True)` does not |
| `AlignmentFile.fetch(reference=c)` | `OpenedAlignment::fetch_reference(id)` | seeks with a `.bai`/`.csi` index; without one reads forward, where pysam refuses |
| `AlignmentFile(..., template=infile)` | `open_output_alignments(.., header_text, ..)` | SAM and BAM headers are the input's text, as htslib carries it over |
| `pysam.sort(...)` | `alignment::CoordinateSorter` | samtools' order and memory budget, spilling sorted runs to temporary files |

### `optparse` → `clap`

Behaviours that had to be reproduced explicitly, each found by a failing case:

- **Unambiguous long-option prefixes.** `optparse` accepts `--reference-file`
  for `--reference-filename`, and the test suite relies on it.
  `Command::infer_long_args(true)`.
- **Repeated flags.** `store_true` accepts repeats silently; `clap` errors. One
  `group` case passes `--out-sam` twice. Every option gets
  `overrides_with(itself)`, making repeats last-wins.
- **`-v/--verbose` takes an integer**, not a bare flag.
- **Two spellings, one destination.** `--buffer-whole-contig` and
  `--whole-contig` both set `whole_contig`; the second is an alias.
- **`--help-extended`** is declared with `action='help'` — a second spelling of
  `--help`, not a longer output.
- **Real defaults**: `chimeric_pairs="use"` and `unpaired_reads="use"`, not
  `"discard"`. Read off a live run's own parameter log, which prints every
  resolved option — a far more reliable source than the `add_option` calls
  scattered across four blocks.

---

## API shape decisions

- **`UmiClusterer::cluster_ordered`** was added when ordering proved
  observable. `cluster` remains, delegating to it with the count map's own
  order, so a library user gets a deterministic default and a caller with the
  original observation order can supply it.
- **`ExtractFilterAndUpdate::process` rewrites in place and returns `bool`**
  rather than returning the mutated read or `None`. Python returns the read;
  Rust would have to clone it. The discard case becomes a `bool` rather than an
  `Option` the caller must destructure.
- **`TwoPassPairWriter::close` consumes `self`.** The Python docstring warns
  that "if close() is not called, at least a contig's worth of mates will be
  missing"; making it consuming turns that warning into a compile error.
- **`BundleReads` is an enum**, not a `Vec`. Python's three retention modes are
  distinguished only by which dictionary keys it populates, so reading
  `bundle[umi]["read"]` in `only_count_reads` mode is a `KeyError`; the enum
  makes that a compile-time distinction.
- **`CellUmiCounts` and `CellCounts` are ordered containers**, not maps. Both
  became so after ordering proved observable — the third and fourth time, which
  is why the order lives in the type rather than at each call site.

## Deliberate non-goals

- **Matplotlib plots** are not reproduced (deviation D5); the data behind them
  is, including the `_cell_thresholds.tsv` sidecar the test suite compares.
- **`crates/cli` is not a library.** The seven `main()` functions and the
  option parser live there precisely so `crates/core` has no notion of a
  command line — a split the Python package cannot make, since every module
  imports `Utilities` and its process globals.
