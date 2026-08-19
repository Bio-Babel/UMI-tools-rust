# umi_tools-rs

A Rust port of [UMI-tools](https://github.com/CGATOxford/UMI-tools) 1.1.6 —
deduplicating, grouping and counting sequencing reads by Unique Molecular
Identifier, and extracting and whitelisting cell barcodes.

UMI-tools identifies PCR duplicates in sequencing data by grouping reads that
share a mapping position and a UMI, correcting for sequencing errors in the UMI
itself. This port reimplements it in Rust: a library crate with no Python
dependency, a drop-in `umi_tools` binary, and optional Python bindings.

## Status

The port is validated against the original by running the **upstream test
suite unchanged against the Rust binary**: 67 of its 70 end-to-end cases
produce byte-identical output. The three that do not are two documented
capability gaps, listed under [Known gaps](#known-gaps).

| Check | Result |
|---|---|
| `cargo test --workspace` | 396 passed |
| Upstream CLI acceptance suite | 67/70 byte-identical |
| Same-process differential harnesses | 40/40 comparisons pass |
| `cargo clippy --workspace --all-targets -- -D warnings` | clean |
| `cargo fmt --check` | clean |
| `cargo doc --no-deps --workspace` | clean |

## Install

The command-line tool:

```sh
cargo install umi_tools-cli
```

The library:

```sh
cargo add umi_tools-core
```

Or from a checkout:

```sh
cargo build --release          # target/release/umi_tools
cargo install --path crates/cli
```

## Quick start

Deduplicate an aligned BAM whose UMIs are in the read names — the pipeline's
usual final step:

```sh
umi_tools dedup --stdin=aligned.bam --stdout=deduplicated.bam --log=dedup.log
```

The seven subcommands are the original's:

```sh
umi_tools extract  --bc-pattern=NNNNNN --stdin=reads.fastq.gz --stdout=extracted.fastq.gz
umi_tools whitelist --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --stdin=R1.fastq.gz > whitelist.tsv
umi_tools group    --stdin=aligned.bam --group-out=groups.tsv --output-bam --stdout=tagged.bam
umi_tools dedup    --stdin=aligned.bam --stdout=deduplicated.bam
umi_tools count    --per-gene --gene-tag=XT --stdin=assigned.bam --stdout=counts.tsv
umi_tools count_tab --stdin=assigned.tsv --stdout=counts.tsv
umi_tools prepare-for-em --stdin=name_sorted.bam --stdout=for_rsem.bam
```

`umi_tools <tool> --help` prints the original's own help text.

### As a library

```rust
use std::collections::BTreeMap;
use umi_tools_core::{ClusterMethod, UmiClusterer};

// The UMIs seen at one alignment position, and how often each occurred.
let counts: BTreeMap<Vec<u8>, u64> = [
    (b"ACGT".to_vec(), 456),
    (b"ACGA".to_vec(), 90),   // one substitution from ACGT
    (b"TTTT".to_vec(), 200),
].into_iter().collect();

let mut clusterer = UmiClusterer::new(ClusterMethod::Directional);
let groups = clusterer.cluster(&counts, 1)?;

// ACGA is absorbed into ACGT as a likely sequencing error; TTTT stands alone.
assert_eq!(groups.len(), 2);
assert_eq!(groups[0], vec![b"ACGT".to_vec(), b"ACGA".to_vec()]);
# Ok::<(), umi_tools_core::PortError>(())
```

## Runnable examples

Each reproduces a step of the original's documented tutorials against the
checked-in fixtures:

| Example | What it shows |
|---|---|
| [`cluster_umis.rs`](crates/core/examples/cluster_umis.rs) | all five clustering methods on one position, and why they differ |
| [`extract_barcodes.rs`](crates/core/examples/extract_barcodes.rs) | positional and regex barcode extraction, plus quality masking |
| [`deduplicate_bam.rs`](crates/core/examples/deduplicate_bam.rs) | the full dedup pipeline over a real BAM |
| [`cell_whitelist.rs`](crates/core/examples/cell_whitelist.rs) | both knee estimators and error correction |

```sh
cargo run --example deduplicate_bam
cargo run --example deduplicate_bam -- path/to/your.bam
```

## Python bindings

The clustering, extraction and whitelisting APIs are also callable from Python:

```sh
cd crates/pyo3-bindings && maturin develop --release
```

```python
import umi_tools_rs

clusterer = umi_tools_rs.UMIClusterer("directional")
clusterer.cluster({b"ACGT": 10, b"ACGA": 2}, 1)
# [[b'ACGT', b'ACGA']]
```

The extension installs as a **separate top-level package**, `umi_tools_rs`,
rather than as `umi_tools._native`. That is deliberate: installing it as a
submodule of a package that is separately `pip install -e .`-ed shadows the
original, which is what makes the two importable side by side here — and what
lets `validation/` compare them in one interpreter with no serialisation
round-trip.

## Version mapping

The crate version tracks the Python package version it was ported from.

| Python `umi_tools` | Rust crates |
|---|---|
| 1.1.6 | `umi_tools-core` 1.1.6, `umi_tools-cli` 1.1.6, `umi_tools_rs` 1.1.6 |

A future port of a newer upstream release takes the matching version number.
Patch releases of the crates that do not correspond to an upstream release use
a fourth component or a pre-release suffix, so the correspondence stays
readable.

## Known gaps

Two capability gaps, both detected and reported with an actionable message
rather than failing obscurely. They account for the 3 upstream cases that do
not pass.

- **Remote reference sequences.** htslib fetches a `--reference-filename` given
  as an `http(s)`/`ftp` URL; this port reads local files only. Download the
  reference and its `.fai` index, and pass the local path. With the reference
  local, `dedup` from CRAM produces byte-identical output.
- **CRAM output of records with no sequence.** A SAM record may carry `*` for
  its sequence, and `noodles-cram`'s writer cannot encode that. CRAM *input* of
  such records works, as does SAM and BAM output.

Every intentional behaviour difference is enumerated in
[`05_design.md`](../port_reports/umi_tools/05_design.md) § Intentional
deviations, with the reasoning for each.

## Layout

```
crates/core/            the library: algorithms, no CLI, no Python
crates/cli/             the `umi_tools` binary
crates/pyo3-bindings/   the `umi_tools_rs` Python extension
validation/             differential harnesses comparing against live Python
```

## Licence

MIT, as the original.
