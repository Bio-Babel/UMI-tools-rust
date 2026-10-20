# 02 — Feature Checklist (Step 3)

Mode: **full-package** (default). Every public export is in scope. No items are
deferred; no user agreement to reduce scope was sought or given.

Derived from the Step 2 audit (`01_audit.md`, `_audit_raw.md`), not from a
template. Names are qualified by module because five modules each define a
`main`.

Status vocabulary: `not_started` → `in_progress` → `implemented` (compiles) →
`validated` (differential parity demonstrated + deviations documented).

**Final status: every one of the 133 items is `validated`.** Evidence per item
is in `08_validation.md`; the three kinds are (a) a same-process differential
comparison in `validation/results_*.csv`, (b) a fixture parity test in
`crates/core/tests/*.rs` against values dumped from the live Python package, or
(c) the upstream CLI acceptance suite run unchanged against the Rust binary,
which exercises every entry point end to end and compares its output byte for
byte. Nothing is `deferred`.

Three items are `validated` in the sense of "the port's behaviour is the
documented and verified one", where the Python behaviour is a defect rather
than something to reproduce: `edit_distance` on unequal lengths (D1),
`CellClusterer.__call__` (D3), and `errorDetectAboveThreshold`'s
never-incremented counter (mirrored exactly, defect and all).

Risk markers: `dyn` = dynamic-dispatch/duck-typed signature needing explicit
Step-7 design; `rng` = non-deterministic; `float` = numerically sensitive;
`io` = alignment/FASTQ I/O; `defect` = upstream bug mirrored or corrected.

---

## Native core (`_dedup_umi.pyx`)

| Python symbol | Module | Rust target | Risk | Status |
|---|---|---|---|---|
| `edit_distance` | `_dedup_umi.pyx` | `core::distance::edit_distance` | defect | validated |

## `network.py` — clustering core

| Python symbol | Module | Rust target | Risk | Status |
|---|---|---|---|---|
| `breadth_first_search` | `network.py` | `core::network::breadth_first_search` | | validated |
| `recursive_search` | `network.py` | `core::network::recursive_search` | | validated |
| `breadth_first_search_recursive` | `network.py` | `core::network::breadth_first_search_recursive` | | validated |
| `remove_umis` | `network.py` | `core::network::remove_umis` | | validated |
| `get_substr_slices` | `network.py` | `core::network::get_substr_slices` | | validated |
| `build_substr_idx` | `network.py` | `core::network::build_substr_idx` | | validated |
| `iter_nearest_neighbours` | `network.py` | `core::network::iter_nearest_neighbours` | | validated |
| `UMIClusterer` | `network.py` | `core::network::UmiClusterer` | dyn | validated |
| `UMIClusterer.__init__` | `network.py` | `UmiClusterer::new` | dyn | validated |
| `UMIClusterer.__call__` | `network.py` | `UmiClusterer::cluster` | | validated |
| `UMIClusterer._get_best_min_account` | `network.py` | `UmiClusterer::get_best_min_account` | | validated |
| `UMIClusterer._get_best_percentile` | `network.py` | `UmiClusterer::get_best_percentile` | float | validated |
| `UMIClusterer._get_adj_list_adjacency` | `network.py` | `UmiClusterer::adj_list_adjacency` | | validated |
| `UMIClusterer._get_adj_list_directional` | `network.py` | `UmiClusterer::adj_list_directional` | | validated |
| `UMIClusterer._get_adj_list_null` | `network.py` | `AdjList::None` variant | | validated |
| `UMIClusterer._get_connected_components_adjacency` | `network.py` | `UmiClusterer::connected_components_adjacency` | | validated |
| `UMIClusterer._get_connected_components_null` | `network.py` | `Clusters::Flat` variant | | validated |
| `UMIClusterer._group_unique` | `network.py` | `UmiClusterer::group_unique` | | validated |
| `UMIClusterer._group_directional` | `network.py` | `UmiClusterer::group_directional` | | validated |
| `UMIClusterer._group_adjacency` | `network.py` | `UmiClusterer::group_adjacency` | | validated |
| `UMIClusterer._group_cluster` | `network.py` | `UmiClusterer::group_cluster` | | validated |
| `UMIClusterer._group_percentile` | `network.py` | `UmiClusterer::group_percentile` | float | validated |
| `ReadDeduplicator` | `network.py` | `core::network::ReadDeduplicator` | dyn | validated |
| `ReadDeduplicator.__init__` | `network.py` | `ReadDeduplicator::new` | | validated |
| `ReadDeduplicator.__call__` | `network.py` | `ReadDeduplicator::dedup` | | validated |
| `CellClusterer` | `network.py` | `core::network::CellClusterer` | defect | validated |
| `CellClusterer.__init__` | `network.py` | `CellClusterer::new` | | validated |
| `CellClusterer.__call__` | `network.py` | `CellClusterer::call` (returns `Err(UpstreamDefect)`) | defect | validated |
| `CellClusterer._get_adj_list_directional` | `network.py` | `CellClusterer::adj_list_directional` | | validated |
| `CellClusterer._get_connected_components_adjacency` | `network.py` | `CellClusterer::connected_components_adjacency` | defect | validated |
| `CellClusterer._get_best_min_account` | `network.py` | `CellClusterer::get_best_min_account` | | validated |

## `umi_methods.py`

| Python symbol | Module | Rust target | Risk | Status |
|---|---|---|---|---|
| `RANGES` | `umi_methods.py` | `core::umi_methods::RANGES` | | validated |
| `Record` | `umi_methods.py` | `core::umi_methods::Record` | | validated |
| `Record.guessFormat` | `umi_methods.py` | `Record::guess_format` | | validated |
| `Record.__str__` | `umi_methods.py` | `impl Display for Record` | | validated |
| `fastqIterate` | `umi_methods.py` | `core::umi_methods::FastqIterator` | io | validated |
| `joinedFastqIterate` | `umi_methods.py` | `core::umi_methods::JoinedFastqIterator` | io | validated |
| `get_average_umi_distance` | `umi_methods.py` | `core::umi_methods::get_average_umi_distance` | float | validated |
| `random_read_generator` | `umi_methods.py` | `core::umi_methods::RandomReadGenerator` | rng, io | validated |
| `random_read_generator.fill` | `umi_methods.py` | `cli::dedup::umi_frequencies` | io | validated |
| `random_read_generator.refill_random` | `umi_methods.py` | `RandomReadGenerator::refill_random` | rng | validated |
| `random_read_generator.getUmis` | `umi_methods.py` | `RandomReadGenerator::get_umis` | rng | validated |

## `extract_methods.py`

| Python symbol | Module | Rust target | Risk | Status |
|---|---|---|---|---|
| `addBarcodesToIdentifier` | `extract_methods.py` | `core::extract::add_barcodes_to_identifier` | | validated |
| `extractSeqAndQuals` | `extract_methods.py` | `core::extract::extract_seq_and_quals` | | validated |
| `get_below_threshold` | `extract_methods.py` | `core::extract::get_below_threshold` | | validated |
| `umi_below_threshold` | `extract_methods.py` | `core::extract::umi_below_threshold` | | validated |
| `mask_umi` | `extract_methods.py` | `core::extract::mask_umi` | | validated |
| `ExtractBarcodes` | `extract_methods.py` | `core::extract::extract_barcodes` | | validated |
| `ExtractFilterAndUpdate` | `extract_methods.py` | `core::extract::ExtractFilterAndUpdate` | dyn | validated |
| `ExtractFilterAndUpdate.__init__` | `extract_methods.py` | `ExtractFilterAndUpdate::new` | dyn | validated |
| `ExtractFilterAndUpdate.__call__` | `extract_methods.py` | `ExtractFilterAndUpdate::process` | | validated |
| `ExtractFilterAndUpdate.getBarcodes` | `extract_methods.py` | `ExtractFilterAndUpdate::get_barcodes` | dyn | validated |
| `ExtractFilterAndUpdate.getCellBarcode` | `extract_methods.py` | `ExtractFilterAndUpdate::get_cell_barcode` | dyn | validated |
| `ExtractFilterAndUpdate.filterQuality` | `extract_methods.py` | `ExtractFilterAndUpdate::filter_quality` | | validated |
| `ExtractFilterAndUpdate.maskQuality` | `extract_methods.py` | `ExtractFilterAndUpdate::mask_quality` | | validated |
| `ExtractFilterAndUpdate.filterCellBarcode` | `extract_methods.py` | `ExtractFilterAndUpdate::filter_cell_barcode` | | validated |
| `ExtractFilterAndUpdate.filterUMIBarcode` | `extract_methods.py` | `ExtractFilterAndUpdate::filter_umi_barcode` | | validated |
| `ExtractFilterAndUpdate.getReadCounts` | `extract_methods.py` | `ExtractFilterAndUpdate::read_counts` | | validated |

## `whitelist_methods.py`

| Python symbol | Module | Rust target | Risk | Status |
|---|---|---|---|---|
| `getKneeEstimateDensity` | `whitelist_methods.py` | `core::whitelist::knee_estimate_density` | float | validated |
| `getKneeEstimateDistance` | `whitelist_methods.py` | `core::whitelist::knee_estimate_distance` | float | validated |
| `getErrorCorrectMapping` | `whitelist_methods.py` | `core::whitelist::get_error_correct_mapping` | | validated |
| `getCellWhitelist` | `whitelist_methods.py` | `core::whitelist::get_cell_whitelist` | float | validated |
| `getUserDefinedBarcodes` | `whitelist_methods.py` | `core::whitelist::get_user_defined_barcodes` | io | validated |
| `checkError` | `whitelist_methods.py` | `core::whitelist::check_error` | | validated |
| `errorDetectAboveThreshold` | `whitelist_methods.py` | `core::whitelist::error_detect_above_threshold` | defect | validated |

## `sam_methods.py`

| Python symbol | Module | Rust target | Risk | Status |
|---|---|---|---|---|
| `get_barcode_read_id` | `sam_methods.py` | `core::sam::get_barcode_read_id` | | validated |
| `get_barcode_tag` | `sam_methods.py` | `core::sam::get_barcode_tag` | | validated |
| `get_umi_read_string` | `sam_methods.py` | `core::sam::get_umi_read_string` | | validated |
| `get_cell_umi_read_string` | `sam_methods.py` | `core::sam::get_cell_umi_read_string` | | validated |
| `get_barcode_umis` | `sam_methods.py` | `core::sam::get_barcode_umis` | | validated |
| `get_gene_count_tab` | `sam_methods.py` | `core::sam::get_gene_count_tab` | io | validated |
| `getMetaContig2contig` | `sam_methods.py` | `core::sam::get_metacontig2contig` | io | validated |
| `metafetcher` | `sam_methods.py` | `core::sam::metafetch` (`MetaFetcher`) | io | validated |
| `find_splice` | `sam_methods.py` | `core::sam::find_splice` | | validated |
| `get_read_position` | `sam_methods.py` | `core::sam::get_read_position` | | validated |
| `get_bundles` | `sam_methods.py` | `core::sam::GetBundles` | dyn, rng, io | validated |
| `get_bundles.__init__` | `sam_methods.py` | `GetBundles::new` | dyn | validated |
| `get_bundles.update_dicts` | `sam_methods.py` | `GetBundles::update_dicts` | rng | validated |
| `get_bundles.check_output` | `sam_methods.py` | `GetBundles::check_output` | | validated |
| `get_bundles.__call__` | `sam_methods.py` | `GetBundles::bundles` | io | validated |
| `TwoPassPairWriter` | `sam_methods.py` | `core::sam::TwoPassPairWriter` | io | validated |
| `TwoPassPairWriter.write` | `sam_methods.py` | `TwoPassPairWriter::write` | io | validated |
| `TwoPassPairWriter.write_mates` | `sam_methods.py` | `TwoPassPairWriter::write_mates` | io | validated |
| `TwoPassPairWriter.close` | `sam_methods.py` | `TwoPassPairWriter::close` | io | validated |

## `Utilities.py`

| Python symbol | Module | Rust target | Risk | Status |
|---|---|---|---|---|
| `DefaultOptions` | `Utilities.py` | `core::utilities::GlobalOptions` | | validated |
| `BetterFormatter` | `Utilities.py` | `cli::help_fmt::BetterFormatter` | | validated |
| `AppendCommaOption` | `Utilities.py` | `cli::options::split_append_value` | | validated |
| `OptionParser` | `Utilities.py` | `cli::parser::build_parser` | dyn | validated |
| `OptionGroup` | `Utilities.py` | `clap::Command` group (alias) | | validated |
| `MultiLineFormatter` | `Utilities.py` | `core::utilities::format_multiline` | | validated |
| `callbackShortHelp` | `Utilities.py` | `cli::parser::short_help` | | validated |
| `openFile` | `Utilities.py` | `core::utilities::open_file` / `create_file` | io | validated |
| `getHeader` | `Utilities.py` | `core::utilities::get_header` | | validated |
| `getParams` | `Utilities.py` | `core::utilities::get_params` | | validated |
| `getFooter` | `Utilities.py` | `core::utilities::get_footer` | | validated |
| `Start` | `Utilities.py` | `cli::start::start` | dyn | validated |
| `Stop` | `Utilities.py` | `cli::start::stop` | | validated |
| `validateExtractOptions` | `Utilities.py` | `core::options::validate_extract_options` | | validated |
| `validateSamOptions` | `Utilities.py` | `core::options::validate_sam_options` | | validated |
| `log` | `Utilities.py` | `core::utilities::log` | | validated |
| `info` | `Utilities.py` | `core::utilities::info` | | validated |
| `warning` | `Utilities.py` | `core::utilities::warning` | | validated |
| `warn` | `Utilities.py` | `core::utilities::warn` | | validated |
| `debug` | `Utilities.py` | `core::utilities::debug` | | validated |
| `error` | `Utilities.py` | `core::utilities::error` | | validated |
| `critical` | `Utilities.py` | `core::utilities::critical` | | validated |
| `getTempFile` | `Utilities.py` | `core::utilities::get_temp_file` | io | validated |
| `getTempFilename` | `Utilities.py` | `core::utilities::get_temp_filename` | io | validated |
| `determine_format` | `Utilities.py` | `core::utilities::determine_format` | | validated |
| `output_names_and_formats` | `Utilities.py` | `core::utilities::output_names_and_formats` | | validated |
| `open_input_alignments` | `Utilities.py` | `core::alignment::open_input_alignments` | io | validated |
| `open_output_alignments` | `Utilities.py` | `core::alignment::open_output_alignments` | io | validated |
| `sort_output` | `Utilities.py` | `core::alignment::sort_output` | io | validated |

## Subcommand entry points

| Python symbol | Module | Rust target | Risk | Status |
|---|---|---|---|---|
| `main` | `umi_tools.py` | `cli::dispatch::main` | | validated |
| `main` | `extract.py` | `cli::extract::main` | io | validated |
| `main` | `whitelist.py` | `cli::whitelist::main` | float, io | validated |
| `main` | `group.py` | `cli::group::main` | io | validated |
| `main` | `dedup.py` | `cli::dedup::main` | rng, io | validated |
| `main` | `count.py` | `cli::count::main` | io | validated |
| `main` | `count_tab.py` | `cli::count_tab::main` | io | validated |
| `main` | `prepare_for_em.py` | `cli::prepare_for_em::main` | io | validated |
| `detect_bam_features` | `dedup.py` | `cli::dedup::detect_bam_features` | io | validated |
| `aggregateStatsDF` | `dedup.py` | `core::stats::aggregate_stats` | float | validated |
| `chunk_bam` | `prepare_for_em.py` | `cli::prepare_for_em::chunk_bam` | io | validated |
| `copy_tags` | `prepare_for_em.py` | `cli::prepare_for_em::copy_tags` | | validated |
| `pick_mate` | `prepare_for_em.py` | `cli::prepare_for_em::pick_mate` | | validated |

## `Documentation.py`

| Python symbol | Module | Rust target | Risk | Status |
|---|---|---|---|---|
| `GENERIC_DOCSTRING_GDC` | `Documentation.py` | `cli::documentation::GENERIC_DOCSTRING_GDC` | | validated |
| `GROUP_DEDUP_GENERIC_OPTIONS` | `Documentation.py` | `cli::documentation::GROUP_DEDUP_GENERIC_OPTIONS` | | validated |
| `GENERIC_DOCSTRING_WE` | `Documentation.py` | `cli::documentation::GENERIC_DOCSTRING_WE` | | validated |
| `GENERIC_DOCSTRING_SBCRAM_INPUT` | `Documentation.py` | `cli::documentation::GENERIC_DOCSTRING_SBCRAM_INPUT` | | validated |
| `GENERIC_DOCSTRING_SBCRAM_OUTPUT` | `Documentation.py` | `cli::documentation::GENERIC_DOCSTRING_SBCRAM_OUTPUT` | | validated |

## `version.py`

| Python symbol | Module | Rust target | Risk | Status |
|---|---|---|---|---|
| `__version__` | `version.py` | `core::VERSION` | | validated |
