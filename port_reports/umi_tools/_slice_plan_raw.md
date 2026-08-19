Nodes: 107
Edges: 66
Components: 47
SCCs: 107 (0 cyclic)
Max depth: 4
Density: 0.0058

# Implementation Plan (auto-generated)

Algorithm: **cohesion_aware**
Total slices: **5**

## Slice Definitions

### Slice 1: Classes + Functions

- **Exports** (30):
  - `AppendCommaOption` (CLASS)
  - `BetterFormatter` (CLASS)
  - `CellClusterer` (CLASS)
  - `DefaultOptions` (CLASS)
  - `ExtractFilterAndUpdate` (CLASS)
  - `MultiLineFormatter` (CLASS)
  - `OptionGroup` (CLASS)
  - `OptionParser` (CLASS)
  - `ReadDeduplicator` (CLASS)
  - `Record` (CLASS)
  - `TwoPassPairWriter` (CLASS)
  - `UMIClusterer` (CLASS)
  - `addBarcodesToIdentifier` (FUNCTION)
  - `aggregateStatsDF` (FUNCTION)
  - `breadth_first_search` (FUNCTION)
  - `callbackShortHelp` (FUNCTION)
  - `checkError` (FUNCTION)
  - `chunk_bam` (FUNCTION)
  - `copy_tags` (FUNCTION)
  - `critical` (FUNCTION)
  - `debug` (FUNCTION)
  - `detect_bam_features` (FUNCTION)
  - `determine_format` (FUNCTION)
  - `error` (FUNCTION)
  - `extractSeqAndQuals` (FUNCTION)
  - `find_splice` (FUNCTION)
  - `getErrorCorrectMapping` (FUNCTION)
  - `getFooter` (FUNCTION)
  - `getHeader` (FUNCTION)
  - `getKneeEstimateDensity` (FUNCTION)
- **Depends on**: (none)
- **Tests**: (to be filled)
- **Staged-architecture stage**: (to be filled — A/B/C/D, see staged-architecture.md)
- **Done when**: (to be filled)

### Slice 2: Classes + Functions

- **Exports** (30):
  - `getKneeEstimateDistance` (FUNCTION)
  - `getMetaContig2contig` (FUNCTION)
  - `getParams` (FUNCTION)
  - `getTempFile` (FUNCTION)
  - `getUserDefinedBarcodes` (FUNCTION)
  - `get_average_umi_distance` (FUNCTION)
  - `get_barcode_read_id` (FUNCTION)
  - `get_barcode_tag` (FUNCTION)
  - `get_barcode_umis` (FUNCTION)
  - `get_below_threshold` (FUNCTION)
  - `get_bundles` (CLASS)
  - `get_cell_umi_read_string` (FUNCTION)
  - `get_gene_count_tab` (FUNCTION)
  - `get_substr_slices` (FUNCTION)
  - `get_umi_read_string` (FUNCTION)
  - `info` (FUNCTION)
  - `iter_nearest_neighbours` (FUNCTION)
  - `joinedFastqIterate` (FUNCTION)
  - `log` (FUNCTION)
  - `main` (FUNCTION)
  - `metafetcher` (FUNCTION)
  - `openFile` (FUNCTION)
  - `open_input_alignments` (FUNCTION)
  - `pick_mate` (FUNCTION)
  - `random_read_generator` (CLASS)
  - `recursive_search` (FUNCTION)
  - `remove_umis` (FUNCTION)
  - `validateExtractOptions` (FUNCTION)
  - `warn` (FUNCTION)
  - `warning` (FUNCTION)
- **Depends on**: (none)
- **Tests**: (to be filled)
- **Staged-architecture stage**: (to be filled — A/B/C/D, see staged-architecture.md)
- **Done when**: (to be filled)

### Slice 3: Methods + Functions

- **Exports** (30):
  - `AppendCommaOption.convert_value` (METHOD)
  - `AppendCommaOption.take_action` (METHOD)
  - `BetterFormatter.__init__` (METHOD)
  - `BetterFormatter.format_description` (METHOD)
  - `BetterFormatter.format_epilog` (METHOD)
  - `BetterFormatter.format_option` (METHOD)
  - `BetterFormatter.format_usage` (METHOD)
  - `CellClusterer.__init__` (METHOD)
  - `ExtractBarcodes` (FUNCTION)
  - `ExtractFilterAndUpdate.__init__` (METHOD)
  - `ExtractFilterAndUpdate.filterCellBarcode` (METHOD)
  - `ExtractFilterAndUpdate.filterUMIBarcode` (METHOD)
  - `ExtractFilterAndUpdate.getReadCounts` (METHOD)
  - `MultiLineFormatter.format` (METHOD)
  - `OptionParser.__init__` (METHOD)
  - `ReadDeduplicator.__init__` (METHOD)
  - `Record.__init__` (METHOD)
  - `Record.guessFormat` (METHOD)
  - `Start` (FUNCTION)
  - `Stop` (FUNCTION)
  - `TwoPassPairWriter.__init__` (METHOD)
  - `TwoPassPairWriter.write_mates` (METHOD)
  - `UMIClusterer.__init__` (METHOD)
  - `breadth_first_search_recursive` (FUNCTION)
  - `build_substr_idx` (FUNCTION)
  - `errorDetectAboveThreshold` (FUNCTION)
  - `fastqIterate` (FUNCTION)
  - `getCellWhitelist` (FUNCTION)
  - `getTempFilename` (FUNCTION)
  - `get_bundles.__init__` (METHOD)
- **Depends on**: Slice 1: Classes + Functions, Slice 2: Classes + Functions
- **Tests**: (to be filled)
- **Staged-architecture stage**: (to be filled — A/B/C/D, see staged-architecture.md)
- **Done when**: (to be filled)

### Slice 4: Methods + Functions

- **Exports** (9):
  - `get_bundles.check_output` (METHOD)
  - `get_bundles.update_dicts` (METHOD)
  - `get_read_position` (FUNCTION)
  - `mask_umi` (FUNCTION)
  - `open_output_alignments` (FUNCTION)
  - `random_read_generator.refill_random` (METHOD)
  - `sort_output` (FUNCTION)
  - `umi_below_threshold` (FUNCTION)
  - `validateSamOptions` (FUNCTION)
- **Depends on**: Slice 1: Classes + Functions, Slice 2: Classes + Functions
- **Tests**: (to be filled)
- **Staged-architecture stage**: (to be filled — A/B/C/D, see staged-architecture.md)
- **Done when**: (to be filled)

### Slice 5: Methods + Functions

- **Exports** (8):
  - `ExtractFilterAndUpdate.filterQuality` (METHOD)
  - `ExtractFilterAndUpdate.maskQuality` (METHOD)
  - `TwoPassPairWriter.close` (METHOD)
  - `TwoPassPairWriter.write` (METHOD)
  - `output_names_and_formats` (FUNCTION)
  - `random_read_generator.__init__` (METHOD)
  - `random_read_generator.fill` (METHOD)
  - `random_read_generator.getUmis` (METHOD)
- **Depends on**: Slice 1: Classes + Functions, Slice 2: Classes + Functions, Slice 3: Methods + Functions, Slice 4: Methods + Functions
- **Tests**: (to be filled)
- **Staged-architecture stage**: (to be filled — A/B/C/D, see staged-architecture.md)
- **Done when**: (to be filled)

