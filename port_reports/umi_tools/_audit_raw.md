# Audit: umi_tools

## Package

- Path: `$UMI_TOOLS_SRC/umi_tools`
- Import name: `umi_tools`
- Version: `1.1.6`
- Description: None
- `__all__` declared: no (public surface inferred from leading underscore)

## Dependencies

- Runtime: (none)

## Public API Summary

- Public functions: 69
- Public classes: 14

### Public Functions

- `callbackShortHelp(option, opt, value, parser)` (`Utilities.py`) _untyped_
- `openFile(filename, mode, create_dir)` (`Utilities.py`) _untyped_
- `getHeader()` (`Utilities.py`) _untyped_
- `getParams(options)` (`Utilities.py`) _untyped_
- `getFooter()` (`Utilities.py`) _untyped_
- `Start(parser, argv, quiet, add_pipe_options, add_extract_options, add_group_dedup_options, add_sam_options, add_dedup_count_sam_options, add_group_sam_options, add_umi_grouping_options, add_s_b_cram_format_options, exclude_output_format, return_parser)` (`Utilities.py`) _untyped_
- `validateExtractOptions(options)` (`Utilities.py`) _untyped_
- `validateSamOptions(options, group)` (`Utilities.py`) _untyped_
- `Stop()` (`Utilities.py`) _untyped_
- `log(loglevel, message)` (`Utilities.py`) _untyped_
- `info(message)` (`Utilities.py`) _untyped_
- `warning(message)` (`Utilities.py`) _untyped_
- `warn(message)` (`Utilities.py`) _untyped_
- `debug(message)` (`Utilities.py`) _untyped_
- `error(message)` (`Utilities.py`) _untyped_
- `critical(message)` (`Utilities.py`) _untyped_
- `getTempFile(dir, shared, suffix)` (`Utilities.py`) _untyped_
- `getTempFilename(dir, shared, suffix)` (`Utilities.py`) _untyped_
- `determine_format(filename, sam, out_format)` (`Utilities.py`) _untyped_
- `output_names_and_formats(outfile, sam, outformat, no_sort_out, tmpdir)` (`Utilities.py`) _untyped_
- `open_input_alignments(in_name, format, options)` (`Utilities.py`) _untyped_
- `open_output_alignments(out_name, format, infile, options)` (`Utilities.py`) _untyped_
- `sort_output(sorted_out_name, infile, format, format_options, reference_filename)` (`Utilities.py`) _untyped_
- `main(argv)` (`count.py`) _untyped_
- `main(argv)` (`count_tab.py`) _untyped_
- `detect_bam_features(bamfile, n_entries)` (`dedup.py`) _untyped_
- `aggregateStatsDF(stats_df)` (`dedup.py`) _untyped_
- `main(argv)` (`dedup.py`) _untyped_
- `main(argv)` (`extract.py`) _untyped_
- `addBarcodesToIdentifier(read, UMI, cell, umi_separator)` (`extract_methods.py`) _untyped_
- `extractSeqAndQuals(seq, quals, umi_bases, cell_bases, discard_bases, retain_umi)` (`extract_methods.py`) _untyped_
- `get_below_threshold(umi_quals, quality_encoding, quality_filter_threshold)` (`extract_methods.py`) _untyped_
- `umi_below_threshold(umi_quals, quality_encoding, quality_filter_threshold)` (`extract_methods.py`) _untyped_
- `mask_umi(umi, umi_quals, quality_encoding, quality_filter_threshold)` (`extract_methods.py`) _untyped_
- `ExtractBarcodes(read, match, extract_umi, extract_cell, discard, retain_umi)` (`extract_methods.py`) _untyped_
- `main(argv)` (`group.py`) _untyped_
- `breadth_first_search(node, adj_list)` (`network.py`) _untyped_
- `recursive_search(node, adj_list)` (`network.py`) _untyped_
- `breadth_first_search_recursive(node, adj_list)` (`network.py`) _untyped_
- `remove_umis(adj_list, cluster, nodes)` (`network.py`) _untyped_
- `get_substr_slices(umi_length, idx_size)` (`network.py`) _untyped_
- `build_substr_idx(umis, umi_length, min_edit)` (`network.py`) _untyped_
- `iter_nearest_neighbours(umis, substr_idx)` (`network.py`) _untyped_
- `chunk_bam(bamfile)` (`prepare_for_em.py`) _untyped_
- `copy_tags(tags, read1, read2)` (`prepare_for_em.py`) _untyped_
- `pick_mate(read, template_dict, mate_key)` (`prepare_for_em.py`) _untyped_
- `main(argv)` (`prepare_for_em.py`) _untyped_
- `get_barcode_read_id(read, cell_barcode, sep)` (`sam_methods.py`) _untyped_
- `get_barcode_tag(read, cell_barcode, umi_tag, cell_tag, umi_tag_split, umi_tag_delim, cell_tag_split, cell_tag_delim)` (`sam_methods.py`) _untyped_
- `get_umi_read_string(read_id, sep)` (`sam_methods.py`) _untyped_
- `get_cell_umi_read_string(read_id, sep)` (`sam_methods.py`) _untyped_
- `get_barcode_umis(read, cell_barcode)` (`sam_methods.py`) _untyped_
- `get_gene_count_tab(infile, bc_getter)` (`sam_methods.py`) _untyped_
- `getMetaContig2contig(bamfile, gene_transcript_map)` (`sam_methods.py`) _untyped_
- `metafetcher(bamfile, metacontig2contig, metatag)` (`sam_methods.py`) _untyped_
- `find_splice(cigar)` (`sam_methods.py`) _untyped_
- `get_read_position(read, soft_clip_threshold)` (`sam_methods.py`) _untyped_
- `fastqIterate(infile, remove_suffix)` (`umi_methods.py`) _untyped_
- `joinedFastqIterate(fastq_iterator1, fastq_iterator2, strict)` (`umi_methods.py`) _untyped_
- `get_average_umi_distance(umis)` (`umi_methods.py`) _untyped_
- `main()` (`umi_tools.py`) _untyped_
- `main(argv)` (`whitelist.py`) _untyped_
- `getKneeEstimateDensity(cell_barcode_counts, expect_cells, cell_number, plotfile_prefix)` (`whitelist_methods.py`) _untyped_
- `getKneeEstimateDistance(cell_barcode_counts, cell_number, plotfile_prefix)` (`whitelist_methods.py`) _untyped_
- `getErrorCorrectMapping(cell_barcodes, whitelist, threshold)` (`whitelist_methods.py`) _untyped_
- `getCellWhitelist(cell_barcode_counts, knee_method, expect_cells, cell_number, error_correct_threshold, plotfile_prefix)` (`whitelist_methods.py`) _untyped_
- `getUserDefinedBarcodes(whitelist_tsv, whitelist_tsv2, getErrorCorrection, deriveErrorCorrection, threshold)` (`whitelist_methods.py`) _untyped_
- `checkError(barcode, whitelist, errors)` (`whitelist_methods.py`) _untyped_
- `errorDetectAboveThreshold(cell_barcode_counts, cell_whitelist, true_to_false_map, errors, resolution_method)` (`whitelist_methods.py`) _untyped_

### Public Classes

- **`DefaultOptions`** (`Utilities.py`), bases: (none)
- **`BetterFormatter`** (`Utilities.py`), bases: `IndentedHelpFormatter`
  - `.__init__(self)`
  - `.format_description(self, description)`
  - `.format_epilog(self, epilog)`
  - `.format_usage(self, usage)`
  - `.format_option(self, option)`
- **`AppendCommaOption`** (`Utilities.py`), bases: `Option`
  - `.convert_value(self, opt, value)`
  - `.take_action(self, action, dest, opt, value, values, parser)`
- **`OptionParser`** (`Utilities.py`), bases: `OptionParser`
  - `.__init__(self)`
- **`OptionGroup`** (`Utilities.py`), bases: `OptionGroup`
- **`MultiLineFormatter`** (`Utilities.py`), bases: `Formatter`
  - `.format(self, record)`
- **`ExtractFilterAndUpdate`** (`extract_methods.py`), bases: (none)
  - `.filterQuality(self, umi_quals)`
  - `.maskQuality(self, umi, umi_quals)`
  - `.filterCellBarcode(self, cell)`
  - `.filterUMIBarcode(self, umi)`
  - `.__init__(self, method, pattern, pattern2, prime3, extract_cell, quality_encoding, quality_filter_threshold, quality_filter_mask, filter_umi_barcode, filter_cell_barcode, retain_umi, either_read, either_read_resolve, umi_separator)`
  - `.getReadCounts(self)`
- **`UMIClusterer`** (`network.py`), bases: (none)
  - `.__init__(self, cluster_method)`
- **`ReadDeduplicator`** (`network.py`), bases: (none)
  - `.__init__(self, options)`
- **`CellClusterer`** (`network.py`), bases: (none)
  - `.__init__(self, cluster_method, dir_threshold, fuzzy_match)`
- **`get_bundles`** (`sam_methods.py`), bases: (none)
  - `.__init__(self, options, only_count_reads, all_reads, return_unmapped, return_read2, metacontig_contig)`
  - `.update_dicts(self, read, pos, key, umi)`
  - `.check_output(self)`
- **`TwoPassPairWriter`** (`sam_methods.py`), bases: (none)
  - `.__init__(self, infile, outfile, tags)`
  - `.write(self, read, unique_id, umi, unmapped)`
  - `.write_mates(self)`
  - `.close(self)`
- **`Record`** (`umi_methods.py`), bases: (none)
  - `.__init__(self, identifier, seq, quals, entry_format)`
  - `.guessFormat(self)`
- **`random_read_generator`** (`umi_methods.py`), bases: (none)
  - `.__init__(self, bamfile, chrom, barcode_getter)`
  - `.refill_random(self)`
  - `.fill(self)`
  - `.getUmis(self, n)`

## Dynamic Behavior Risk

- Total dynamic-behavior markers: 0
- Functions with `*args`/`**kwargs`: 0
- Functions with non-standard decorators: 0
- Untyped public functions: 69

## Native Code

- Cython/native source files (`.pyx`/`.pxd`/`.pxi`): `_dedup_umi.pyx`
- Compiled extension files (`.so`/`.pyd`/`.dylib`): `_dedup_umi.cpython-311-darwin.so`
- `ext_modules`/`ext-modules` declared in build config: no
- Native FFI imports (cffi/ctypes/Cython/pybind11): (none)
- Native submodule usage (highest-priority Rust-port targets):
  - `_dedup_umi` imported in `network.py`: edit_distance
  - `_dedup_umi` imported in `umi_methods.py`: edit_distance
  - `_dedup_umi` imported in `whitelist_methods.py`: edit_distance

## Source Layout

- Modules: 17
- Test files: 0
- Example/notebook files: 0
