//! Command-line option definitions and the `Start`/`Stop` lifecycle
//! (Python source: `umi_tools/Utilities.py`'s `optparse` half).
//!
//! Python builds one `OptionParser` per subcommand out of shared
//! `add_*_options` blocks, then hands every `main()` the resulting untyped
//! `Values` object. This module builds the equivalent `clap` command and
//! translates the parsed result into the concrete structs in
//! `umi_tools_core::options`, so no `clap` type reaches the library.
//!
//! # Deviation: help text layout
//!
//! `optparse`'s `BetterFormatter` wraps and indents help text with its own
//! rules, which `clap` does not reproduce byte for byte. The option *names,
//! arities, defaults and validation* are identical; only the layout differs.
//! See `port_reports/umi_tools/05_design.md` deviation D10.

use std::path::PathBuf;

use clap::{Arg, ArgAction, ArgMatches, Command};

use umi_tools_core::extract::{EitherReadResolve, ExtractMethod};
use umi_tools_core::network::ClusterMethod;
use umi_tools_core::options::{
    BundleOptions, ChimericPairs, DedupOptions, DetectionMethod, ExtractOptions, IoOptions,
    SamValidationInput, UmiMethod, UnmappedReads, UnpairedReads,
};
use umi_tools_core::umi_methods::QualityEncoding;
use umi_tools_core::utilities::{AlignmentFormat, GlobalOptions, LogLevel, LogSink};
use umi_tools_core::{PortError, Result};

/// Splits a comma-separated value into its parts
/// (Python source: `Utilities.AppendCommaOption.convert_value`).
///
/// `optparse` options with `action="append"` treat a comma in the value as a
/// list separator, dropping empty parts — the behaviour the docstring
/// attributes to Galaxy's conventions. An empty string stays a single empty
/// value rather than becoming an empty list, which is the one case a naive
/// `split(',')` gets wrong.
///
/// # Examples
///
/// ```
/// use umi_tools_cli::options::split_append_value;
/// assert_eq!(split_append_value("a,b"), vec!["a", "b"]);
/// assert_eq!(split_append_value("a,,b"), vec!["a", "b"]);
/// assert_eq!(split_append_value(""), vec![""]);
/// assert_eq!(split_append_value("solo"), vec!["solo"]);
/// ```
#[must_use]
pub fn split_append_value(value: &str) -> Vec<String> {
    if !value.contains(',') {
        return vec![value.to_string()];
    }
    value
        .split(',')
        .filter(|v| !v.is_empty())
        .map(str::to_string)
        .collect()
}

/// Which option groups a subcommand's parser should include
/// (Python source: `Utilities.Start`'s `add_*_options` boolean parameters).
#[derive(Debug, Clone, Copy)]
pub struct OptionGroups {
    /// `add_extract_options`: barcode-extraction options.
    pub extract: bool,
    /// `add_group_dedup_options`: `--output-stats` and friends.
    pub group_dedup: bool,
    /// `add_sam_options`: alignment-input options.
    pub sam: bool,
    /// `add_dedup_count_sam_options`: the `dedup`/`count` alignment options.
    pub dedup_count_sam: bool,
    /// `add_group_sam_options`: the `group` alignment options.
    pub group_sam: bool,
    /// `add_umi_grouping_options`: `--method` and `--edit-distance-threshold`.
    pub umi_grouping: bool,
    /// `add_s_b_cram_format_options`: SAM/BAM/CRAM format switches.
    pub sbcram: bool,
    /// `exclude_output_format`: suppress the output-format switches.
    pub exclude_output_format: bool,
}

impl Default for OptionGroups {
    /// Python's `Start` defaults: `add_group_dedup_options=True`,
    /// `add_umi_grouping_options=True`, everything else off.
    fn default() -> Self {
        Self {
            extract: false,
            group_dedup: true,
            sam: false,
            dedup_count_sam: false,
            group_sam: false,
            umi_grouping: true,
            sbcram: false,
            exclude_output_format: false,
        }
    }
}

/// A boolean flag.
///
/// `overrides_with` makes a repeated flag last-wins rather than an error.
/// `optparse`'s `store_true` silently accepts repeats, and the package's own
/// test suite relies on it — one `group` case passes `--out-sam` twice.
fn flag(name: &'static str, long: &'static str, help: &'static str) -> Arg {
    Arg::new(name)
        .long(long)
        .action(ArgAction::SetTrue)
        .overrides_with(name)
        .help(help)
}

/// An option taking one value, last-wins on repeats as `optparse`'s `store`
/// action is.
fn value(name: &'static str, long: &'static str, help: &'static str) -> Arg {
    Arg::new(name)
        .long(long)
        .num_args(1)
        .overrides_with(name)
        .help(help)
}

/// Adds the options every subcommand shares
/// (Python source: the `add_pipe_options` block of `Utilities.Start`).
fn add_pipe_options(command: Command) -> Command {
    command
        .arg(value("stdin", "stdin", "file to read stdin from").short('I'))
        .arg(value("stdout", "stdout", "file where output is to go").short('S'))
        .arg(value("stdlog", "log", "file with logging information").short('L'))
        .arg(value("stderr", "error", "file with error information").short('E'))
        .arg(flag(
            "log2stderr",
            "log2stderr",
            "send logging information to stderr",
        ))
        // `-v/--verbose` takes an integer in the original, not a bare flag.
        .arg(
            value(
                "loglevel",
                "verbose",
                "loglevel. The higher, the more output",
            )
            .short('v'),
        )
        .arg(value(
            "timeit_file",
            "timeit",
            "store timing information in file",
        ))
        .arg(value("timeit_name", "timeit-name", "name in timing file"))
        .arg(flag(
            "timeit_header",
            "timeit-header",
            "add header for timing information",
        ))
        .arg(value(
            "random_seed",
            "random-seed",
            "random seed to initialize number generator with",
        ))
        .arg(value("tmpdir", "temp-dir", "directory for temporary files"))
        .arg(value(
            "compresslevel",
            "compresslevel",
            "compression level for gzip output",
        ))
        // `Utilities.Start` declares this with `action='help'`, so it is a
        // second spelling of `--help`; `crate::wants_help` intercepts both
        // before parsing, and it is declared here so the option surface
        // matches the help text.
        .arg(flag(
            "help_extended",
            "help-extended",
            "output full documentation",
        ))
}

/// Adds the barcode-extraction options
/// (Python source: the `add_extract_options` block).
fn add_extract_options(command: Command) -> Command {
    command
        .arg(value(
            "extract_method",
            "extract-method",
            "how to extract the umi +/- cell barcodes: string or regex",
        ))
        .arg(value("pattern", "bc-pattern", "barcode pattern").short('p'))
        .arg(value(
            "pattern2",
            "bc-pattern2",
            "barcode pattern for read 2",
        ))
        .arg(flag("prime3", "3prime", "barcode is on 3' end of read"))
        .arg(value("read2_in", "read2-in", "file name for read pairs"))
        .arg(value(
            "filtered_out",
            "filtered-out",
            "write out reads not matching regex pattern or cell barcode whitelist to this file",
        ))
        .arg(value(
            "filtered_out2",
            "filtered-out2",
            "write out read 2s not matching regex pattern or cell barcode whitelist to this file",
        ))
        .arg(flag(
            "ignore_suffix",
            "ignore-read-pair-suffixes",
            "ignore SRR/ENA read name suffixes '/1' and '/2'",
        ))
        .arg(value(
            "quality_filter_threshold",
            "quality-filter-threshold",
            "remove reads where any UMI base quality score falls below this threshold",
        ))
        .arg(value(
            "quality_filter_mask",
            "quality-filter-mask",
            "if a UMI base has a quality below this threshold, replace the base with 'N'",
        ))
        .arg(value(
            "quality_encoding",
            "quality-encoding",
            "quality score encoding: phred33, phred64 or solexa",
        ))
}

/// Adds the UMI grouping options
/// (Python source: the `add_umi_grouping_options` block).
fn add_umi_grouping_options(command: Command) -> Command {
    command
        .arg(value(
            "method",
            "method",
            "method to use for umi grouping: unique, percentile, cluster, adjacency or directional",
        ))
        .arg(value(
            "threshold",
            "edit-distance-threshold",
            "edit distance threshold at which to join two UMIs when grouping",
        ))
}

/// Adds the alignment-input options shared by `group`, `dedup` and `count`
/// (Python source: the `add_sam_options` block).
fn add_sam_options(command: Command) -> Command {
    command
        .arg(value(
            "get_umi_method",
            "extract-umi-method",
            "how is the read UMI encoded: read_id, tag or umis",
        ))
        .arg(value(
            "umi_sep",
            "umi-separator",
            "separator between read id and UMI",
        ))
        .arg(value("umi_tag", "umi-tag", "tag containing umi"))
        .arg(value(
            "umi_tag_split",
            "umi-tag-split",
            "split UMI in tag and take the first element",
        ))
        .arg(value(
            "umi_tag_delim",
            "umi-tag-delimiter",
            "concatenate UMI in tag separated by delimiter",
        ))
        .arg(value("cell_tag", "cell-tag", "tag containing cell barcode"))
        .arg(value(
            "cell_tag_split",
            "cell-tag-split",
            "split cell barcode in tag and take the first element",
        ))
        .arg(value(
            "cell_tag_delim",
            "cell-tag-delimiter",
            "concatenate cell barcode in tag separated by delimiter",
        ))
        .arg(flag("paired", "paired", "BAM is paired end"))
        .arg(flag(
            "per_cell",
            "per-cell",
            "group/dedup/count per cell barcode",
        ))
        .arg(flag("per_gene", "per-gene", "group/dedup/count per gene"))
        .arg(flag(
            "per_contig",
            "per-contig",
            "group/dedup/count per contig",
        ))
        .arg(value(
            "gene_tag",
            "gene-tag",
            "tag which contains the gene assignment",
        ))
        .arg(value(
            "assigned_tag",
            "assigned-status-tag",
            "tag which contains the assignment status",
        ))
        .arg(value(
            "skip_regex",
            "skip-tags-regex",
            "skip any reads where the assignment tag matches this regex",
        ))
        .arg(value(
            "gene_transcript_map",
            "gene-transcript-map",
            "file mapping transcripts to genes",
        ))
        .arg(value(
            "subset",
            "subset",
            "randomly select this fraction of reads",
        ))
        .arg(value("chrom", "chrom", "restrict to one chromosome"))
        .arg(value(
            "mapping_quality",
            "mapping-quality",
            "minimum mapping quality for a read to be retained",
        ))
        .arg(flag(
            "ignore_umi",
            "ignore-umi",
            "ignore UMI and dedup only on position",
        ))
        .arg(flag(
            "ignore_tlen",
            "ignore-tlen",
            "ignore the template length when grouping",
        ))
        .arg(flag(
            "read_length",
            "read-length",
            "use read length in addition to position",
        ))
        .arg(flag(
            "spliced",
            "spliced-is-unique",
            "spliced reads are treated as different to unspliced reads",
        ))
        .arg(value(
            "soft_clip_threshold",
            "soft-clip-threshold",
            "number of bases clipped from 5' end before a read is deemed spliced",
        ))
        .arg(value(
            "detection_method",
            "multimapping-detection-method",
            "tag to use to identify multimapping reads: NH, X0 or XT",
        ))
        // Both spellings set the same destination in the original.
        .arg(
            flag(
                "whole_contig",
                "whole-contig",
                "consider all alignments to a single contig together",
            )
            .alias("buffer-whole-contig"),
        )
        .arg(value(
            "unmapped_reads",
            "unmapped-reads",
            "how to handle unmapped reads: discard, use or output",
        ))
        .arg(value(
            "chimeric_pairs",
            "chimeric-pairs",
            "how to handle chimeric read pairs: discard, use or output",
        ))
        .arg(value(
            "unpaired_reads",
            "unpaired-reads",
            "how to handle unpaired reads: discard, use or output",
        ))
        .arg(flag(
            "output_unmapped",
            "output-unmapped",
            "deprecated; use --unmapped-reads=output",
        ))
}

/// Adds the SAM/BAM/CRAM format switches
/// (Python source: the `add_s_b_cram_format_options` block).
fn add_sbcram_options(command: Command, exclude_output_format: bool) -> Command {
    let command = command
        .arg(flag("in_sam", "in-sam", "input is SAM").short('i'))
        .arg(value(
            "in_format",
            "in-format",
            "input format: sam, bam or cram",
        ))
        .arg(value(
            "input_options",
            "input-options",
            "comma-separated format options for the input file",
        ))
        .arg(value(
            "reference_filename",
            "reference-filename",
            "reference sequence FASTA, required for CRAM",
        ));

    if exclude_output_format {
        return command;
    }
    command
        .arg(flag("out_sam", "out-sam", "output is SAM").short('o'))
        .arg(value(
            "out_format",
            "out-format",
            "output format: sam, bam or cram",
        ))
        .arg(value(
            "output_options",
            "output-options",
            "comma-separated format options for the output file",
        ))
        .arg(flag(
            "no_sort_output",
            "no-sort-output",
            "do not sort the output",
        ))
}

/// Builds the `clap` command for a subcommand
/// (Python source: `Utilities.OptionParser` plus the `add_*_options` blocks
/// `Utilities.Start` selects).
#[must_use]
pub fn build_parser(name: &'static str, about: &'static str, groups: OptionGroups) -> Command {
    let mut command = Command::new(name)
        .about(about)
        .disable_help_flag(false)
        .disable_version_flag(true)
        // `optparse` accepts any unambiguous prefix of a long option, and the
        // package's own test suite relies on it — several cases pass
        // `--reference-file` for `--reference-filename`. `clap` needs to be
        // told to do the same.
        .infer_long_args(true)
        .override_usage(format!("umi_tools {name} [OPTIONS]"));

    command = add_pipe_options(command);
    if groups.extract {
        command = add_extract_options(command);
    }
    if groups.umi_grouping {
        command = add_umi_grouping_options(command);
    }
    if groups.sam || groups.dedup_count_sam || groups.group_sam {
        command = add_sam_options(command);
    }
    if groups.sbcram {
        command = add_sbcram_options(command, groups.exclude_output_format);
    }
    if groups.group_dedup {
        command = command
            .arg(value(
                "umi_whitelist",
                "umi-whitelist",
                "whitelist of accepted UMIs",
            ))
            .arg(value(
                "umi_whitelist_paired",
                "umi-whitelist-paired",
                "whitelist of accepted UMIs for read 2",
            ))
            .arg(flag(
                "filter_umi",
                "filter-umi",
                "filter UMIs against the whitelist",
            ));
    }
    command
}

/// Reads a string option that may not be registered on this subcommand's parser.
///
/// `clap`'s `get_one` panics for an id the command never defined, and the
/// shared translation functions here are called from parsers built with
/// different option groups — so every access has to be the fallible form.
fn opt_str(matches: &ArgMatches, name: &str) -> Option<String> {
    matches.try_get_one::<String>(name).ok().flatten().cloned()
}

/// Reads a boolean flag that may not be registered on this subcommand's parser.
fn opt_flag(matches: &ArgMatches, name: &str) -> bool {
    matches
        .try_get_one::<bool>(name)
        .ok()
        .flatten()
        .copied()
        .unwrap_or(false)
}

fn parse_opt<T: std::str::FromStr<Err = PortError>>(
    matches: &ArgMatches,
    name: &str,
) -> Result<Option<T>> {
    opt_str(matches, name).map(|v| v.parse::<T>()).transpose()
}

fn parse_num<T: std::str::FromStr>(
    matches: &ArgMatches,
    name: &str,
    label: &str,
) -> Result<Option<T>>
where
    T::Err: std::fmt::Display,
{
    opt_str(matches, name)
        .map(|v| {
            v.parse::<T>()
                .map_err(|e| PortError::value(format!("{label}: {e}")))
        })
        .transpose()
}

fn two_char_tag(matches: &ArgMatches, name: &str) -> Result<Option<[u8; 2]>> {
    match opt_str(matches, name) {
        None => Ok(None),
        Some(v) => {
            let bytes = v.as_bytes();
            if bytes.len() != 2 {
                return Err(PortError::value(format!(
                    "SAM tags are two characters; got '{v}'"
                )));
            }
            Ok(Some([bytes[0], bytes[1]]))
        }
    }
}

/// Translates parsed arguments into [`IoOptions`].
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Value`] for an unparseable numeric or format value.
pub fn io_options(matches: &ArgMatches) -> Result<IoOptions> {
    let format = |name: &str| -> Result<Option<AlignmentFormat>> {
        opt_str(matches, name)
            .map(|v| AlignmentFormat::parse(&v))
            .transpose()
    };
    Ok(IoOptions {
        stdin: opt_str(matches, "stdin").map(PathBuf::from),
        stdout: opt_str(matches, "stdout").map(PathBuf::from),
        stdlog: opt_str(matches, "stdlog").map(PathBuf::from),
        log2stderr: opt_flag(matches, "log2stderr"),
        loglevel: parse_num::<i32>(matches, "loglevel", "--verbose")?.unwrap_or(1),
        tmpdir: opt_str(matches, "tmpdir").map(PathBuf::from),
        in_sam: opt_flag(matches, "in_sam"),
        out_sam: opt_flag(matches, "out_sam"),
        in_format: format("in_format")?,
        out_format: format("out_format")?,
        no_sort_output: opt_flag(matches, "no_sort_output"),
        input_options: opt_str(matches, "input_options"),
        output_options: opt_str(matches, "output_options"),
        reference_filename: opt_str(matches, "reference_filename").map(PathBuf::from),
        chrom: opt_str(matches, "chrom"),
        compresslevel: parse_num::<u32>(matches, "compresslevel", "--compresslevel")?.unwrap_or(6),
    })
}

/// Translates parsed arguments into [`BundleOptions`].
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Value`] for an unparseable value or an invalid tag.
pub fn bundle_options(matches: &ArgMatches) -> Result<BundleOptions> {
    let get_flag = |name: &str| opt_flag(matches, name);
    Ok(BundleOptions {
        paired: get_flag("paired"),
        per_cell: get_flag("per_cell"),
        per_gene: get_flag("per_gene"),
        per_contig: get_flag("per_contig"),
        gene_tag: two_char_tag(matches, "gene_tag")?,
        assigned_tag: two_char_tag(matches, "assigned_tag")?,
        skip_regex: opt_str(matches, "skip_regex")
            .map(|p| {
                regex::Regex::new(&p)
                    .map_err(|_| PortError::value(format!("skip-regex '{p}' is not a valid regex")))
            })
            .transpose()?,
        spliced: get_flag("spliced"),
        ignore_tlen: get_flag("ignore_tlen"),
        ignore_umi: get_flag("ignore_umi"),
        read_length: get_flag("read_length"),
        soft_clip_threshold: parse_num::<usize>(
            matches,
            "soft_clip_threshold",
            "--soft-clip-threshold",
        )?
        .unwrap_or(4),
        subset: parse_num::<f64>(matches, "subset", "--subset")?,
        mapping_quality: parse_num::<u32>(matches, "mapping_quality", "--mapping-quality")?,
        detection_method: parse_opt::<DetectionMethod>(matches, "detection_method")?,
        whole_contig: get_flag("whole_contig"),
        unmapped_reads: parse_opt::<UnmappedReads>(matches, "unmapped_reads")?.unwrap_or_default(),
        unpaired_reads: parse_opt::<UnpairedReads>(matches, "unpaired_reads")?.unwrap_or_default(),
        chimeric_pairs: parse_opt::<ChimericPairs>(matches, "chimeric_pairs")?.unwrap_or_default(),
        random_seed: parse_num::<u64>(matches, "random_seed", "--random-seed")?,
    })
}

/// Translates parsed arguments into [`DedupOptions`].
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::UnknownClusterMethod`] for an unknown `--method`, or
/// [`umi_tools_core::PortError::Value`] for an unparseable threshold.
pub fn dedup_options(matches: &ArgMatches) -> Result<DedupOptions> {
    Ok(DedupOptions {
        method: parse_opt::<ClusterMethod>(matches, "method")?
            .unwrap_or(ClusterMethod::Directional),
        threshold: parse_num::<u32>(matches, "threshold", "--edit-distance-threshold")?
            .unwrap_or(1),
        filter_umi: opt_flag(matches, "filter_umi"),
        umi_whitelist: opt_str(matches, "umi_whitelist").map(PathBuf::from),
        umi_whitelist_paired: opt_str(matches, "umi_whitelist_paired").map(PathBuf::from),
    })
}

/// Translates parsed arguments into [`ExtractOptions`].
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Value`] for an unknown `--extract-method`.
pub fn extract_options(matches: &ArgMatches) -> Result<ExtractOptions> {
    let method = opt_str(matches, "extract_method").unwrap_or_else(|| "string".to_string());
    // Reject an unknown method here rather than deep inside extraction.
    method.parse::<ExtractMethod>()?;
    Ok(ExtractOptions {
        extract_method: method,
        pattern: opt_str(matches, "pattern"),
        pattern2: opt_str(matches, "pattern2"),
        prime3: opt_flag(matches, "prime3"),
        read2_in: opt_str(matches, "read2_in").map(PathBuf::from),
        read2_out: opt_str(matches, "read2_out").map(PathBuf::from),
        read2_only: opt_flag(matches, "read2_only"),
        filtered_out: opt_str(matches, "filtered_out").map(PathBuf::from),
        filtered_out2: opt_str(matches, "filtered_out2").map(PathBuf::from),
    })
}

/// Builds the input [`SamValidationInput`] from parsed arguments.
///
/// `argv` is inspected for which paired-only flags the user actually typed,
/// because Python's validator tests `" ".join(sys.argv)` for them — a flag left
/// at its default is treated differently from one explicitly set to that
/// default.
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Value`] for an unknown `--extract-umi-method`.
pub fn sam_validation_input(matches: &ArgMatches, argv: &[String]) -> Result<SamValidationInput> {
    let get_flag = |name: &str| opt_flag(matches, name);
    Ok(SamValidationInput {
        per_gene: get_flag("per_gene"),
        per_contig: get_flag("per_contig"),
        gene_tag: opt_str(matches, "gene_tag"),
        assigned_tag: opt_str(matches, "assigned_tag"),
        gene_transcript_map: opt_str(matches, "gene_transcript_map").map(PathBuf::from),
        get_umi_method: parse_opt::<UmiMethod>(matches, "get_umi_method")?.unwrap_or_default(),
        umi_tag: opt_str(matches, "umi_tag"),
        cell_tag: opt_str(matches, "cell_tag"),
        per_cell: get_flag("per_cell"),
        skip_regex: opt_str(matches, "skip_regex"),
        paired: get_flag("paired"),
        unmapped_reads: parse_opt::<UnmappedReads>(matches, "unmapped_reads")?.unwrap_or_default(),
        unpaired_reads: parse_opt::<UnpairedReads>(matches, "unpaired_reads")?.unwrap_or_default(),
        chimeric_pairs: parse_opt::<ChimericPairs>(matches, "chimeric_pairs")?.unwrap_or_default(),
        output_unmapped: get_flag("output_unmapped"),
        flags_present: argv
            .iter()
            .filter(|a| a.starts_with("--"))
            .map(|a| a.split('=').next().unwrap_or(a).to_string())
            .collect(),
    })
}

/// The quality encoding named on the command line, if any.
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Value`] for a name that is not one of the three known
/// encodings.
pub fn quality_encoding(matches: &ArgMatches) -> Result<Option<QualityEncoding>> {
    parse_opt::<QualityEncoding>(matches, "quality_encoding")
}

/// The either-read resolution named on the command line, if any.
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Value`] for anything but `discard` or `quality`.
pub fn either_read_resolve(matches: &ArgMatches) -> Result<Option<EitherReadResolve>> {
    parse_opt::<EitherReadResolve>(matches, "either_read_resolve")
}

/// Configures logging and writes the run header
/// (Python source: `Utilities.Start`).
///
/// Returns the run identifier, which [`stop`] needs for the footer.
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Io`] if the log file cannot be opened.
pub fn start(io: &IoOptions, argv: &[String]) -> Result<String> {
    // Python chooses the `#`-prefixed log format only when the log and the
    // data share a stream, so the log does not corrupt the data.
    let sharing_stream = io.stdlog.is_none() && !io.log2stderr && io.stdout.is_none();
    let sink = match (&io.stdlog, io.log2stderr) {
        (Some(path), _) => LogSink::File(path.clone()),
        (None, true) => LogSink::Stderr,
        (None, false) => LogSink::Stdout,
    };
    umi_tools_core::utilities::set_global_options(GlobalOptions {
        loglevel: LogLevel::from_verbosity(io.loglevel),
        sink,
        comment_prefix: sharing_stream,
        compresslevel: io.compresslevel,
    });

    let run_id = run_identifier();
    if io.loglevel >= 1 {
        let host = hostname();
        let header = umi_tools_core::utilities::get_header(
            argv,
            umi_tools_core::VERSION,
            &run_id,
            &host,
            &uname(),
        );
        write_log_lines(io, &header)?;
    }
    Ok(run_id)
}

/// Writes the run footer (Python source: `Utilities.Stop`).
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Io`] if the log file cannot be written.
pub fn stop(io: &IoOptions, run_id: &str, elapsed_secs: u64) -> Result<()> {
    if io.loglevel >= 1 {
        let footer = umi_tools_core::utilities::get_footer(elapsed_secs, run_id, [0.0; 4]);
        write_log_lines(io, &footer)?;
    }
    Ok(())
}

fn write_log_lines(io: &IoOptions, text: &str) -> Result<()> {
    use std::io::Write as _;
    match (&io.stdlog, io.log2stderr) {
        (Some(path), _) => {
            let mut file = umi_tools_core::utilities::append_file(path)?;
            writeln!(file, "{text}").map_err(|e| PortError::io(path, e))
        }
        (None, true) => {
            let mut out = std::io::stderr().lock();
            writeln!(out, "{text}").map_err(|e| PortError::io(std::path::Path::new("<stderr>"), e))
        }
        (None, false) => {
            let mut out = std::io::stdout().lock();
            writeln!(out, "{text}").map_err(|e| PortError::io(std::path::Path::new("<stdout>"), e))
        }
    }
}

/// A per-run identifier, standing in for Python's `uuid.uuid4()`.
///
/// It only has to be unique per run and stable within one, so it is derived
/// from the process id and the wall clock rather than pulling in a UUID
/// dependency for a log field.
fn run_identifier() -> String {
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default();
    format!("{:08x}-{:08x}", std::process::id(), now.as_nanos() as u32)
}

fn hostname() -> String {
    std::env::var("HOSTNAME").unwrap_or_else(|_| "localhost".to_string())
}

fn uname() -> String {
    format!("{} {}", std::env::consts::OS, std::env::consts::ARCH)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(args: &[&str], groups: OptionGroups) -> ArgMatches {
        build_parser("test", "test command", groups)
            .no_binary_name(true)
            .try_get_matches_from(args)
            .expect("arguments parse")
    }

    #[test]
    fn comma_splitting_matches_optparse() {
        assert_eq!(split_append_value("a,b"), vec!["a", "b"]);
        assert_eq!(split_append_value("a,,b"), vec!["a", "b"]);
        assert_eq!(split_append_value(",a,"), vec!["a"]);
        assert_eq!(split_append_value(""), vec![""]);
        assert_eq!(split_append_value("solo"), vec!["solo"]);
    }

    #[test]
    fn io_defaults_match_python() {
        let m = parse(&[], OptionGroups::default());
        let io = io_options(&m).unwrap();
        assert_eq!(io.loglevel, 1);
        assert_eq!(io.compresslevel, 6);
        assert!(!io.log2stderr);
        assert!(io.stdin.is_none());
    }

    /// `-v/--verbose` takes an integer in the original, not a bare flag.
    #[test]
    fn verbose_takes_a_numeric_level() {
        let m = parse(&["--verbose=2"], OptionGroups::default());
        assert_eq!(io_options(&m).unwrap().loglevel, 2);
        let m = parse(&["-v", "0"], OptionGroups::default());
        assert_eq!(io_options(&m).unwrap().loglevel, 0);
    }

    /// `optparse` accepts unambiguous long-option prefixes, and the package's
    /// own test suite passes `--reference-file` for `--reference-filename`.
    #[test]
    fn unambiguous_long_option_prefixes_are_accepted() {
        let groups = OptionGroups {
            sbcram: true,
            ..OptionGroups::default()
        };
        let m = parse(&["--reference-file=ref.fa"], groups);
        assert_eq!(
            io_options(&m).unwrap().reference_filename.as_deref(),
            Some(std::path::Path::new("ref.fa"))
        );
    }

    /// Both spellings of the buffering flag set the same option, as in the
    /// original.
    #[test]
    fn whole_contig_has_both_spellings() {
        let groups = OptionGroups {
            sam: true,
            ..OptionGroups::default()
        };
        assert!(
            bundle_options(&parse(&["--whole-contig"], groups))
                .unwrap()
                .whole_contig
        );
        assert!(
            bundle_options(&parse(&["--buffer-whole-contig"], groups))
                .unwrap()
                .whole_contig
        );
    }

    #[test]
    fn bundle_defaults_match_python() {
        let groups = OptionGroups {
            sam: true,
            ..OptionGroups::default()
        };
        let m = parse(&[], groups);
        let b = bundle_options(&m).unwrap();
        assert_eq!(b.soft_clip_threshold, 4);
        assert!(!b.paired);
        assert_eq!(b.unmapped_reads, UnmappedReads::Discard);
    }

    #[test]
    fn bundle_options_are_read_from_the_command_line() {
        let groups = OptionGroups {
            sam: true,
            ..OptionGroups::default()
        };
        let m = parse(
            &[
                "--paired",
                "--per-cell",
                "--mapping-quality=20",
                "--soft-clip-threshold=7",
                "--unmapped-reads=use",
                "--multimapping-detection-method=NH",
                "--gene-tag=XT",
            ],
            groups,
        );
        let b = bundle_options(&m).unwrap();
        assert!(b.paired && b.per_cell);
        assert_eq!(b.mapping_quality, Some(20));
        assert_eq!(b.soft_clip_threshold, 7);
        assert_eq!(b.unmapped_reads, UnmappedReads::Use);
        assert_eq!(b.detection_method, Some(DetectionMethod::Nh));
        assert_eq!(b.gene_tag, Some(*b"XT"));
    }

    #[test]
    fn dedup_defaults_and_overrides() {
        let m = parse(&[], OptionGroups::default());
        let d = dedup_options(&m).unwrap();
        assert_eq!(d.method, ClusterMethod::Directional);
        assert_eq!(d.threshold, 1);

        let m = parse(
            &["--method=adjacency", "--edit-distance-threshold=2"],
            OptionGroups::default(),
        );
        let d = dedup_options(&m).unwrap();
        assert_eq!(d.method, ClusterMethod::Adjacency);
        assert_eq!(d.threshold, 2);
    }

    #[test]
    fn an_unknown_method_is_rejected() {
        let m = parse(&["--method=bogus"], OptionGroups::default());
        assert!(matches!(
            dedup_options(&m),
            Err(PortError::UnknownClusterMethod(_))
        ));
    }

    #[test]
    fn a_malformed_tag_is_rejected() {
        let groups = OptionGroups {
            sam: true,
            ..OptionGroups::default()
        };
        let m = parse(&["--gene-tag=TOOLONG"], groups);
        assert!(bundle_options(&m).is_err());
    }

    #[test]
    fn an_invalid_skip_regex_is_rejected() {
        let groups = OptionGroups {
            sam: true,
            ..OptionGroups::default()
        };
        let m = parse(&["--skip-tags-regex=(unclosed"], groups);
        assert!(bundle_options(&m).is_err());
    }

    #[test]
    fn extract_options_are_read_from_the_command_line() {
        let groups = OptionGroups {
            extract: true,
            ..OptionGroups::default()
        };
        let m = parse(
            &[
                "--extract-method=regex",
                "--bc-pattern=^(?P<umi_1>.{3})",
                "--3prime",
            ],
            groups,
        );
        let e = extract_options(&m).unwrap();
        assert_eq!(e.extract_method, "regex");
        assert_eq!(e.pattern.as_deref(), Some("^(?P<umi_1>.{3})"));
        assert!(e.prime3);
    }

    #[test]
    fn an_unknown_extract_method_is_rejected() {
        let groups = OptionGroups {
            extract: true,
            ..OptionGroups::default()
        };
        let m = parse(&["--extract-method=magic"], groups);
        assert!(extract_options(&m).is_err());
    }

    #[test]
    fn format_switches_parse() {
        let groups = OptionGroups {
            sbcram: true,
            ..OptionGroups::default()
        };
        let m = parse(
            &["--in-sam", "--out-format=cram", "--no-sort-output"],
            groups,
        );
        let io = io_options(&m).unwrap();
        assert!(io.in_sam);
        assert_eq!(io.out_format, Some(AlignmentFormat::Cram));
        assert!(io.no_sort_output);
    }

    #[test]
    fn output_format_options_can_be_excluded() {
        let groups = OptionGroups {
            sbcram: true,
            exclude_output_format: true,
            ..OptionGroups::default()
        };
        let command = build_parser("test", "t", groups).no_binary_name(true);
        assert!(command
            .clone()
            .try_get_matches_from(["--out-format=bam"])
            .is_err());
        assert!(command.try_get_matches_from(["--in-format=bam"]).is_ok());
    }

    #[test]
    fn typed_flags_are_recorded_for_validation() {
        let groups = OptionGroups {
            sam: true,
            ..OptionGroups::default()
        };
        let m = parse(&["--paired", "--ignore-tlen"], groups);
        let argv = vec!["--paired".to_string(), "--ignore-tlen".to_string()];
        let input = sam_validation_input(&m, &argv).unwrap();
        assert!(input.flags_present.contains(&"--ignore-tlen".to_string()));
        assert!(!input
            .flags_present
            .contains(&"--chimeric-pairs".to_string()));
    }

    #[test]
    fn typed_flags_strip_their_values() {
        let groups = OptionGroups {
            sam: true,
            ..OptionGroups::default()
        };
        let m = parse(&["--unmapped-reads=discard"], groups);
        let argv = vec!["--unmapped-reads=discard".to_string()];
        let input = sam_validation_input(&m, &argv).unwrap();
        assert!(input
            .flags_present
            .contains(&"--unmapped-reads".to_string()));
    }
}
