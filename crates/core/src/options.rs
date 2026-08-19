//! Concrete option structs, replacing Python's `optparse.Values` god-object
//! (Python source: `umi_tools/Utilities.py`'s `Start` and the option
//! validation it feeds).
//!
//! `sam_methods.get_bundles`, `network.ReadDeduplicator` and all seven
//! `main()` functions take a single untyped `options` object and read some
//! fifty attributes off it, defined across four separate `add_*_options`
//! blocks. There is no class, no schema, and no way to discover the valid
//! attribute set short of reading `Start`'s thousand-line body together with
//! every subcommand's `parser.set_defaults` call.
//!
//! This module resolves that into three structs, split by who actually reads
//! them, so the library never sees a `clap` type and the CLI does the
//! translation. See `port_reports/umi_tools/05_design.md` § 2.2.

use regex::Regex;

use crate::error::{PortError, Result};
use crate::network::ClusterMethod;
use crate::utilities::AlignmentFormat;

/// Which SAM tag identifies a multimapping read
/// (Python source: the `--multimapping-detection-method` choice).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DetectionMethod {
    /// The `NH` tag: number of reported alignments.
    Nh,
    /// The `X0` tag: number of best hits (BWA).
    X0,
    /// The `XT` tag: `U` for a unique hit (BWA).
    Xt,
}

impl std::str::FromStr for DetectionMethod {
    type Err = PortError;

    /// # Errors
    ///
    /// Returns [`PortError::Value`] for anything but `NH`, `X0` or `XT`.
    fn from_str(s: &str) -> Result<Self> {
        match s {
            "NH" => Ok(DetectionMethod::Nh),
            "X0" => Ok(DetectionMethod::X0),
            "XT" => Ok(DetectionMethod::Xt),
            other => Err(PortError::value(format!(
                "unknown multimapping detection method: {other}"
            ))),
        }
    }
}

/// What to do with unmapped reads (Python source: `--unmapped-reads`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum UnmappedReads {
    /// Drop them.
    #[default]
    Discard,
    /// Write them out without grouping (`group` only).
    Output,
    /// Group the pair using read 1's coordinates alone.
    Use,
}

/// What to do with a read whose mate is unpaired (Python source:
/// `--unpaired-reads`).
///
/// The default is `Use`, not `Discard`: `Utilities.Start` declares
/// `default="use"` for this option in both the dedup/count and group blocks.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum UnpairedReads {
    /// Drop them.
    Discard,
    /// Write them out without grouping (`group` only).
    Output,
    /// Group them using read 1's coordinates alone.
    #[default]
    Use,
}

/// What to do with a read pair mapped to two contigs (Python source:
/// `--chimeric-pairs`).
///
/// The default is `Use`, not `Discard`: `Utilities.Start` declares
/// `default="use"` for this option in both the dedup/count and group blocks.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ChimericPairs {
    /// Drop them.
    Discard,
    /// Write them out without grouping (`group` only).
    Output,
    /// Group them using read 1's coordinates alone.
    #[default]
    Use,
}

macro_rules! impl_three_way_from_str {
    ($ty:ty, $flag:literal, $discard:expr, $output:expr, $use_:expr) => {
        impl std::str::FromStr for $ty {
            type Err = PortError;

            /// # Errors
            ///
            /// Returns [`PortError::Value`] for anything but `discard`,
            /// `output` or `use`.
            fn from_str(s: &str) -> Result<Self> {
                match s {
                    "discard" => Ok($discard),
                    "output" => Ok($output),
                    "use" => Ok($use_),
                    other => Err(PortError::value(format!(
                        concat!($flag, " must be discard, output or use, not {}"),
                        other
                    ))),
                }
            }
        }
    };
}

impl_three_way_from_str!(
    UnmappedReads,
    "--unmapped-reads",
    UnmappedReads::Discard,
    UnmappedReads::Output,
    UnmappedReads::Use
);
impl_three_way_from_str!(
    UnpairedReads,
    "--unpaired-reads",
    UnpairedReads::Discard,
    UnpairedReads::Output,
    UnpairedReads::Use
);
impl_three_way_from_str!(
    ChimericPairs,
    "--chimeric-pairs",
    ChimericPairs::Discard,
    ChimericPairs::Output,
    ChimericPairs::Use
);

/// How the UMI is stored in an alignment record
/// (Python source: the `--extract-umi-method` choice).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum UmiMethod {
    /// Appended to the read name, as `umi_tools extract` writes it.
    #[default]
    ReadId,
    /// In a SAM tag.
    Tag,
    /// In the `umis` toolkit's `:UMI_<seq>:` read-name convention.
    Umis,
}

impl std::str::FromStr for UmiMethod {
    type Err = PortError;

    /// # Errors
    ///
    /// Returns [`PortError::Value`] with Python's message
    /// `"Unknown UMI extraction method"`.
    fn from_str(s: &str) -> Result<Self> {
        match s {
            "read_id" => Ok(UmiMethod::ReadId),
            "tag" => Ok(UmiMethod::Tag),
            "umis" => Ok(UmiMethod::Umis),
            _ => Err(PortError::value("Unknown UMI extraction method")),
        }
    }
}

/// Everything [`crate::bundles::GetBundles`] reads off Python's `options`
/// object.
///
/// [`Default`] gives Python's defaults, so a caller sets only what it means to
/// change.
#[derive(Debug, Clone)]
pub struct BundleOptions {
    /// Input is paired-end (`--paired`).
    pub paired: bool,
    /// Group per cell barcode (`--per-cell`).
    pub per_cell: bool,
    /// Group per gene rather than per coordinate (`--per-gene`).
    pub per_gene: bool,
    /// Treat each contig as a gene (`--per-contig`).
    pub per_contig: bool,
    /// The SAM tag holding the gene assignment (`--gene-tag`).
    pub gene_tag: Option<[u8; 2]>,
    /// The SAM tag holding the assignment status (`--assigned-status-tag`);
    /// defaults to `gene_tag`.
    pub assigned_tag: Option<[u8; 2]>,
    /// Skip reads whose assignment tag matches this pattern (`--skip-tags-regex`).
    pub skip_regex: Option<Regex>,
    /// Treat a spliced read as distinct from an unspliced one at the same
    /// position (`--spliced-is-unique`).
    pub spliced: bool,
    /// Ignore the template length when grouping (`--ignore-tlen`).
    pub ignore_tlen: bool,
    /// Ignore the UMI entirely (`--ignore-umi`).
    pub ignore_umi: bool,
    /// Include the read length in the grouping key (`--read-length`).
    pub read_length: bool,
    /// Soft clips longer than this count as splicing (`--soft-clip-threshold`).
    pub soft_clip_threshold: usize,
    /// Randomly retain this fraction of reads (`--subset`).
    pub subset: Option<f64>,
    /// Drop reads below this mapping quality (`--mapping-quality`).
    pub mapping_quality: Option<u32>,
    /// How multimapping reads are identified
    /// (`--multimapping-detection-method`).
    pub detection_method: Option<DetectionMethod>,
    /// Buffer a whole contig before flushing (`--whole-contig`).
    pub whole_contig: bool,
    /// What to do with unmapped reads.
    pub unmapped_reads: UnmappedReads,
    /// What to do with unpaired reads.
    pub unpaired_reads: UnpairedReads,
    /// What to do with chimeric read pairs.
    pub chimeric_pairs: ChimericPairs,
    /// Seed for the tie-breaking sampler (`--random-seed`).
    pub random_seed: Option<u64>,
}

impl Default for BundleOptions {
    /// Python's defaults, from `Utilities.Start`'s `add_option` calls:
    /// `soft_clip_threshold=4`, `chimeric_pairs="use"`, `unpaired_reads="use"`,
    /// `unmapped_reads="discard"`, everything else off or absent.
    fn default() -> Self {
        Self {
            paired: false,
            per_cell: false,
            per_gene: false,
            per_contig: false,
            gene_tag: None,
            assigned_tag: None,
            skip_regex: None,
            spliced: false,
            ignore_tlen: false,
            ignore_umi: false,
            read_length: false,
            soft_clip_threshold: 4,
            subset: None,
            mapping_quality: None,
            detection_method: None,
            whole_contig: false,
            unmapped_reads: UnmappedReads::default(),
            unpaired_reads: UnpairedReads::default(),
            chimeric_pairs: ChimericPairs::default(),
            random_seed: None,
        }
    }
}

/// Everything [`crate::network::ReadDeduplicator`] reads off Python's
/// `options` object.
#[derive(Debug, Clone)]
pub struct DedupOptions {
    /// Which clustering method to apply (`--method`).
    pub method: ClusterMethod,
    /// The edit-distance threshold (`--edit-distance-threshold`).
    pub threshold: u32,
    /// Filter grouped UMIs against a whitelist (`--filter-umi`).
    pub filter_umi: bool,
    /// The UMI whitelist file (`--umi-whitelist`).
    pub umi_whitelist: Option<std::path::PathBuf>,
    /// The read-2 UMI whitelist file (`--umi-whitelist-paired`).
    pub umi_whitelist_paired: Option<std::path::PathBuf>,
}

impl Default for DedupOptions {
    /// Python's defaults: `method="directional"`, `threshold=1`.
    fn default() -> Self {
        Self {
            method: ClusterMethod::Directional,
            threshold: 1,
            filter_umi: false,
            umi_whitelist: None,
            umi_whitelist_paired: None,
        }
    }
}

/// Stream, format and temporary-file settings shared by every subcommand
/// (Python source: the `add_pipe_options` and
/// `add_s_b_cram_format_options` blocks of `Utilities.Start`).
#[derive(Debug, Clone, Default)]
pub struct IoOptions {
    /// Input file, or standard input when absent (`-I/--stdin`).
    pub stdin: Option<std::path::PathBuf>,
    /// Output file, or standard output when absent (`-S/--stdout`).
    pub stdout: Option<std::path::PathBuf>,
    /// Log file, or the log stream when absent (`-L/--log`).
    pub stdlog: Option<std::path::PathBuf>,
    /// Send the log to standard error (`--log2stderr`).
    pub log2stderr: bool,
    /// Verbosity (`-v/--verbose`).
    pub loglevel: i32,
    /// Directory for temporary files (`--temp-dir`).
    pub tmpdir: Option<std::path::PathBuf>,
    /// Input is SAM (`--in-sam`).
    pub in_sam: bool,
    /// Output is SAM (`--out-sam`).
    pub out_sam: bool,
    /// Explicit input format (`--in-format`).
    pub in_format: Option<AlignmentFormat>,
    /// Explicit output format (`--out-format`).
    pub out_format: Option<AlignmentFormat>,
    /// Skip the final coordinate sort (`--no-sort-output`).
    pub no_sort_output: bool,
    /// Extra format options for the reader (`--input-format-options`).
    pub input_options: Option<String>,
    /// Extra format options for the writer (`--output-format-options`).
    pub output_options: Option<String>,
    /// CRAM reference FASTA (`--reference-fasta`).
    pub reference_filename: Option<std::path::PathBuf>,
    /// Restrict processing to one contig (`--chrom`).
    pub chrom: Option<String>,
    /// gzip compression level.
    pub compresslevel: u32,
}

/// The barcode-extraction options `extract` and `whitelist` share
/// (Python source: the `add_extract_options` block of `Utilities.Start`).
#[derive(Debug, Clone, Default)]
pub struct ExtractOptions {
    /// `--extract-method`.
    pub extract_method: String,
    /// `--bc-pattern`.
    pub pattern: Option<String>,
    /// `--bc-pattern2`.
    pub pattern2: Option<String>,
    /// `--3prime`.
    pub prime3: bool,
    /// `--read2-in`.
    pub read2_in: Option<std::path::PathBuf>,
    /// `--read2-out`.
    pub read2_out: Option<std::path::PathBuf>,
    /// Extract from read 2 only.
    pub read2_only: bool,
    /// `--filtered-out`.
    pub filtered_out: Option<std::path::PathBuf>,
    /// `--filtered-out2`.
    pub filtered_out2: Option<std::path::PathBuf>,
}

/// Whether the configured patterns extract a cell barcode, a UMI, or both
/// (Python source: the `(extract_cell, extract_umi)` tuple
/// `validateExtractOptions` returns).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExtractRoles {
    /// The patterns mark cell-barcode bases.
    pub extract_cell: bool,
    /// The patterns mark UMI bases.
    pub extract_umi: bool,
}

/// Validates the barcode-extraction options and reports what they extract
/// (Python source: `Utilities.validateExtractOptions`).
///
/// Runs the same checks in the same order, so the first failure a user sees is
/// the one Python would report. Note that Python *mutates* the options object
/// here (compiling the regexes and copying `pattern` into `pattern2`); this
/// takes `&mut` for the same reason, and the mutation is documented on the
/// parameter rather than left implicit.
///
/// # Errors
///
/// Returns [`PortError::Value`] with Python's exact message for each of:
/// `--read2-only` without `--bc-pattern2` or with `--bc-pattern`; no pattern
/// at all; `--bc-pattern2` without `--read2-in`; `--filtered-out2` without
/// `--read2-in`; only one of `--filtered-out`/`--filtered-out2` for paired
/// input; an invalid regex; and patterns with no UMI bases or groups.
pub fn validate_extract_options(options: &mut ExtractOptions) -> Result<ExtractRoles> {
    if options.read2_only {
        if options.pattern2.is_none() {
            return Err(PortError::value(
                "Must supply --bc-pattern2 if extracting from just read2",
            ));
        }
        if options.pattern.is_some() {
            return Err(PortError::value(
                "Don't supply --bc-pattern if extracting from just read2",
            ));
        }
    }

    if options.pattern.is_none() && options.pattern2.is_none() {
        return Err(PortError::value(if options.read2_in.is_none() {
            "Must supply --bc-pattern for single-end"
        } else {
            "Must supply --bc-pattern and/or --bc-pattern2 if paired-end "
        }));
    }

    if options.pattern2.is_some() && options.read2_in.is_none() {
        return Err(PortError::value(
            "must specify a paired fastq ``--read2-in``",
        ));
    }

    if options.filtered_out2.is_some() && options.read2_in.is_none() {
        return Err(PortError::value(
            "Cannot use --filtered-out2 without read2 input (--read2-in)",
        ));
    }

    if (options.read2_in.is_some()
        && options.filtered_out.is_some()
        && options.filtered_out2.is_none())
        || (options.filtered_out2.is_some() && options.filtered_out.is_none())
    {
        return Err(PortError::value(
            "Must supply both --filtered-out and --filtered-out2\
             to write out filtered reads for paired end",
        ));
    }

    let mut extract_cell = false;
    let mut extract_umi = false;

    match options.extract_method.as_str() {
        "regex" => {
            for (pattern, label) in [
                (options.pattern.as_deref(), "--bc-pattern"),
                (options.pattern2.as_deref(), "--bc-pattern2"),
            ] {
                let Some(pattern) = pattern else { continue };
                let compiled = crate::pattern::Pattern::parse(pattern).map_err(|_| {
                    PortError::value(format!("{label} '{pattern}' is not a valid regex"))
                })?;
                for group in compiled.group_names() {
                    if group.starts_with("cell_") {
                        extract_cell = true;
                    } else if group.starts_with("umi_") {
                        extract_umi = true;
                    }
                }
            }
        }
        "string" => {
            for pattern in [options.pattern.as_deref(), options.pattern2.as_deref()]
                .into_iter()
                .flatten()
            {
                if pattern.contains('C') {
                    extract_cell = true;
                }
                if pattern.contains('N') {
                    extract_umi = true;
                }
            }
        }
        other => return Err(PortError::value(format!("unknown extract method: {other}"))),
    }

    if !extract_umi {
        let (p1, p2) = (
            options.pattern.as_deref().unwrap_or("None"),
            options.pattern2.as_deref().unwrap_or("None"),
        );
        return Err(PortError::value(if options.extract_method == "string" {
            format!("barcode pattern(s) do not include any umi bases (marked with 'Ns') {p1}, {p2}")
        } else {
            format!(
                "barcode regex(es) do not include any umi groups (starting with 'umi_') {p1}, {p2}"
            )
        }));
    }

    Ok(ExtractRoles {
        extract_cell,
        extract_umi,
    })
}

/// The SAM-related options `validate_sam_options` checks.
///
/// A view over the fields of [`BundleOptions`] and [`IoOptions`] the
/// validation actually reads, plus the two flags that live only on the command
/// line.
#[derive(Debug, Clone, Default)]
pub struct SamValidationInput {
    /// From [`BundleOptions`].
    pub per_gene: bool,
    /// From [`BundleOptions`].
    pub per_contig: bool,
    /// From [`BundleOptions`].
    pub gene_tag: Option<String>,
    /// From [`BundleOptions`].
    pub assigned_tag: Option<String>,
    /// `--gene-transcript-map`.
    pub gene_transcript_map: Option<std::path::PathBuf>,
    /// `--extract-umi-method`.
    pub get_umi_method: UmiMethod,
    /// `--umi-tag`.
    pub umi_tag: Option<String>,
    /// `--cell-tag`.
    pub cell_tag: Option<String>,
    /// `--per-cell`.
    pub per_cell: bool,
    /// `--skip-tags-regex`.
    pub skip_regex: Option<String>,
    /// `--paired`.
    pub paired: bool,
    /// `--unmapped-reads`.
    pub unmapped_reads: UnmappedReads,
    /// `--unpaired-reads`.
    pub unpaired_reads: UnpairedReads,
    /// `--chimeric-pairs`.
    pub chimeric_pairs: ChimericPairs,
    /// The legacy `--output-unmapped` flag.
    pub output_unmapped: bool,
    /// Which of the paired-only flags the user actually typed.
    ///
    /// Python inspects `" ".join(sys.argv)` for these, so a flag left at its
    /// default is treated differently from one explicitly set to that default.
    /// Passing the set of flags that appeared reproduces that without
    /// re-parsing the command line.
    pub flags_present: Vec<String>,
}

/// Validates the SAM-related option combinations
/// (Python source: `Utilities.validateSamOptions`).
///
/// `group` relaxes three checks: only the `group` command may write unmapped,
/// chimeric or unpaired reads through without grouping them.
///
/// Returns any warnings raised, so the caller can log them; Python calls
/// `warn()` directly from inside the validator, which a library should not do.
///
/// # Errors
///
/// Returns [`PortError::Value`] with Python's exact message for each invalid
/// combination — see the Python source for the full list, which this mirrors
/// check for check and in order.
pub fn validate_sam_options(options: &mut SamValidationInput, group: bool) -> Result<Vec<String>> {
    let mut warnings = Vec::new();

    if options.per_gene {
        if options.gene_tag.is_some() && options.per_contig {
            return Err(PortError::value(
                "need to use either --per-contig OR --gene-tag, please do not provide both",
            ));
        }
        if !options.per_contig && options.gene_tag.is_none() {
            return Err(PortError::value(
                "for per-gene applications, must supply --per-contig or --gene-tag",
            ));
        }
    }
    if options.per_contig && !options.per_gene {
        return Err(PortError::value("need to use --per-gene with --per-contig"));
    }
    if options.gene_tag.is_some() && !options.per_gene {
        return Err(PortError::value("need to use --per-gene with --gene_tag"));
    }
    if options.gene_transcript_map.is_some() && !options.per_contig {
        return Err(PortError::value(
            "need to use --per-contig and --per-genewith --gene-transcript-map",
        ));
    }

    if options.get_umi_method == UmiMethod::Tag {
        if options.umi_tag.is_none() {
            return Err(PortError::value("Need to supply the --umi-tag option"));
        }
        if options.per_cell && options.cell_tag.is_none() {
            return Err(PortError::value("Need to supply the --cell-tag option"));
        }
    }

    // Python assigns `options.assigned_tag = options.gene_tag` when unset.
    if options.assigned_tag.is_none() {
        options.assigned_tag = options.gene_tag.clone();
    }

    if let Some(pattern) = &options.skip_regex {
        Regex::new(pattern).map_err(|_| {
            PortError::value(format!("skip-regex '{pattern}' is not a valid regex"))
        })?;
    }

    if !group {
        if options.unmapped_reads == UnmappedReads::Output {
            return Err(PortError::value(
                "Cannot use --unmapped-reads=output. If you want to retain unmapped \
                 without deduplicating them, use the group command",
            ));
        }
        if options.chimeric_pairs == ChimericPairs::Output {
            return Err(PortError::value(
                "Cannot use --chimeric-pairs=output. If you want to retain chimeric \
                 read pairs without deduplicating them, use the group command",
            ));
        }
        if options.unpaired_reads == UnpairedReads::Output {
            return Err(PortError::value(
                "Cannot use --unpaired-reads=output. If you want to retain unmapped \
                 without deduplicating them, use the group command",
            ));
        }
    }

    if options.paired {
        if options.chimeric_pairs == ChimericPairs::Use {
            warnings.push(
                "Chimeric read pairs are being used. Some read pair UMIs may be \
                 grouped/deduplicated using just the mapping coordinates from read1."
                    .to_string(),
            );
        }
        if options.unpaired_reads == UnpairedReads::Use {
            warnings.push(
                "Unpaired read pairs are being used. Some read pair UMIs may be \
                 grouped/deduplicated using just the mapping coordinates from read1."
                    .to_string(),
            );
        }
        if options.unmapped_reads == UnmappedReads::Use {
            warnings.push(
                "Unmapped read pairs are being used. Some read pair UMIs may be \
                 grouped/deduplicated using just the mapping coordinates from read1."
                    .to_string(),
            );
        }
    }

    let typed = |flag: &str| options.flags_present.iter().any(|f| f == flag);

    if (typed("--umi-tag") || typed("--cell-tag")) && options.get_umi_method != UmiMethod::Tag {
        return Err(PortError::value(
            "--umi-tag and/or --cell-tag options provided. \
             Need to set --extract-umi-method=tag",
        ));
    }
    if options.unmapped_reads == UnmappedReads::Use && !options.paired {
        return Err(PortError::value(
            "--unmapped-reads=use is only compatible with paired end reads (--paired)",
        ));
    }
    for (flag, message) in [
        (
            "--chimeric-pairs",
            "--chimeric-pairs is only compatible with paired end reads (--paired)",
        ),
        (
            "--unpaired-reads",
            "--unpaired-reads is only compatible with paired end reads (--paired)",
        ),
        (
            "--ignore-tlen",
            "--ignore-tlen is only compatible with paired end reads (--paired)",
        ),
    ] {
        if typed(flag) && !options.paired {
            return Err(PortError::value(message));
        }
    }

    if options.output_unmapped {
        warnings.push(
            "--output-unmapped will be removed in the near future. \
             Use --unmapped-reads=output instead"
                .to_string(),
        );
        if typed("--unmapped_reads") || typed("--unmapped-reads") {
            return Err(PortError::value(
                "Do not use --output-unmapped in combination with\
                 --unmapped-reads. Just use --unmapped-reads",
            ));
        }
        options.unmapped_reads = UnmappedReads::Output;
    }

    Ok(warnings)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn extract_options() -> ExtractOptions {
        ExtractOptions {
            extract_method: "string".to_string(),
            ..ExtractOptions::default()
        }
    }

    #[test]
    fn three_way_options_parse_and_reject() {
        assert_eq!(
            "discard".parse::<UnmappedReads>().unwrap(),
            UnmappedReads::Discard
        );
        assert_eq!(
            "output".parse::<UnpairedReads>().unwrap(),
            UnpairedReads::Output
        );
        assert_eq!("use".parse::<ChimericPairs>().unwrap(), ChimericPairs::Use);
        assert!("maybe".parse::<UnmappedReads>().is_err());
    }

    #[test]
    fn detection_and_umi_method_parsing() {
        assert_eq!(
            "NH".parse::<DetectionMethod>().unwrap(),
            DetectionMethod::Nh
        );
        assert!("nh".parse::<DetectionMethod>().is_err());
        assert_eq!("tag".parse::<UmiMethod>().unwrap(), UmiMethod::Tag);
        let err = "other".parse::<UmiMethod>().unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m == "Unknown UMI extraction method"));
    }

    /// Verified against a live `umi_tools dedup` run's own parameter log,
    /// which prints every resolved option: `chimeric_pairs: use`,
    /// `unpaired_reads: use`, `unmapped_reads: discard`.
    #[test]
    fn bundle_defaults_match_python() {
        let d = BundleOptions::default();
        assert_eq!(d.soft_clip_threshold, 4);
        assert!(!d.paired);
        assert_eq!(d.unmapped_reads, UnmappedReads::Discard);
        assert_eq!(d.chimeric_pairs, ChimericPairs::Use);
        assert_eq!(d.unpaired_reads, UnpairedReads::Use);
        assert!(d.detection_method.is_none());
    }

    #[test]
    fn dedup_defaults_match_python() {
        let d = DedupOptions::default();
        assert_eq!(d.method, ClusterMethod::Directional);
        assert_eq!(d.threshold, 1);
    }

    #[test]
    fn extract_validation_needs_a_pattern() {
        let mut o = extract_options();
        let err = validate_extract_options(&mut o).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m == "Must supply --bc-pattern for single-end"));

        let mut o = ExtractOptions {
            read2_in: Some("r2.fastq".into()),
            ..extract_options()
        };
        let err = validate_extract_options(&mut o).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.starts_with("Must supply --bc-pattern and/or --bc-pattern2")));
    }

    #[test]
    fn extract_validation_reports_roles_for_string_patterns() {
        let mut o = ExtractOptions {
            pattern: Some("CCCNNN".to_string()),
            ..extract_options()
        };
        let roles = validate_extract_options(&mut o).unwrap();
        assert!(roles.extract_cell);
        assert!(roles.extract_umi);

        let mut o = ExtractOptions {
            pattern: Some("NNN".to_string()),
            ..extract_options()
        };
        assert!(!validate_extract_options(&mut o).unwrap().extract_cell);
    }

    #[test]
    fn extract_validation_reports_roles_for_regex_patterns() {
        let mut o = ExtractOptions {
            extract_method: "regex".to_string(),
            pattern: Some("^(?P<cell_1>.{6})(?P<umi_1>.{10})".to_string()),
            ..ExtractOptions::default()
        };
        let roles = validate_extract_options(&mut o).unwrap();
        assert!(roles.extract_cell && roles.extract_umi);
    }

    #[test]
    fn extract_validation_rejects_patterns_without_a_umi() {
        let mut o = ExtractOptions {
            pattern: Some("CCC".to_string()),
            ..extract_options()
        };
        let err = validate_extract_options(&mut o).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.starts_with("barcode pattern(s) do not include any umi bases")));

        let mut o = ExtractOptions {
            extract_method: "regex".to_string(),
            pattern: Some("^(?P<cell_1>.{6})".to_string()),
            ..ExtractOptions::default()
        };
        let err = validate_extract_options(&mut o).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.starts_with("barcode regex(es) do not include any umi groups")));
    }

    #[test]
    fn extract_validation_rejects_an_invalid_regex() {
        let mut o = ExtractOptions {
            extract_method: "regex".to_string(),
            pattern: Some("(?P<umi_1>.{3}".to_string()),
            ..ExtractOptions::default()
        };
        let err = validate_extract_options(&mut o).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.contains("is not a valid regex")));
    }

    #[test]
    fn extract_validation_checks_read2_only() {
        let mut o = ExtractOptions {
            read2_only: true,
            ..extract_options()
        };
        assert!(validate_extract_options(&mut o)
            .unwrap_err()
            .to_string()
            .contains("Must supply --bc-pattern2"));

        let mut o = ExtractOptions {
            read2_only: true,
            pattern: Some("NNN".to_string()),
            pattern2: Some("NNN".to_string()),
            ..extract_options()
        };
        assert!(validate_extract_options(&mut o)
            .unwrap_err()
            .to_string()
            .contains("Don't supply --bc-pattern"));
    }

    #[test]
    fn extract_validation_checks_filtered_out_pairs() {
        let mut o = ExtractOptions {
            pattern: Some("NNN".to_string()),
            filtered_out2: Some("f2".into()),
            ..extract_options()
        };
        assert!(validate_extract_options(&mut o)
            .unwrap_err()
            .to_string()
            .contains("Cannot use --filtered-out2 without read2 input"));

        let mut o = ExtractOptions {
            pattern: Some("NNN".to_string()),
            read2_in: Some("r2".into()),
            filtered_out: Some("f1".into()),
            ..extract_options()
        };
        assert!(validate_extract_options(&mut o)
            .unwrap_err()
            .to_string()
            .contains("Must supply both --filtered-out and --filtered-out2"));
    }

    fn sam_input() -> SamValidationInput {
        SamValidationInput::default()
    }

    #[test]
    fn sam_validation_per_gene_requires_a_source() {
        let mut o = SamValidationInput {
            per_gene: true,
            ..sam_input()
        };
        assert!(validate_sam_options(&mut o, false)
            .unwrap_err()
            .to_string()
            .contains("must supply --per-contig or --gene-tag"));

        let mut o = SamValidationInput {
            per_gene: true,
            per_contig: true,
            gene_tag: Some("XT".to_string()),
            ..sam_input()
        };
        assert!(validate_sam_options(&mut o, false)
            .unwrap_err()
            .to_string()
            .contains("do not provide both"));
    }

    #[test]
    fn sam_validation_per_contig_requires_per_gene() {
        let mut o = SamValidationInput {
            per_contig: true,
            ..sam_input()
        };
        assert!(validate_sam_options(&mut o, false)
            .unwrap_err()
            .to_string()
            .contains("need to use --per-gene with --per-contig"));
    }

    #[test]
    fn sam_validation_defaults_assigned_tag_to_gene_tag() {
        let mut o = SamValidationInput {
            per_gene: true,
            gene_tag: Some("XT".to_string()),
            ..sam_input()
        };
        validate_sam_options(&mut o, false).unwrap();
        assert_eq!(o.assigned_tag.as_deref(), Some("XT"));
    }

    #[test]
    fn sam_validation_tag_method_requires_tags() {
        let mut o = SamValidationInput {
            get_umi_method: UmiMethod::Tag,
            ..sam_input()
        };
        assert!(validate_sam_options(&mut o, false)
            .unwrap_err()
            .to_string()
            .contains("Need to supply the --umi-tag option"));

        let mut o = SamValidationInput {
            get_umi_method: UmiMethod::Tag,
            umi_tag: Some("RX".to_string()),
            per_cell: true,
            ..sam_input()
        };
        assert!(validate_sam_options(&mut o, false)
            .unwrap_err()
            .to_string()
            .contains("Need to supply the --cell-tag option"));
    }

    #[test]
    fn sam_validation_rejects_output_modes_outside_group() {
        for (label, mut o) in [
            (
                "unmapped",
                SamValidationInput {
                    unmapped_reads: UnmappedReads::Output,
                    ..sam_input()
                },
            ),
            (
                "chimeric",
                SamValidationInput {
                    chimeric_pairs: ChimericPairs::Output,
                    ..sam_input()
                },
            ),
            (
                "unpaired",
                SamValidationInput {
                    unpaired_reads: UnpairedReads::Output,
                    ..sam_input()
                },
            ),
        ] {
            assert!(
                validate_sam_options(&mut o, false).is_err(),
                "{label} should be rejected outside group"
            );
            assert!(
                validate_sam_options(&mut o, true).is_ok(),
                "{label} should be accepted by group"
            );
        }
    }

    #[test]
    fn sam_validation_paired_only_flags() {
        let mut o = SamValidationInput {
            flags_present: vec!["--chimeric-pairs".to_string()],
            ..sam_input()
        };
        assert!(validate_sam_options(&mut o, false)
            .unwrap_err()
            .to_string()
            .contains("only compatible with paired end"));

        let mut o = SamValidationInput {
            paired: true,
            flags_present: vec!["--chimeric-pairs".to_string()],
            ..sam_input()
        };
        assert!(validate_sam_options(&mut o, false).is_ok());
    }

    #[test]
    fn sam_validation_rejects_tags_without_the_tag_method() {
        let mut o = SamValidationInput {
            flags_present: vec!["--umi-tag".to_string()],
            ..sam_input()
        };
        assert!(validate_sam_options(&mut o, false)
            .unwrap_err()
            .to_string()
            .contains("Need to set --extract-umi-method=tag"));
    }

    #[test]
    fn sam_validation_warns_on_use_modes() {
        let mut o = SamValidationInput {
            paired: true,
            chimeric_pairs: ChimericPairs::Use,
            unpaired_reads: UnpairedReads::Use,
            unmapped_reads: UnmappedReads::Use,
            ..sam_input()
        };
        let warnings = validate_sam_options(&mut o, false).unwrap();
        assert_eq!(warnings.len(), 3);
    }

    #[test]
    fn sam_validation_migrates_the_legacy_output_unmapped_flag() {
        let mut o = SamValidationInput {
            output_unmapped: true,
            ..sam_input()
        };
        let warnings = validate_sam_options(&mut o, true).unwrap();
        assert!(warnings.iter().any(|w| w.contains("will be removed")));
        assert_eq!(o.unmapped_reads, UnmappedReads::Output);

        let mut o = SamValidationInput {
            output_unmapped: true,
            flags_present: vec!["--unmapped-reads".to_string()],
            ..sam_input()
        };
        assert!(validate_sam_options(&mut o, true)
            .unwrap_err()
            .to_string()
            .contains("Do not use --output-unmapped in combination"));
    }

    #[test]
    fn sam_validation_rejects_an_invalid_skip_regex() {
        let mut o = SamValidationInput {
            skip_regex: Some("(unclosed".to_string()),
            ..sam_input()
        };
        assert!(validate_sam_options(&mut o, false)
            .unwrap_err()
            .to_string()
            .contains("is not a valid regex"));
    }
}
