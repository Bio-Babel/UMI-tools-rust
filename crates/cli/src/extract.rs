//! `umi_tools extract` — move the UMI into the read name
//! (Python source: `umi_tools/extract.py`).
//!
//! Reads FASTQ, pulls the UMI (and optionally the cell barcode) out of each
//! read's sequence, appends them to the read identifier, and writes the
//! shortened read back out. Everything downstream — the aligner, then `dedup`
//! or `count` — carries the barcode in the read name from here on.

use std::io::{BufWriter, Write};
use std::path::Path;

use umi_tools_core::extract::{ExtractConfig, ExtractFilterAndUpdate, ExtractMethod};
use umi_tools_core::options::validate_extract_options;
use umi_tools_core::umi_methods::{FastqIterator, JoinedFastqIterator, Record};
use umi_tools_core::utilities::{self, InputFile, OutputFile};
use umi_tools_core::whitelist::{get_user_defined_barcodes, ErrorCorrectionMode};
use umi_tools_core::{PortError, Result};

use crate::documentation;
use crate::options::{self, OptionGroups};

/// How often the progress line is printed
/// (Python source: `displayMax = 100000`).
const DISPLAY_EVERY: u64 = 100_000;

/// The subcommand's short description, shown in `--help`.
pub const ABOUT: &str = "extract - Extract UMI barcode from a read and add it to the read name";

/// Runs `umi_tools extract` (Python source: `extract.main`).
///
/// `args` is the argument list without the subcommand name; `argv` includes it,
/// for the log header and for the typed-flag inspection the SAM validator
/// needs.
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Value`] for an invalid option combination, and
/// [`umi_tools_core::PortError::Io`] or [`umi_tools_core::PortError::Format`] for input problems.
pub fn main(args: &[String], argv: &[String]) -> Result<i32> {
    if crate::wants_help(args) {
        // `optparse`'s layout, reproduced verbatim; see `crate::help`.
        print!("{}", crate::help::EXTRACT);
        return Ok(0);
    }

    let started = std::time::Instant::now();
    let groups = OptionGroups {
        extract: true,
        group_dedup: false,
        umi_grouping: false,
        sam: false,
        ..OptionGroups::default()
    };

    let command = options::build_parser("extract", ABOUT, groups)
        .after_help(documentation::GENERIC_DOCSTRING_WE)
        .arg(
            clap::Arg::new("whitelist")
                .overrides_with("whitelist")
                .long("whitelist")
                .num_args(1)
                .help("a whitelist of accepted cell barcodes"),
        )
        .arg(
            clap::Arg::new("blacklist")
                .overrides_with("blacklist")
                .long("blacklist")
                .num_args(1)
                .help("a blacklist of rejected cell barcodes"),
        )
        .arg(
            clap::Arg::new("error_correct_cell")
                .overrides_with("error_correct_cell")
                .long("error-correct-cell")
                .action(clap::ArgAction::SetTrue)
                .help("correct errors in the cell barcode to the whitelist"),
        )
        .arg(
            clap::Arg::new("read2_out")
                .overrides_with("read2_out")
                .long("read2-out")
                .num_args(1)
                .help("file to output processed paired read to"),
        )
        .arg(
            clap::Arg::new("read2_stdout")
                .overrides_with("read2_stdout")
                .long("read2-stdout")
                .action(clap::ArgAction::SetTrue)
                .help("send read2 to stdout, discarding read1"),
        )
        .arg(
            clap::Arg::new("read2_only")
                .overrides_with("read2_only")
                .long("read2-only")
                .action(clap::ArgAction::SetTrue)
                .help("extract from read 2 only"),
        )
        .arg(
            clap::Arg::new("retain_umi")
                .overrides_with("retain_umi")
                .long("retain-umi")
                .action(clap::ArgAction::SetTrue)
                .help("retain the UMI sequence on the read"),
        )
        .arg(
            clap::Arg::new("reads_subset")
                .overrides_with("reads_subset")
                .long("subset-reads")
                .alias("reads-subset")
                .num_args(1)
                .help("only extract from the first N reads"),
        )
        .arg(
            clap::Arg::new("reconcile")
                .overrides_with("reconcile")
                .long("reconcile-pairs")
                .action(clap::ArgAction::SetTrue)
                .help("allow read2 input containing reads absent from read1"),
        )
        .arg(
            clap::Arg::new("umi_separator")
                .overrides_with("umi_separator")
                .long("umi-separator")
                .num_args(1)
                .help("separator to use to add the UMI to the read name"),
        )
        .arg(
            clap::Arg::new("either_read")
                .overrides_with("either_read")
                .long("either-read")
                .action(clap::ArgAction::SetTrue)
                .help("the UMI may be on either read"),
        )
        .arg(
            clap::Arg::new("either_read_resolve")
                .overrides_with("either_read_resolve")
                .long("either-read-resolve")
                .num_args(1)
                .help("how to resolve a UMI found on both reads: discard or quality"),
        )
        .arg(
            clap::Arg::new("umi_correct_log")
                .overrides_with("umi_correct_log")
                .long("umi-correct-log")
                .num_args(1)
                .help("file logging UMI error correction"),
        );

    let matches = match command
        .try_get_matches_from(std::iter::once("extract".to_string()).chain(args.iter().cloned()))
    {
        Ok(m) => m,
        Err(e) => {
            e.print().ok();
            return Ok(i32::from(e.use_stderr()));
        }
    };

    let io = options::io_options(&matches)?;
    let run_id = options::start(&io, argv)?;

    let mut extract_opts = options::extract_options(&matches)?;
    let roles = validate_extract_options(&mut extract_opts)?;

    let whitelist_path = matches.get_one::<String>("whitelist").map(Path::new);
    let blacklist_path = matches.get_one::<String>("blacklist").map(Path::new);

    if whitelist_path.is_some() && !roles.extract_cell {
        return Err(PortError::value(
            if extract_opts.extract_method == "string" {
                "barcode pattern(s) do not include any cell bases (marked with 'Cs')"
            } else {
                "barcode regex(es) do not include any cell groups (starting with 'cell_')"
            },
        ));
    }

    let retain_umi = matches.get_flag("retain_umi");
    if retain_umi && extract_opts.extract_method != "regex" {
        return Err(PortError::value(
            "option --retain-umi only works with --extract-method=regex",
        ));
    }

    let either_read = matches.get_flag("either_read");
    if either_read {
        if roles.extract_cell {
            return Err(PortError::value(
                "Option to extract from either read (--either-read) \
                 is not currently compatible with cell barcode extraction",
            ));
        }
        if extract_opts.extract_method != "regex" {
            return Err(PortError::value(
                "Option to extract from either read (--either-read)\
                 requires --extract-method=regex",
            ));
        }
        if extract_opts.pattern.is_none() || extract_opts.pattern2.is_none() {
            return Err(PortError::value(
                "Option to extract from either read (--either-read)\
                 requires --bc-pattern=[PATTERN1] and--bc-pattern2=[PATTERN2]",
            ));
        }
    }

    let config = ExtractConfig {
        method: extract_opts.extract_method.parse::<ExtractMethod>()?,
        pattern: extract_opts.pattern.clone(),
        pattern2: extract_opts.pattern2.clone(),
        prime3: extract_opts.prime3,
        extract_cell: roles.extract_cell,
        quality_encoding: options::quality_encoding(&matches)?,
        quality_filter_threshold: matches
            .get_one::<String>("quality_filter_threshold")
            .map(|v| v.parse::<i32>())
            .transpose()
            .map_err(|e| PortError::value(format!("--quality-filter-threshold: {e}")))?,
        quality_filter_mask: matches
            .get_one::<String>("quality_filter_mask")
            .map(|v| v.parse::<i32>())
            .transpose()
            .map_err(|e| PortError::value(format!("--quality-filter-mask: {e}")))?,
        filter_umi_barcode: false,
        filter_cell_barcode: whitelist_path.is_some(),
        retain_umi,
        either_read,
        either_read_resolve: options::either_read_resolve(&matches)?.unwrap_or_default(),
        umi_separator: matches
            .get_one::<String>("umi_separator")
            .cloned()
            .unwrap_or_else(|| "_".to_string()),
    };

    if (config.quality_filter_threshold.is_some() || config.quality_filter_mask.is_some())
        && config.quality_encoding.is_none()
    {
        return Err(PortError::value(
            "must provide a quality encoding (--quality-encoding) to filter UMIs by quality \
             (--quality-filter-threshold) or mask low quality bases with (--quality-filter-mask)",
        ));
    }

    let mut extractor = ExtractFilterAndUpdate::new(config)?;

    if let Some(path) = whitelist_path {
        let mode = if matches.get_flag("error_correct_cell") {
            ErrorCorrectionMode::FromFile
        } else {
            ErrorCorrectionMode::None
        };
        let (whitelist, map) = get_user_defined_barcodes(path, None, mode, 1)?;
        extractor.cell_whitelist = Some(
            whitelist
                .iter()
                .map(|b| String::from_utf8_lossy(b).into_owned())
                .collect(),
        );
        extractor.false_to_true_map = map.map(|m| {
            m.into_iter()
                .filter_map(|(k, v)| {
                    v.map(|v| {
                        (
                            String::from_utf8_lossy(&k).into_owned(),
                            String::from_utf8_lossy(&v).into_owned(),
                        )
                    })
                })
                .collect()
        });
    }

    if let Some(path) = blacklist_path {
        use std::io::BufRead as _;
        let reader = utilities::open_file(path)?;
        let mut blacklist = std::collections::BTreeSet::new();
        for line in reader.lines() {
            let line = line.map_err(|e| PortError::io(path, e))?;
            if let Some(first) = line.trim().split('\t').next() {
                if !first.is_empty() {
                    blacklist.insert(first.to_string());
                }
            }
        }
        extractor.cell_blacklist = Some(blacklist);
    }

    let ignore_suffix = matches.get_flag("ignore_suffix");
    let reads_subset: Option<u64> = matches
        .get_one::<String>("reads_subset")
        .map(|v| v.parse::<u64>())
        .transpose()
        .map_err(|e| PortError::value(format!("--subset-reads: {e}")))?;

    let mut out = open_output(io.stdout.as_deref())?;
    let mut filtered_out = extract_opts
        .filtered_out
        .as_deref()
        .map(|p| utilities::create_file(p, false))
        .transpose()?;
    let mut filtered_out2 = extract_opts
        .filtered_out2
        .as_deref()
        .map(|p| utilities::create_file(p, false))
        .transpose()?;
    let mut read2_out = matches
        .get_one::<String>("read2_out")
        .map(|p| utilities::create_file(Path::new(p), false))
        .transpose()?;
    let read2_stdout = matches.get_flag("read2_stdout");

    let stdin_path = io.stdin.as_deref().ok_or_else(|| {
        PortError::value("reading FASTQ from standard input is not supported; use -I/--stdin")
    })?;
    let read1s = FastqIterator::new(utilities::open_file(stdin_path)?, ignore_suffix);

    utilities::info("Starting barcode extraction");
    let mut progress = 0u64;

    match extract_opts.read2_in.as_deref() {
        None => {
            for read in read1s {
                let mut read = read?;
                progress += 1;
                if progress.is_multiple_of(DISPLAY_EVERY) {
                    utilities::info(&format!("Parsed {progress} reads"));
                }

                let kept = extractor.process(&mut read, None)?;
                if let Some(limit) = reads_subset {
                    if extractor
                        .read_counts()
                        .get("Input Reads")
                        .copied()
                        .unwrap_or(0)
                        > limit
                    {
                        break;
                    }
                }
                if kept {
                    writeln!(out, "{read}").map_err(io_error)?;
                } else if let Some(f) = filtered_out.as_mut() {
                    writeln!(f, "{read}").map_err(io_error)?;
                }
            }
        }
        Some(read2_path) => {
            let read2s = FastqIterator::new(utilities::open_file(read2_path)?, ignore_suffix);
            let strict = !matches.get_flag("reconcile");
            for pair in JoinedFastqIterator::new(read1s, read2s, strict) {
                let (mut r1, mut r2) = pair?;
                progress += 1;
                if progress.is_multiple_of(DISPLAY_EVERY) {
                    utilities::info(&format!("Parsed {progress} reads"));
                }

                let kept = extractor.process(&mut r1, Some(&mut r2))?;
                if let Some(limit) = reads_subset {
                    if extractor
                        .read_counts()
                        .get("Input Reads")
                        .copied()
                        .unwrap_or(0)
                        > limit
                    {
                        break;
                    }
                }
                if !kept {
                    if let Some(f) = filtered_out.as_mut() {
                        writeln!(f, "{r1}").map_err(io_error)?;
                    }
                    if let Some(f) = filtered_out2.as_mut() {
                        writeln!(f, "{r2}").map_err(io_error)?;
                    }
                    continue;
                }

                if read2_stdout {
                    writeln!(out, "{r2}").map_err(io_error)?;
                } else {
                    writeln!(out, "{r1}").map_err(io_error)?;
                    if let Some(f) = read2_out.as_mut() {
                        writeln!(f, "{r2}").map_err(io_error)?;
                    }
                }
            }
        }
    }

    out.finish()?;
    for file in [read2_out.take(), filtered_out.take(), filtered_out2.take()]
        .into_iter()
        .flatten()
    {
        file.finish()?;
    }

    // Python prints the tallies in descending count order.
    let mut events: Vec<(&String, &u64)> = extractor.read_counts().iter().collect();
    events.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
    for (event, count) in events {
        utilities::info(&format!("{event}: {count}"));
    }

    if let Some(path) = matches.get_one::<String>("umi_correct_log") {
        let mut log = utilities::create_file(Path::new(path), false)?;
        writeln!(log, "umi\tcount_no_errors\tcount_errors").map_err(io_error)?;
        for (umi, counts) in &extractor.umi_whitelist_counts {
            writeln!(log, "{umi}\t{}\t{}", counts.no_error, counts.error).map_err(io_error)?;
        }
        log.finish()?;
    }

    options::stop(&io, &run_id, started.elapsed().as_secs())?;
    Ok(0)
}

/// A destination for extracted reads: a file, or standard output.
enum Output {
    File(Box<OutputFile>),
    Stdout(BufWriter<std::io::Stdout>),
}

impl Write for Output {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Output::File(f) => f.write(buf),
            Output::Stdout(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Output::File(f) => f.flush(),
            Output::Stdout(w) => w.flush(),
        }
    }
}

impl Output {
    fn finish(self) -> Result<()> {
        match self {
            Output::File(f) => f.finish(),
            Output::Stdout(mut w) => w.flush().map_err(io_error),
        }
    }
}

fn open_output(path: Option<&Path>) -> Result<Output> {
    match path {
        Some(p) => Ok(Output::File(Box::new(utilities::create_file(p, false)?))),
        None => Ok(Output::Stdout(BufWriter::new(std::io::stdout()))),
    }
}

fn io_error(e: std::io::Error) -> PortError {
    PortError::io(Path::new("<output>"), e)
}

/// Reads a FASTQ into memory, for callers that need random access.
///
/// # Errors
///
/// Propagates the iterator's errors.
pub fn read_fastq(reader: InputFile, remove_suffix: bool) -> Result<Vec<Record>> {
    FastqIterator::new(reader, remove_suffix).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn about_names_the_command() {
        assert!(ABOUT.starts_with("extract - "));
    }

    #[test]
    fn a_pattern_is_required() {
        let err = main(&[], &["extract".to_string()]).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.contains("Must supply --bc-pattern")));
    }

    #[test]
    fn retain_umi_requires_the_regex_method() {
        let args: Vec<String> = ["--bc-pattern=NNN", "--retain-umi"]
            .iter()
            .map(|s| (*s).to_string())
            .collect();
        let err = main(&args, &["extract".to_string()]).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.contains("--retain-umi only works with --extract-method=regex")));
    }

    #[test]
    fn quality_filtering_requires_an_encoding() {
        let args: Vec<String> = ["--bc-pattern=NNN", "--quality-filter-threshold=30"]
            .iter()
            .map(|s| (*s).to_string())
            .collect();
        let err = main(&args, &["extract".to_string()]).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.contains("must provide a quality encoding")));
    }

    #[test]
    fn either_read_requires_two_regex_patterns() {
        let args: Vec<String> = [
            "--extract-method=regex",
            "--bc-pattern=^(?P<umi_1>.{3})",
            "--either-read",
        ]
        .iter()
        .map(|s| (*s).to_string())
        .collect();
        let err = main(&args, &["extract".to_string()]).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.contains("--bc-pattern2")));
    }
}
