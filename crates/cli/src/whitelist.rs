//! `umi_tools whitelist` — pick the real cell barcodes
//! (Python source: `umi_tools/whitelist.py`).
//!
//! Scans a FASTQ, tallies how often each candidate cell barcode occurs, and
//! decides which of them correspond to real cells. Droplet protocols produce
//! far more distinct barcodes than there were cells — most are sequencing
//! errors or ambient RNA — so the tally has a characteristic "knee" separating
//! real cells from background, and this command finds it.

use std::collections::BTreeMap;
use std::io::{BufWriter, Write};
use std::path::Path;

use umi_tools_core::extract::{ExtractConfig, ExtractFilterAndUpdate, ExtractMethod};
use umi_tools_core::options::validate_extract_options;
use umi_tools_core::umi_methods::FastqIterator;
use umi_tools_core::utilities;
use umi_tools_core::whitelist::{
    error_detect_above_threshold, get_cell_whitelist, write_cell_thresholds, CellCounts,
    KneeMethod, ResolutionMethod,
};
use umi_tools_core::{PortError, Result};

use crate::documentation;
use crate::options::{self, OptionGroups};

/// How often the progress line is printed.
const DISPLAY_EVERY: u64 = 100_000;

/// Python's default cap on how many reads are scanned
/// (Python source: `subset_reads=100000000` in `parser.set_defaults`).
const DEFAULT_SUBSET_READS: u64 = 100_000_000;

/// The subcommand's short description.
pub const ABOUT: &str = "whitelist - Generates a whitelist of accepted cell barcodes";

/// Whether barcodes are ranked by read count or by distinct UMI count
/// (Python source: the `--method` choice, `reads` or `umis`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RankBy {
    Reads,
    Umis,
}

/// Runs `umi_tools whitelist` (Python source: `whitelist.main`).
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Value`] for an invalid option combination, and
/// [`umi_tools_core::PortError::Io`] or [`umi_tools_core::PortError::Format`] for input problems.
pub fn main(args: &[String], argv: &[String]) -> Result<i32> {
    if crate::wants_help(args) {
        // `optparse`'s layout, reproduced verbatim; see `crate::help`.
        print!("{}", crate::help::WHITELIST);
        return Ok(0);
    }

    let started = std::time::Instant::now();
    let groups = OptionGroups {
        extract: true,
        group_dedup: false,
        umi_grouping: false,
        ..OptionGroups::default()
    };

    let command = options::build_parser("whitelist", ABOUT, groups)
        .after_help(documentation::GENERIC_DOCSTRING_WE)
        .arg(value_arg(
            "plot_prefix",
            "plot-prefix",
            "prefix for the files describing the automated threshold detection",
        ))
        .arg(value_arg(
            "subset_reads",
            "subset-reads",
            "use the first N reads to identify the true cell barcodes",
        ))
        .arg(value_arg(
            "error_correct_threshold",
            "error-correct-threshold",
            "hamming distance for correction of barcodes to whitelist barcodes",
        ))
        .arg(value_arg(
            "method",
            "method",
            "use reads or unique umi counts per cell",
        ))
        .arg(value_arg(
            "knee_method",
            "knee-method",
            "use the distance or density method for detection of the knee",
        ))
        .arg(value_arg(
            "expect_cells",
            "expect-cells",
            "prior expectation on the upper limit on the number of cells sequenced",
        ))
        .arg(
            clap::Arg::new("allow_threshold_error")
                .overrides_with("allow_threshold_error")
                .long("allow-threshold-error")
                .action(clap::ArgAction::SetTrue)
                .help("don't fail when no threshold is selected"),
        )
        .arg(value_arg(
            "cell_number",
            "set-cell-number",
            "specify the number of cell barcodes to accept",
        ))
        .arg(value_arg(
            "ed_above_threshold",
            "ed-above-threshold",
            "detect CBs above the threshold which may be errors: discard or correct",
        ));

    let matches = match command
        .try_get_matches_from(std::iter::once("whitelist".to_string()).chain(args.iter().cloned()))
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
    if extract_opts.filtered_out.is_some() && extract_opts.extract_method != "regex" {
        return Err(PortError::value(
            "Reads will not be filtered unless extract method isset to regex \
             (--extract-method=regex)",
        ));
    }

    let knee_method: KneeMethod = matches
        .get_one::<String>("knee_method")
        .map_or(Ok(KneeMethod::Distance), |v| v.parse())?;
    let expect_cells: Option<u64> = parse_num(&matches, "expect_cells", "--expect-cells")?;
    let cell_number: Option<usize> = parse_num(&matches, "cell_number", "--set-cell-number")?;

    if expect_cells.is_some() {
        if knee_method == KneeMethod::Distance {
            return Err(PortError::value(
                "Cannot use --expect-cells with 'distance' knee method. Switch to \
                 --knee-method=density if you want to provide an expectation for the number \
                 of cells. Alternatively, if you know the number of cell barcodes, use \
                 --cell-number",
            ));
        }
        if cell_number.is_some() {
            return Err(PortError::value(
                "Cannot supply both --expect-cells and --cell-number options",
            ));
        }
    }

    let roles = validate_extract_options(&mut extract_opts)?;
    if !roles.extract_cell {
        return Err(PortError::value(
            if extract_opts.extract_method == "string" {
                "barcode pattern(s) do not include any cell bases (marked with 'Cs')"
            } else {
                "barcode regex(es) do not include any cell groups (starting with 'cell_')"
            },
        ));
    }

    let mut extractor = ExtractFilterAndUpdate::new(ExtractConfig {
        method: extract_opts.extract_method.parse::<ExtractMethod>()?,
        pattern: extract_opts.pattern.clone(),
        pattern2: extract_opts.pattern2.clone(),
        prime3: extract_opts.prime3,
        extract_cell: true,
        ..ExtractConfig::new()
    })?;

    let rank_by = match matches.get_one::<String>("method").map(String::as_str) {
        None | Some("reads") => RankBy::Reads,
        Some("umis") => RankBy::Umis,
        Some(other) => {
            return Err(PortError::value(format!(
                "--method must be reads or umis, not {other}"
            )))
        }
    };

    let stdin_path = io.stdin.as_deref().ok_or_else(|| {
        PortError::value("reading FASTQ from standard input is not supported; use -I/--stdin")
    })?;
    let subset_reads: u64 =
        parse_num(&matches, "subset_reads", "--subset-reads")?.unwrap_or(DEFAULT_SUBSET_READS);

    let mut counts = CellCounts::new();
    let mut cell_umis = CellUmiSets::default();
    let mut n_reads = 0u64;
    let mut n_cell_barcodes = 0u64;

    utilities::info("Starting barcode extraction");
    let mut filtered_out = extract_opts
        .filtered_out
        .as_deref()
        .map(|p| utilities::create_file(p, false))
        .transpose()?;

    let read1s = FastqIterator::new(utilities::open_file(stdin_path)?, false);
    match extract_opts.read2_in.as_deref() {
        None => {
            for read in read1s {
                let read = read?;
                if n_reads.is_multiple_of(DISPLAY_EVERY) {
                    utilities::info(&format!("Parsed {n_reads} reads"));
                }
                n_reads += 1;
                match extractor.get_barcodes(&read, None)? {
                    None => {
                        if let Some(f) = filtered_out.as_mut() {
                            writeln!(f, "{read}").map_err(io_error)?;
                        }
                    }
                    Some(values) => {
                        record(
                            &mut counts,
                            &mut cell_umis,
                            rank_by,
                            &values.cell,
                            &values.umi,
                        );
                        n_cell_barcodes += 1;
                    }
                }
                if n_cell_barcodes > subset_reads {
                    break;
                }
            }
        }
        Some(read2_path) => {
            let read2s = FastqIterator::new(utilities::open_file(read2_path)?, false);
            let mut filtered_out2 = extract_opts
                .filtered_out2
                .as_deref()
                .map(|p| utilities::create_file(p, false))
                .transpose()?;
            // Python zips the two iterators here rather than using
            // joinedFastqIterate, so it does not check the identifiers match.
            for (r1, r2) in read1s.zip(read2s) {
                let (r1, r2) = (r1?, r2?);
                if n_reads.is_multiple_of(DISPLAY_EVERY) {
                    utilities::info(&format!("Parsed {n_reads} reads"));
                }
                n_reads += 1;
                match extractor.get_barcodes(&r1, Some(&r2))? {
                    None => {
                        if let Some(f) = filtered_out.as_mut() {
                            writeln!(f, "{r1}").map_err(io_error)?;
                        }
                        if let Some(f) = filtered_out2.as_mut() {
                            writeln!(f, "{r2}").map_err(io_error)?;
                        }
                    }
                    Some(values) => {
                        record(
                            &mut counts,
                            &mut cell_umis,
                            rank_by,
                            &values.cell,
                            &values.umi,
                        );
                        n_cell_barcodes += 1;
                    }
                }
                if n_reads > subset_reads {
                    break;
                }
            }
            if let Some(f) = filtered_out2.take() {
                f.finish()?;
            }
        }
    }
    if let Some(f) = filtered_out.take() {
        f.finish()?;
    }

    utilities::info("Starting - whitelist determination");
    if rank_by == RankBy::Umis {
        for (cell, n) in cell_umis.iter() {
            counts.set(cell, n as u64);
        }
    }

    if let Some(n) = cell_number {
        if n > counts.len() {
            return Err(PortError::value(format!(
                "--set-cell-barcode option specifies more cell barcodes than the number of \
                 observed cell barcodes. This may be because --subset-reads was set to a value \
                 too low to capture reads from all cells. {} cell barcodes observed from {} \
                 parsed reads. Expected>= {n} cell barcodes",
                counts.len(),
                subset_reads
            )));
        }
    }

    let error_correct_threshold: u32 = parse_num(
        &matches,
        "error_correct_threshold",
        "--error-correct-threshold",
    )?
    .unwrap_or(1);

    let (knee, true_to_false) = get_cell_whitelist(
        &counts,
        knee_method,
        expect_cells,
        cell_number,
        error_correct_threshold,
    )?;

    if let Some(prefix) = matches.get_one::<String>("plot_prefix") {
        write_cell_thresholds(&knee, knee_method, Path::new(prefix), cell_number)?;
    }

    let mut whitelist = knee.barcodes.clone();
    let mut mapping = true_to_false;

    if let Some(mode) = matches.get_one::<String>("ed_above_threshold") {
        let method: ResolutionMethod = mode.parse()?;
        if whitelist.is_some() {
            // Most-abundant-first: the tie order decides which of two
            // mutually-erroneous barcodes is dropped.
            let (new_whitelist, new_map) = error_detect_above_threshold(
                &counts,
                &knee.ordered_barcodes,
                mapping.as_ref(),
                error_correct_threshold,
                method,
            )?;
            whitelist = Some(new_whitelist);
            mapping = Some(new_map);
        }
    }

    let mut out = BufWriter::new(std::io::stdout());
    let mut total_correct = 0u64;
    let mut total_corrected = 0u64;

    if let Some(whitelist) = whitelist.as_ref().filter(|w| !w.is_empty()) {
        utilities::info(&format!(
            "Top {} cell barcodes passed the selected threshold",
            whitelist.len()
        ));
        utilities::info("Writing out whitelist");
        for barcode in whitelist {
            total_correct += counts.get(barcode).unwrap_or(0);
            let (corrected, corrected_counts) = match mapping.as_ref() {
                Some(map) => {
                    let children: Vec<&Vec<u8>> = map
                        .get(barcode)
                        .map(|s| s.iter().collect())
                        .unwrap_or_default();
                    let names: Vec<String> = children
                        .iter()
                        .map(|c| String::from_utf8_lossy(c).into_owned())
                        .collect();
                    let child_counts: Vec<u64> = children
                        .iter()
                        .map(|c| counts.get(c).unwrap_or(0))
                        .collect();
                    total_corrected += child_counts.iter().sum::<u64>();
                    (
                        names.join(","),
                        child_counts
                            .iter()
                            .map(u64::to_string)
                            .collect::<Vec<_>>()
                            .join(","),
                    )
                }
                None => (String::new(), String::new()),
            };
            writeln!(
                out,
                "{}\t{corrected}\t{}\t{corrected_counts}",
                String::from_utf8_lossy(barcode),
                counts.get(barcode).unwrap_or(0)
            )
            .map_err(io_error)?;
        }
    } else {
        let msg = "No local minima was accepted. Recommend checking the plot output and \
                   counts per local minima (requires `--plot-prefix`option) and then \
                   re-running with manually selected threshold (`--set-cell-number` option)";
        if matches.get_flag("allow_threshold_error") {
            utilities::info(msg);
        } else {
            return Err(PortError::value(msg));
        }
    }
    out.flush().map_err(io_error)?;

    utilities::info(&format!("Parsed {n_reads} reads"));
    utilities::info(&format!(
        "{n_cell_barcodes} reads matched the barcode pattern"
    ));
    utilities::info(&format!("Found {} unique cell barcodes", counts.len()));
    if whitelist.as_ref().is_some_and(|w| !w.is_empty()) {
        utilities::info(&format!(
            "Found {total_correct} total reads matching the selected cell barcodes"
        ));
        utilities::info(&format!(
            "Found {total_corrected} total reads which can be error corrected to the \
             selected cell barcodes"
        ));
    }

    options::stop(&io, &run_id, started.elapsed().as_secs())?;
    Ok(0)
}

/// Records one observed barcode.
///
/// The `umis` ranking needs the distinct UMIs per cell, which are only turned
/// into counts after the scan; the `reads` ranking counts directly. Either
/// way the *first-seen order* is what breaks later ties, so `cell_umis` keeps
/// its own order alongside the set.
fn record(
    counts: &mut CellCounts,
    cell_umis: &mut CellUmiSets,
    rank_by: RankBy,
    cell: &str,
    umi: &str,
) {
    match rank_by {
        RankBy::Umis => cell_umis.observe(cell.as_bytes(), umi),
        RankBy::Reads => counts.observe(cell.as_bytes()),
    }
}

/// Distinct UMIs per cell, in the order the cells were first seen.
#[derive(Debug, Default)]
struct CellUmiSets {
    order: Vec<Vec<u8>>,
    sets: BTreeMap<Vec<u8>, std::collections::BTreeSet<String>>,
}

impl CellUmiSets {
    fn observe(&mut self, cell: &[u8], umi: &str) {
        if !self.sets.contains_key(cell) {
            self.order.push(cell.to_vec());
            self.sets
                .insert(cell.to_vec(), std::collections::BTreeSet::new());
        }
        if let Some(set) = self.sets.get_mut(cell) {
            set.insert(umi.to_string());
        }
    }

    fn iter(&self) -> impl Iterator<Item = (&Vec<u8>, usize)> + '_ {
        self.order.iter().map(move |c| {
            (
                c,
                self.sets.get(c).map_or(0, std::collections::BTreeSet::len),
            )
        })
    }
}

fn value_arg(name: &'static str, long: &'static str, help: &'static str) -> clap::Arg {
    clap::Arg::new(name).long(long).num_args(1).help(help)
}

fn parse_num<T: std::str::FromStr>(
    matches: &clap::ArgMatches,
    name: &str,
    label: &str,
) -> Result<Option<T>>
where
    T::Err: std::fmt::Display,
{
    matches
        .get_one::<String>(name)
        .map(|v| {
            v.parse::<T>()
                .map_err(|e| PortError::value(format!("{label}: {e}")))
        })
        .transpose()
}

fn io_error(e: std::io::Error) -> PortError {
    PortError::io(Path::new("<output>"), e)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(args: &[&str]) -> Result<i32> {
        let args: Vec<String> = args.iter().map(|s| (*s).to_string()).collect();
        main(&args, &["whitelist".to_string()])
    }

    #[test]
    fn about_names_the_command() {
        assert!(ABOUT.starts_with("whitelist - "));
    }

    #[test]
    fn a_cell_pattern_is_required() {
        let err = run(&["--bc-pattern=NNN"]).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.contains("do not include any cell bases")));
    }

    #[test]
    fn expect_cells_is_incompatible_with_the_distance_method() {
        let err = run(&["--bc-pattern=CCCNNN", "--expect-cells=100"]).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.contains("Cannot use --expect-cells with 'distance' knee method")));
    }

    #[test]
    fn expect_cells_and_cell_number_are_mutually_exclusive() {
        let err = run(&[
            "--bc-pattern=CCCNNN",
            "--knee-method=density",
            "--expect-cells=100",
            "--set-cell-number=50",
        ])
        .unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.contains("Cannot supply both --expect-cells and --cell-number")));
    }

    #[test]
    fn an_unknown_ranking_method_is_rejected() {
        let err = run(&["--bc-pattern=CCCNNN", "--method=magic", "-I", "/dev/null"]).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.contains("--method must be reads or umis")));
    }
}
