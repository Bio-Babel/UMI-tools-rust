//! Plumbing shared by the alignment-consuming subcommands
//! (Python source: the near-identical preambles of `dedup.main`,
//! `group.main` and `count.main`).
//!
//! The three commands open their input the same way, build the same barcode
//! getter from the same options, and set up the same bundler. Python repeats
//! that preamble in each module; keeping it in one place here means the three
//! `main` functions differ only where the commands genuinely differ.

use std::path::Path;

use umi_tools_core::alignment::{open_input_alignments, OpenedAlignment};
use umi_tools_core::bundles::{BundleFlags, GetBundles};
use umi_tools_core::options::{BundleOptions, IoOptions, UmiMethod};
use umi_tools_core::sam::{BarcodeGetter, TagConfig};
use umi_tools_core::utilities::{self, determine_format};
use umi_tools_core::{PortError, Result};

/// Builds the barcode getter for the configured `--extract-umi-method`
/// (Python source: the `functools.partial` chain in `get_bundles.__init__`).
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Value`] for an unknown method or a malformed tag.
pub fn barcode_getter(matches: &clap::ArgMatches, per_cell: bool) -> Result<BarcodeGetter> {
    let method: UmiMethod = matches
        .get_one::<String>("get_umi_method")
        .map_or(Ok(UmiMethod::ReadId), |v| v.parse())?;

    let tag = |name: &str, default: Option<[u8; 2]>| -> Result<Option<[u8; 2]>> {
        match matches.get_one::<String>(name) {
            None => Ok(default),
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
    };

    Ok(match method {
        UmiMethod::ReadId => BarcodeGetter::ReadId {
            cell_barcode: per_cell,
            sep: matches
                .get_one::<String>("umi_sep")
                .cloned()
                .unwrap_or_else(|| "_".to_string()),
        },
        UmiMethod::Tag => BarcodeGetter::Tag(TagConfig {
            cell_barcode: per_cell,
            umi_tag: tag("umi_tag", Some(*b"RX"))?.unwrap_or(*b"RX"),
            cell_tag: tag("cell_tag", None)?,
            umi_tag_split: matches.get_one::<String>("umi_tag_split").cloned(),
            umi_tag_delim: matches.get_one::<String>("umi_tag_delim").cloned(),
            cell_tag_split: Some(
                matches
                    .get_one::<String>("cell_tag_split")
                    .cloned()
                    .unwrap_or_else(|| "-".to_string()),
            ),
            cell_tag_delim: matches.get_one::<String>("cell_tag_delim").cloned(),
        }),
        UmiMethod::Umis => BarcodeGetter::Umis {
            cell_barcode: per_cell,
        },
    })
}

/// Opens the alignment named by `-I/--stdin`
/// (Python source: the `options.stdin` handling every alignment command
/// repeats, including its refusal to read an alignment from a pipe).
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Value`] when no input file was given — Python raises
/// `ValueError("Input on standard in not currently supported")` — and
/// [`umi_tools_core::PortError::Io`] or [`umi_tools_core::PortError::Alignment`] if the file cannot be read.
pub fn open_input(io: &IoOptions) -> Result<OpenedAlignment> {
    let path = io
        .stdin
        .as_deref()
        .ok_or_else(|| PortError::value("Input on standard in not currently supported"))?;
    let format = determine_format(&path.to_string_lossy(), io.in_sam, io.in_format);
    open_input_alignments(path, format, io.reference_filename.as_deref())
}

/// Builds a bundler over an opened alignment.
///
/// # Errors
///
/// Propagates [`barcode_getter`]'s errors.
pub fn build_bundler(
    matches: &clap::ArgMatches,
    options: &BundleOptions,
    flags: BundleFlags,
) -> Result<GetBundles> {
    let getter = barcode_getter(matches, options.per_cell)?;
    Ok(GetBundles::new(options, flags, getter))
}

/// Logs the read-event tallies the way `dedup` and `group` do
/// (Python source: the `"Reads: %s" % ", ".join(...)` line, ordered by
/// descending count).
pub fn log_read_events(events: &std::collections::BTreeMap<String, u64>) {
    let mut ordered: Vec<(&String, &u64)> = events.iter().collect();
    ordered.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
    let rendered: Vec<String> = ordered
        .iter()
        .map(|(name, count)| format!("{name}: {count}"))
        .collect();
    utilities::info(&format!("Reads: {}", rendered.join(", ")));
}

/// Logs the clusterer's per-position statistics
/// (Python source: the trailing block of `dedup.main` and `group.main`).
pub fn log_cluster_stats(positions: u64, total_umis: u64, max_umis: usize) {
    utilities::info(&format!(
        "Total number of positions deduplicated: {positions}"
    ));
    if positions > 0 {
        let mean = total_umis as f64 / positions as f64;
        utilities::info(&format!(
            "Mean number of unique UMIs per position: {mean:.2}"
        ));
        utilities::info(&format!(
            "Max. number of unique UMIs per position: {max_umis}"
        ));
    } else {
        utilities::warn("The BAM did not contain any valid reads/read pairs for deduplication");
    }
}

/// Turns an I/O error into a [`PortError`] tagged with a generic output path.
pub fn io_error(e: std::io::Error) -> PortError {
    PortError::io(Path::new("<output>"), e)
}
