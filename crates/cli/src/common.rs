//! Plumbing shared by the alignment-consuming subcommands
//! (Python source: the near-identical preambles of `dedup.main`,
//! `group.main` and `count.main`).
//!
//! The three commands open their input the same way, build the same barcode
//! getter from the same options, and set up the same bundler. Python repeats
//! that preamble in each module; keeping it in one place here means the three
//! `main` functions differ only where the commands genuinely differ.

use std::path::Path;

use umi_tools_core::alignment::{
    open_input_alignments, AlignmentSink, AlignmentSource, OpenedAlignment,
};
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

/// Describes the `-I/--stdin` input, for side passes that reopen it.
pub fn input_source(io: &IoOptions) -> Result<AlignmentSource> {
    let path = io
        .stdin
        .clone()
        .ok_or_else(|| PortError::value("Input on standard in not currently supported"))?;
    let format = determine_format(&path.to_string_lossy(), io.in_sam, io.in_format);
    Ok(AlignmentSource {
        path,
        format,
        reference: io.reference_filename.clone(),
    })
}

/// Opens where `dedup` and `group` write their reads, coordinate-sorted
/// unless `--no-sort-output` (Python source: `Utilities.sort_output`).
pub fn open_output(io: &IoOptions, input: &OpenedAlignment) -> Result<AlignmentSink> {
    let eventual = io
        .stdout
        .as_ref()
        .map_or_else(|| "-".to_string(), |p| p.to_string_lossy().into_owned());
    AlignmentSink::open(
        io.stdout.as_deref(),
        determine_format(&eventual, io.out_sam, io.out_format),
        &input.header,
        input.header_text(),
        io.reference_filename.as_deref(),
        !io.no_sort_output,
        io.tmpdir.as_deref(),
    )
}

/// What the main pass reads (Python source: `inreads`).
#[derive(Debug, Clone)]
pub enum Input {
    /// `infile.fetch()`, or `fetch(until_eof=True)`.
    Fetch {
        /// Include the unplaced reads at the end.
        until_eof: bool,
    },
    /// `infile.fetch(reference=options.chrom)`: one reference, by index.
    Chrom(usize),
    /// `metafetcher(infile, metacontig2contig, "MC")`.
    Meta(umi_tools_core::sam::MetaContigMap),
}

/// Chooses the main pass's input as Python does: `--chrom`, then
/// `--gene-transcript-map`, then the whole file.
pub fn input(
    io: &IoOptions,
    gene_transcript_map: Option<&Path>,
    header: &noodles_sam::Header,
    until_eof: bool,
) -> Result<Input> {
    let names = reference_names(header);
    if let Some(chrom) = &io.chrom {
        return names
            .iter()
            .position(|n| n == chrom)
            .map(Input::Chrom)
            .ok_or_else(|| PortError::value(format!("invalid contig `{chrom}`")));
    }
    if let Some(map) = gene_transcript_map {
        let references = names.into_iter().collect();
        return Ok(Input::Meta(umi_tools_core::sam::get_metacontig2contig(
            &references,
            map,
        )?));
    }
    Ok(Input::Fetch { until_eof })
}

/// Feeds every record of `input` to `each`, in order.
pub fn for_each_record(
    opened: &mut OpenedAlignment,
    input: &Input,
    mut each: impl FnMut(umi_tools_core::BamRecord) -> Result<()>,
) -> Result<()> {
    match input {
        Input::Fetch { until_eof } => {
            for record in opened.records(*until_eof) {
                each(record?)?;
            }
        }
        Input::Chrom(id) => {
            for record in opened.fetch_reference(*id)? {
                each(record?)?;
            }
        }
        Input::Meta(map) => {
            for record in umi_tools_core::sam::metafetch(opened, map, *b"MC")? {
                each(record?)?;
            }
        }
    }
    Ok(())
}

/// The header's reference-sequence names, indexed by reference id.
#[must_use]
pub fn reference_names(header: &noodles_sam::Header) -> Vec<String> {
    header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
        .collect()
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

/// A text destination: `-S/--stdout`, or standard output.
pub enum TextOutput {
    /// The `-S` file, gzip-compressed if its name ends in `.gz`.
    File(Box<umi_tools_core::utilities::OutputFile>),
    /// Standard output.
    Stdout(std::io::BufWriter<std::io::Stdout>),
}

impl std::io::Write for TextOutput {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            TextOutput::File(f) => f.write(buf),
            TextOutput::Stdout(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            TextOutput::File(f) => f.flush(),
            TextOutput::Stdout(w) => w.flush(),
        }
    }
}

impl TextOutput {
    /// Flushes and closes the destination.
    pub fn finish(self) -> Result<()> {
        use std::io::Write as _;
        match self {
            TextOutput::File(f) => f.finish(),
            TextOutput::Stdout(mut w) => w.flush().map_err(io_error),
        }
    }
}

/// Opens `path`, or standard output if `None`.
pub fn open_text_output(path: Option<&Path>) -> Result<TextOutput> {
    match path {
        Some(p) => Ok(TextOutput::File(Box::new(utilities::create_file(
            p, false,
        )?))),
        None => Ok(TextOutput::Stdout(std::io::BufWriter::new(
            std::io::stdout(),
        ))),
    }
}

/// Turns an I/O error into a [`PortError`] tagged with a generic output path.
pub fn io_error(e: std::io::Error) -> PortError {
    PortError::io(Path::new("<output>"), e)
}
