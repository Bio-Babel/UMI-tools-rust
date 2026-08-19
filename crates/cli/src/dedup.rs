//! `umi_tools dedup` — keep one read per original molecule
//! (Python source: `umi_tools/dedup.py`).
//!
//! Groups reads by alignment position and UMI, merges UMIs that are within an
//! edit distance of each other (they are sequencing errors of one true UMI),
//! and writes out one representative read per group. That representative is
//! chosen by fewest mapping locations, then highest mapping quality, then at
//! random among the reads still tied.

use std::collections::BTreeMap;
use std::path::Path;

use umi_tools_core::alignment::open_output_alignments;
use umi_tools_core::bundles::{BundleFlags, BundleItem};
use umi_tools_core::network::{ClusterMethod, ReadDeduplicator};
use umi_tools_core::options::validate_sam_options;
use umi_tools_core::stats;
use umi_tools_core::umi_methods::{get_average_umi_distance, RandomReadGenerator};
use umi_tools_core::utilities::{self, output_names_and_formats};
use umi_tools_core::whitelist::{get_user_defined_barcodes, ErrorCorrectionMode};
use umi_tools_core::{PortError, Result};

use crate::common;
use crate::documentation;
use crate::options::{self, OptionGroups};

/// The subcommand's short description.
pub const ABOUT: &str = "dedup - Deduplicate reads using UMI and mapping coordinates";

/// The SAM tags `detect_bam_features` looks for
/// (Python source: `tags = ["NH", "X0", "XT"]`).
pub const MULTIMAPPING_TAGS: [&str; 3] = ["NH", "X0", "XT"];

/// Which multimapping tags the first `n_entries` mapped reads all carry
/// (Python source: `dedup.detect_bam_features`).
///
/// A tag is reported as available only if **every** examined mapped read has
/// it, which is what makes `--multimapping-detection-method` safe to use: a
/// tag present on only some reads would silently change which read wins a tie.
///
/// Python examines `n_entries + 1` records, because its loop breaks on
/// `n > n_entries` after having already processed record `n`; that off-by-one
/// is reproduced.
#[must_use]
pub fn detect_bam_features(
    records: &[noodles_sam::alignment::RecordBuf],
    n_entries: usize,
) -> BTreeMap<String, bool> {
    let mut available: BTreeMap<String, bool> = MULTIMAPPING_TAGS
        .iter()
        .map(|t| ((*t).to_string(), true))
        .collect();

    for (n, record) in records.iter().enumerate() {
        if n > n_entries {
            break;
        }
        if record.flags().is_unmapped() {
            continue;
        }
        for tag in MULTIMAPPING_TAGS {
            let bytes = tag.as_bytes();
            let key = [bytes[0], bytes[1]];
            if umi_tools_core::sam::tag_text(record, key).is_none() {
                available.insert(tag.to_string(), false);
            }
        }
    }
    available
}

/// Runs `umi_tools dedup` (Python source: `dedup.main`).
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Value`] for an invalid option combination, and
/// [`umi_tools_core::PortError::Alignment`] or [`umi_tools_core::PortError::Io`] for input or output problems.
pub fn main(args: &[String], argv: &[String]) -> Result<i32> {
    if crate::wants_help(args) {
        // `optparse`'s layout, reproduced verbatim; see `crate::help`.
        print!("{}", crate::help::DEDUP);
        return Ok(0);
    }

    let started = std::time::Instant::now();
    let groups = OptionGroups {
        dedup_count_sam: true,
        sbcram: true,
        ..OptionGroups::default()
    };

    let command = options::build_parser("dedup", ABOUT, groups)
        .after_help(format!(
            "{}{}{}{}",
            documentation::GENERIC_DOCSTRING_GDC,
            documentation::GROUP_DEDUP_GENERIC_OPTIONS,
            documentation::GENERIC_DOCSTRING_SBCRAM_INPUT,
            documentation::GENERIC_DOCSTRING_SBCRAM_OUTPUT
        ))
        .arg(
            clap::Arg::new("stats")
                .overrides_with("stats")
                .long("output-stats")
                .num_args(1)
                .help("prefix for the statistics output files"),
        );

    let matches = match command
        .try_get_matches_from(std::iter::once("dedup".to_string()).chain(args.iter().cloned()))
    {
        Ok(m) => m,
        Err(e) => {
            e.print().ok();
            return Ok(i32::from(e.use_stderr()));
        }
    };

    let io = options::io_options(&matches)?;
    let run_id = options::start(&io, argv)?;

    let mut bundle_options = options::bundle_options(&matches)?;
    let dedup_options = options::dedup_options(&matches)?;
    let mut validation = options::sam_validation_input(&matches, argv)?;
    for warning in validate_sam_options(&mut validation, false)? {
        utilities::warn(&warning);
    }
    bundle_options.unmapped_reads = validation.unmapped_reads;

    let stats_prefix = matches.get_one::<String>("stats").map(Path::new);
    if stats_prefix.is_some() && bundle_options.ignore_umi {
        return Err(PortError::value(
            "'--output-stats' and '--ignore-umi' options cannot be used together",
        ));
    }

    let mut opened = common::open_input(&io)?;
    let records = opened.read_all()?;
    let header = opened.header.clone();

    if let Some(method) = bundle_options.detection_method {
        let features = detect_bam_features(&records, 1000);
        let wanted = match method {
            umi_tools_core::options::DetectionMethod::Nh => "NH",
            umi_tools_core::options::DetectionMethod::X0 => "X0",
            umi_tools_core::options::DetectionMethod::Xt => "XT",
        };
        if !features.get(wanted).copied().unwrap_or(false) {
            let available: Vec<&str> = features
                .iter()
                .filter(|(_, present)| **present)
                .map(|(t, _)| t.as_str())
                .collect();
            return Err(PortError::value(if available.is_empty() {
                "There are no bam tags available to detect multimapping. \
                 Do not set --multimapping-detection-method"
                    .to_string()
            } else {
                format!(
                    "The chosen method of detection for multimapping ({wanted}) will not work \
                     with this bam. Multimapping can be detected for this bam using any of the \
                     following: {}",
                    available.join(",")
                )
            }));
        }
    }

    // Python names the eventual output "-" when it is standard output, and
    // still sorts: `output_names_and_formats` routes the records through a
    // temporary file whenever `--no-sort-output` is absent, then sorts that
    // file to the pipe. Skipping the sort there yields the right reads in the
    // wrong order, which the package's md5-comparing test suite rejects.
    let to_stdout = io.stdout.is_none();
    let eventual = io
        .stdout
        .as_ref()
        .map_or_else(|| "-".to_string(), |p| p.to_string_lossy().into_owned());
    let names = output_names_and_formats(
        &eventual,
        io.out_sam,
        io.out_format,
        io.no_sort_output,
        io.tmpdir.as_deref(),
    )?;

    // With sorting on, `out_name` is a temporary file whichever way the
    // eventual output goes; only the unsorted, straight-to-stdout case writes
    // to the pipe directly.
    let unsorted_to_stdout = to_stdout && names.sorted_out_name.is_none();
    let raw_writer = if unsorted_to_stdout {
        umi_tools_core::alignment::open_output_alignments_to_stdout(
            names.out_format,
            &header,
            io.reference_filename.as_deref(),
        )?
    } else {
        open_output_alignments(
            Path::new(&names.out_name),
            names.out_format,
            &header,
            io.reference_filename.as_deref(),
        )?
    };

    // With `--paired`, each retained read 1 must be accompanied by its read 2,
    // which `TwoPassPairWriter` finds in a second pass over the same records.
    let reference_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
        .collect();
    let mut writer = if bundle_options.paired {
        Writer::Paired(Box::new(umi_tools_core::sam::TwoPassPairWriter::new(
            &records,
            reference_names,
            raw_writer,
        )))
    } else {
        Writer::Single(Box::new(raw_writer))
    };

    let mut deduplicator = if dedup_options.filter_umi {
        let path = dedup_options.umi_whitelist.as_deref().ok_or_else(|| {
            PortError::value("must provide a UMI whitelist (--umi-whitelist) if using --filter-umi")
        })?;
        let (whitelist, _) = get_user_defined_barcodes(
            path,
            dedup_options.umi_whitelist_paired.as_deref(),
            ErrorCorrectionMode::None,
            1,
        )?;
        ReadDeduplicator::with_umi_whitelist(dedup_options.method, whitelist)
    } else {
        ReadDeduplicator::new(dedup_options.method)
    };

    let mut bundler = common::build_bundler(&matches, &bundle_options, BundleFlags::default())?;
    let getter = common::barcode_getter(&matches, bundle_options.per_cell)?;

    // The null distribution for --output-stats needs UMI frequencies over the
    // whole file. Python makes a second pass over the BAM for this; the
    // records are already in memory here, so one scan suffices.
    let mut random_reads = if stats_prefix.is_none() {
        None
    } else {
        // Python's `random_read_generator` counts UMIs into a
        // `defaultdict(int)`, so the table is in first-seen order — and
        // `numpy.random.choice` indexes into that order, so it has to be
        // preserved rather than sorted.
        let mut order: Vec<Vec<u8>> = Vec::new();
        let mut frequencies: BTreeMap<Vec<u8>, u64> = BTreeMap::new();
        for record in umi_tools_core::alignment::fetch(&records, false) {
            let flags = record.flags();
            // Python skips unmapped reads and read 2s, and swallows a
            // missing tag rather than failing.
            if flags.is_unmapped() || (flags.is_segmented() && flags.is_last_segment()) {
                continue;
            }
            if let Ok(barcode) = getter.get(record) {
                if let Some(count) = frequencies.get_mut(&barcode.umi) {
                    *count += 1;
                } else {
                    order.push(barcode.umi.clone());
                    frequencies.insert(barcode.umi, 1);
                }
            }
        }
        let table: Vec<(Vec<u8>, u64)> = order
            .into_iter()
            .map(|umi| {
                let count = frequencies.get(&umi).copied().unwrap_or(0);
                (umi, count)
            })
            .collect();
        Some(RandomReadGenerator::new(
            &table,
            u32::try_from(bundle_options.random_seed.unwrap_or(0) & 0xffff_ffff).unwrap_or(0),
        )?)
    };

    let mut n_input = 0u64;
    let mut n_output = 0u64;
    let mut input_reported = 0u64;
    let mut output_reported = 0u64;

    let mut stats_pre: Vec<(Vec<u8>, u64)> = Vec::new();
    let mut stats_post: Vec<(Vec<u8>, u64)> = Vec::new();
    let mut pre_cluster: Vec<f64> = Vec::new();
    let mut post_cluster: Vec<f64> = Vec::new();
    let mut pre_null: Vec<f64> = Vec::new();
    let mut post_null: Vec<f64> = Vec::new();

    let mut process = |items: Vec<BundleItem>,
                       writer: &mut Writer<'_>,
                       deduplicator: &mut ReadDeduplicator|
     -> Result<()> {
        for item in items {
            let BundleItem::Bundle { bundle, .. } = item else {
                continue;
            };
            let counts: umi_tools_core::network::CountMap = bundle
                .iter()
                .map(|(umi, entry)| (umi.clone(), entry.count))
                .collect();
            n_input += counts.values().sum::<u64>();
            while n_output >= output_reported + 100_000 {
                output_reported += 100_000;
                utilities::info(&format!("Written out {output_reported} reads"));
            }
            while n_input >= input_reported + 1_000_000 {
                input_reported += 1_000_000;
                utilities::info(&format!("Parsed {input_reported} input reads"));
            }

            if bundle_options.ignore_umi {
                for entry in bundle.values() {
                    if let Some(read) = entry.reads.one() {
                        writer.write(read)?;
                        n_output += 1;
                    }
                }
                continue;
            }

            // The bundle's own insertion order, so the groups — and so the
            // reads written out — come out in Python's order.
            let order = umi_tools_core::bundles::insertion_order(&bundle);
            let result = deduplicator.dedup_ordered(&order, &counts, dedup_options.threshold)?;
            if result.final_umis.is_empty() {
                continue;
            }
            for umi in &result.final_umis {
                if let Some(read) = bundle.get(umi).and_then(|e| e.reads.one()) {
                    writer.write(read)?;
                    n_output += 1;
                }
            }

            if let Some(generator) = random_reads.as_mut() {
                let bundle_umis: Vec<Vec<u8>> = bundle.keys().cloned().collect();
                pre_cluster.push(get_average_umi_distance(&bundle_umis));
                let random = generator.get_umis(bundle_umis.len());
                pre_null.push(get_average_umi_distance(&random));
                for (umi, entry) in &bundle {
                    stats_pre.push((umi.clone(), entry.count));
                }

                for (umi, count) in result.final_umis.iter().zip(result.umi_counts.iter()) {
                    stats_post.push((umi.clone(), *count));
                }
                post_cluster.push(get_average_umi_distance(&result.final_umis));
                let random = generator.get_umis(result.final_umis.len());
                post_null.push(get_average_umi_distance(&random));
            }
        }
        Ok(())
    };

    // `dedup` always uses plain `fetch()`, which skips the unplaced block.
    for record in umi_tools_core::alignment::fetch(&records, false) {
        let reference_name = opened.reference_name(record);
        let items = bundler.push(record.clone(), reference_name.as_deref())?;
        process(items, &mut writer, &mut deduplicator)?;
    }
    process(bundler.finish(), &mut writer, &mut deduplicator)?;
    writer.finish()?;

    if let (Some(sorted_name), Some(sorted_format)) =
        (names.sorted_out_name.as_ref(), names.sorted_out_format)
    {
        umi_tools_core::alignment::sort_output_to(
            Path::new(&names.out_name),
            (!to_stdout).then(|| Path::new(sorted_name.as_str())),
            names.out_format,
            sorted_format,
            io.reference_filename.as_deref(),
        )?;
    }

    if let Some(prefix) = stats_prefix {
        write_stats(
            prefix,
            &stats_pre,
            &stats_post,
            &pre_cluster,
            &post_cluster,
            &pre_null,
            &post_null,
            dedup_options.method,
        )?;
    }

    common::log_read_events(bundler.read_events());
    utilities::info(&format!("Number of reads out: {n_output}"));
    if !bundle_options.ignore_umi {
        let clusterer = deduplicator.clusterer();
        common::log_cluster_stats(
            clusterer.positions(),
            clusterer.total_umis_per_position(),
            clusterer.max_umis_per_position(),
        );
    }
    if dedup_options.filter_umi {
        let discarded = deduplicator
            .umi_whitelist_counts()
            .get(umi_tools_core::network::NON_WHITELIST_UMI)
            .copied()
            .unwrap_or(0);
        utilities::info(&format!(
            "{discarded} UMIs were in a group where the top UMI was not a whitelist UMI and \
             were therefore discarded"
        ));
    }

    options::stop(&io, &run_id, started.elapsed().as_secs())?;
    Ok(0)
}

/// Where deduplicated reads go.
///
/// `--paired` routes them through [`umi_tools_core::sam::TwoPassPairWriter`],
/// which also emits each retained read's mate; everything else writes
/// directly. Python makes the same choice by rebinding `outfile`.
enum Writer<'a> {
    /// Single-end: write the record and nothing else.
    Single(Box<umi_tools_core::alignment::OpenedAlignmentWriter>),
    /// Paired-end: also owe, and later emit, each read's mate.
    Paired(Box<umi_tools_core::sam::TwoPassPairWriter<'a>>),
}

impl Writer<'_> {
    fn write(&mut self, record: &noodles_sam::alignment::RecordBuf) -> Result<()> {
        match self {
            Writer::Single(w) => w.write(record),
            Writer::Paired(w) => w.write(record, false),
        }
    }

    fn finish(self) -> Result<()> {
        match self {
            Writer::Single(w) => w.finish(),
            Writer::Paired(w) => w.close(),
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn write_stats(
    prefix: &Path,
    stats_pre: &[(Vec<u8>, u64)],
    stats_post: &[(Vec<u8>, u64)],
    pre_cluster: &[f64],
    post_cluster: &[f64],
    pre_null: &[f64],
    post_null: &[f64],
    method: ClusterMethod,
) -> Result<()> {
    let sidecar = |suffix: &str| -> std::path::PathBuf {
        prefix.with_file_name(format!(
            "{}{suffix}",
            prefix
                .file_name()
                .map_or_else(String::new, |n| n.to_string_lossy().into_owned())
        ))
    };

    stats::write_per_umi_per_position(
        stats_pre,
        stats_post,
        &sidecar("_per_umi_per_position.tsv"),
    )?;
    stats::write_per_umi(
        &stats::aggregate_stats(stats_pre),
        &stats::aggregate_stats(stats_post),
        &sidecar("_per_umi.tsv"),
    )?;
    let table = stats::bin_edit_distances(pre_cluster, post_cluster, pre_null, post_null)?;
    stats::write_edit_distance(&table, method.as_str(), &sidecar("_edit_distance.tsv"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::alignment::record::data::field::Tag;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record_buf::data::field::Value;
    use noodles_sam::alignment::RecordBuf;

    fn record(tags: &[(&[u8; 2], i32)], unmapped: bool) -> RecordBuf {
        let mut r = RecordBuf::default();
        *r.flags_mut() = if unmapped {
            Flags::UNMAPPED
        } else {
            Flags::empty()
        };
        for (tag, value) in tags {
            r.data_mut().insert(Tag::from(**tag), Value::Int32(*value));
        }
        r
    }

    #[test]
    fn about_names_the_command() {
        assert!(ABOUT.starts_with("dedup - "));
    }

    #[test]
    fn feature_detection_requires_every_read_to_carry_the_tag() {
        let records = vec![record(&[(b"NH", 1)], false), record(&[(b"NH", 1)], false)];
        let features = detect_bam_features(&records, 1000);
        assert!(features["NH"]);
        assert!(!features["X0"]);
        assert!(!features["XT"]);

        // One read without NH is enough to make it unavailable.
        let records = vec![record(&[(b"NH", 1)], false), record(&[], false)];
        assert!(!detect_bam_features(&records, 1000)["NH"]);
    }

    #[test]
    fn feature_detection_skips_unmapped_reads() {
        let records = vec![record(&[(b"NH", 1)], false), record(&[], true)];
        assert!(detect_bam_features(&records, 1000)["NH"]);
    }

    #[test]
    fn feature_detection_honours_the_entry_limit() {
        // Python breaks on `n > n_entries`, so it examines n_entries + 1 reads.
        let mut records = vec![record(&[(b"NH", 1)], false); 3];
        records.push(record(&[], false));
        assert!(detect_bam_features(&records, 2)["NH"]);
        assert!(!detect_bam_features(&records, 3)["NH"]);
    }

    #[test]
    fn stats_and_ignore_umi_are_mutually_exclusive() {
        let args: Vec<String> = ["--output-stats=x", "--ignore-umi", "-I", "/dev/null"]
            .iter()
            .map(|s| (*s).to_string())
            .collect();
        let err = main(&args, &["dedup".to_string()]).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.contains("cannot be used together")));
    }
}
