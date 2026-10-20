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

use umi_tools_core::alignment::{AlignmentSink, BamRecord};
use umi_tools_core::bundles::{BundleFlags, BundleItem, BundleReads};
use umi_tools_core::network::{ClusterMethod, ReadDeduplicator};
use umi_tools_core::options::validate_sam_options;
use umi_tools_core::sam::TwoPassPairWriter;
use umi_tools_core::stats;
use umi_tools_core::umi_methods::{get_average_umi_distance, RandomReadGenerator};
use umi_tools_core::utilities;
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
/// is reproduced. `records` is a fresh `fetch(until_eof=True)` of the input.
pub fn detect_bam_features(
    records: impl IntoIterator<Item = Result<BamRecord>>,
    n_entries: usize,
) -> Result<BTreeMap<String, bool>> {
    let mut available: BTreeMap<String, bool> = MULTIMAPPING_TAGS
        .iter()
        .map(|t| ((*t).to_string(), true))
        .collect();

    for (n, record) in records.into_iter().enumerate() {
        if n > n_entries {
            break;
        }
        let record = record?;
        if record.flags().is_unmapped() {
            continue;
        }
        let view = record.view();
        let data = view.data();
        for tag in MULTIMAPPING_TAGS {
            let bytes = tag.as_bytes();
            if data.get(&[bytes[0], bytes[1]]).is_none() {
                available.insert(tag.to_string(), false);
            }
        }
    }
    Ok(available)
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

    let source = common::input_source(&io)?;
    let mut opened = source.open()?;
    let header = opened.header.clone();
    let reference_names = common::reference_names(&header);

    if let Some(method) = bundle_options.detection_method {
        let features = detect_bam_features(source.open()?.records(true), 1000)?;
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

    // Sorted unless `--no-sort-output`, even to a pipe, as in Python.
    let sink = common::open_output(&io, &opened)?;
    let mut writer = if bundle_options.paired {
        Writer::Paired(Box::new(TwoPassPairWriter::new(source.clone(), sink)))
    } else {
        Writer::Single(Box::new(sink))
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
    let input = common::input(
        &io,
        validation.gene_transcript_map.as_deref(),
        &header,
        false,
    )?;
    if let common::Input::Meta(map) = &input {
        bundler.with_metacontigs(map.clone());
    }

    // Python's `random_read_generator`: its own pass over the input.
    let mut random_reads = if stats_prefix.is_none() {
        None
    } else {
        // Python's `random_read_generator` counts UMIs into a
        // `defaultdict(int)`, so the table is in first-seen order — and
        // `numpy.random.choice` indexes into that order, so it has to be
        // preserved rather than sorted.
        let mut frequencies: indexmap::IndexMap<Vec<u8>, u64> = indexmap::IndexMap::new();
        let whole = match &input {
            common::Input::Chrom(id) => common::Input::Chrom(*id),
            _ => common::Input::Fetch { until_eof: false },
        };
        common::for_each_record(&mut source.open()?, &whole, |record| {
            let flags = record.flags();
            // Python skips unmapped reads and read 2s, and swallows a
            // missing tag rather than failing.
            if flags.is_unmapped() || (flags.is_segmented() && flags.is_last_segment()) {
                return Ok(());
            }
            if let Ok(barcode) = getter.get(&record) {
                *frequencies.entry(barcode.umi).or_insert(0) += 1;
            }
            Ok(())
        })?;
        let table: Vec<(Vec<u8>, u64)> = frequencies.into_iter().collect();
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

    let mut process = |items: &mut dyn Iterator<Item = BundleItem>,
                       writer: &mut Writer,
                       deduplicator: &mut ReadDeduplicator|
     -> Result<()> {
        for item in items {
            let BundleItem::Bundle { mut bundle, .. } = item else {
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
                for (_, entry) in bundle {
                    if let BundleReads::One(read) = entry.reads {
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
                let reads = bundle.get_mut(umi).map(|e| std::mem::take(&mut e.reads));
                if let Some(BundleReads::One(read)) = reads {
                    writer.write(read)?;
                    n_output += 1;
                }
            }

            if let Some(generator) = random_reads.as_mut() {
                pre_cluster.push(get_average_umi_distance(&order));
                let random = generator.get_umis(order.len());
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

    common::for_each_record(&mut opened, &input, |record| {
        let reference_name = record
            .reference_id()
            .and_then(|id| reference_names.get(id))
            .map(String::as_str);
        process(
            &mut bundler.push(record, reference_name)?,
            &mut writer,
            &mut deduplicator,
        )
    })?;
    process(&mut bundler.finish(), &mut writer, &mut deduplicator)?;
    writer.finish()?;

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
/// `--paired` routes them through [`TwoPassPairWriter`], which also emits
/// each retained read's mate.
enum Writer {
    /// Single-end: write the record and nothing else.
    Single(Box<AlignmentSink>),
    /// Paired-end: also owe, and later emit, each read's mate.
    Paired(Box<TwoPassPairWriter>),
}

impl Writer {
    fn write(&mut self, record: BamRecord) -> Result<()> {
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

    fn bam(tags: &[(&[u8; 2], i32)], unmapped: bool) -> BamRecord {
        let mut r = RecordBuf::default();
        *r.flags_mut() = if unmapped {
            Flags::UNMAPPED
        } else {
            Flags::empty()
        };
        for (tag, value) in tags {
            r.data_mut().insert(Tag::from(**tag), Value::Int32(*value));
        }
        BamRecord::from_alignment_record(&noodles_sam::Header::default(), &r).unwrap()
    }

    #[test]
    fn about_names_the_command() {
        assert!(ABOUT.starts_with("dedup - "));
    }

    #[test]
    fn feature_detection_requires_every_read_to_carry_the_tag() {
        let records = vec![Ok(bam(&[(b"NH", 1)], false)), Ok(bam(&[(b"NH", 1)], false))];
        let features = detect_bam_features(records, 1000).unwrap();
        assert!(features["NH"]);
        assert!(!features["X0"]);
        assert!(!features["XT"]);

        // One read without NH is enough to make it unavailable.
        let records = vec![Ok(bam(&[(b"NH", 1)], false)), Ok(bam(&[], false))];
        assert!(!detect_bam_features(records, 1000).unwrap()["NH"]);
    }

    #[test]
    fn feature_detection_skips_unmapped_reads() {
        let records = vec![Ok(bam(&[(b"NH", 1)], false)), Ok(bam(&[], true))];
        assert!(detect_bam_features(records, 1000).unwrap()["NH"]);
    }

    #[test]
    fn feature_detection_honours_the_entry_limit() {
        // Python breaks on `n > n_entries`, so it examines n_entries + 1 reads.
        let records = || {
            let mut v: Vec<Result<BamRecord>> =
                (0..3).map(|_| Ok(bam(&[(b"NH", 1)], false))).collect();
            v.push(Ok(bam(&[], false)));
            v
        };
        assert!(detect_bam_features(records(), 2).unwrap()["NH"]);
        assert!(!detect_bam_features(records(), 3).unwrap()["NH"]);
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
