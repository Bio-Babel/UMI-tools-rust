//! `umi_tools group` — tag reads with their group, keeping all of them
//! (Python source: `umi_tools/group.py`).
//!
//! Does the same grouping as [`crate::dedup`] but discards nothing: each read
//! is tagged with the group it was assigned to (`UG`) and the group's inferred
//! true UMI (`BX`), and optionally described in a flat file. Useful for
//! inspecting what deduplication *would* do before committing to it.

use std::io::Write;
use std::path::Path;

use umi_tools_core::alignment::AlignmentSink;
use umi_tools_core::bundles::{BundleFlags, BundleItem, BundleReads};
use umi_tools_core::network::UmiClusterer;
use umi_tools_core::options::{validate_sam_options, UnmappedReads};
use umi_tools_core::record::TagValue;
use umi_tools_core::sam::get_read_position;
use umi_tools_core::utilities;
use umi_tools_core::{PortError, Result};

use crate::common;
use crate::documentation;
use crate::options::{self, OptionGroups};

/// The subcommand's short description.
pub const ABOUT: &str = "group - Group reads based on their UMI and mapping coordinates";

/// The columns of the `--group-out` flat file
/// (Python source: the header `group.main` writes).
pub const GROUP_TSV_HEADER: &str = "read_id\tcontig\tposition\tgene\tumi\tumi_count\t\
                                    final_umi\tfinal_umi_count\tunique_id";

/// Runs `umi_tools group` (Python source: `group.main`).
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Value`] for an invalid option combination, and
/// [`umi_tools_core::PortError::Alignment`] or [`umi_tools_core::PortError::Io`] for input or output problems.
pub fn main(args: &[String], argv: &[String]) -> Result<i32> {
    if crate::wants_help(args) {
        // `optparse`'s layout, reproduced verbatim; see `crate::help`.
        print!("{}", crate::help::GROUP);
        return Ok(0);
    }

    let started = std::time::Instant::now();
    let groups = OptionGroups {
        group_sam: true,
        sbcram: true,
        ..OptionGroups::default()
    };

    let command = options::build_parser("group", ABOUT, groups)
        .after_help(format!(
            "{}{}{}{}",
            documentation::GENERIC_DOCSTRING_GDC,
            documentation::GROUP_DEDUP_GENERIC_OPTIONS,
            documentation::GENERIC_DOCSTRING_SBCRAM_INPUT,
            documentation::GENERIC_DOCSTRING_SBCRAM_OUTPUT
        ))
        .arg(
            clap::Arg::new("tsv")
                .overrides_with("tsv")
                .long("group-out")
                .num_args(1)
                .help("outfile name for the file mapping read id to read group"),
        )
        .arg(
            clap::Arg::new("output_bam")
                .overrides_with("output_bam")
                .long("output-bam")
                .action(clap::ArgAction::SetTrue)
                .help("output a bam file with read groups tagged using the UG tag"),
        )
        .arg(
            clap::Arg::new("umi_group_tag")
                .overrides_with("umi_group_tag")
                .long("umi-group-tag")
                .num_args(1)
                .help("tag for the outputted umi group"),
        );

    let matches = match command
        .try_get_matches_from(std::iter::once("group".to_string()).chain(args.iter().cloned()))
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
    for warning in validate_sam_options(&mut validation, true)? {
        utilities::warn(&warning);
    }
    // `validate_sam_options` migrates the legacy `--output-unmapped` flag onto
    // `unmapped_reads`; the bundler reads its own copy, so the migrated value
    // has to come back across.
    bundle_options.unmapped_reads = validation.unmapped_reads;

    let output_bam = matches.get_flag("output_bam");
    if io.stdout.is_some() && !output_bam {
        return Err(PortError::assertion(
            "To output a bam you must include --output-bam option",
        ));
    }

    let umi_group_tag: [u8; 2] =
        matches
            .get_one::<String>("umi_group_tag")
            .map_or(Ok(*b"BX"), |v| {
                let b = v.as_bytes();
                if b.len() == 2 {
                    Ok([b[0], b[1]])
                } else {
                    Err(PortError::value(format!(
                        "SAM tags are two characters; got '{v}'"
                    )))
                }
            })?;

    let mut opened = common::open_input(&io)?;
    let header = opened.header.clone();
    let reference_names = common::reference_names(&header);

    // Python: `output_unmapped = options.unmapped_reads in ["use", "output"]`,
    // which also decides whether the reader runs to EOF.
    let output_unmapped = matches!(
        bundle_options.unmapped_reads,
        UnmappedReads::Use | UnmappedReads::Output
    );

    let mut writer: Option<AlignmentSink> = if output_bam || io.out_format.is_some() {
        Some(common::open_output(&io, &opened)?)
    } else {
        None
    };

    let mut mapping_out = matches
        .get_one::<String>("tsv")
        .map(|p| utilities::create_file(Path::new(p), false))
        .transpose()?;
    if let Some(f) = mapping_out.as_mut() {
        writeln!(f, "{GROUP_TSV_HEADER}").map_err(common::io_error)?;
    }

    let mut bundler = common::build_bundler(
        &matches,
        &bundle_options,
        BundleFlags {
            all_reads: true,
            summarise_reads: writer.is_none(),
            return_read2: true,
            return_unmapped: output_unmapped,
            ..BundleFlags::default()
        },
    )?;
    let input = common::input(
        &io,
        validation.gene_transcript_map.as_deref(),
        &header,
        output_unmapped,
    )?;
    if let common::Input::Meta(map) = &input {
        bundler.with_metacontigs(map.clone());
    }
    let mut clusterer = UmiClusterer::new(dedup_options.method);

    let mut n_input = 0u64;
    let mut n_output = 0u64;
    let mut unique_id = 0u64;
    let mut input_reported = 0u64;
    let mut output_reported = 0u64;

    let mut process = |items: &mut dyn Iterator<Item = BundleItem>,
                       clusterer: &mut UmiClusterer,
                       writer: &mut Option<AlignmentSink>,
                       mapping_out: &mut Option<umi_tools_core::utilities::OutputFile>|
     -> Result<()> {
        for item in items {
            match item {
                BundleItem::Single(read) => {
                    n_input += 1;
                    if let Some(w) = writer.as_mut() {
                        w.write(read)?;
                    }
                    n_output += 1;
                }
                BundleItem::Bundle { mut bundle, key } => {
                    let counts: umi_tools_core::network::CountMap = bundle
                        .iter()
                        .map(|(umi, entry)| (umi.clone(), entry.count))
                        .collect();
                    n_input += counts.values().sum::<u64>();
                    while n_output >= output_reported + 10_000 {
                        output_reported += 10_000;
                        utilities::info(&format!("Written out {output_reported} reads"));
                    }
                    while n_input >= input_reported + 1_000_000 {
                        input_reported += 1_000_000;
                        utilities::info(&format!("Parsed {input_reported} input reads"));
                    }

                    let order = umi_tools_core::bundles::insertion_order(&bundle);
                    let umi_groups =
                        clusterer.cluster_ordered(&order, &counts, dedup_options.threshold)?;
                    let gene = match (bundle_options.per_gene, bundle_options.per_contig) {
                        (true, true) => None,
                        (true, false) => Some(key.gene.as_deref().unwrap_or("NA")),
                        (false, _) => Some("NA"),
                    };
                    for umi_group in &umi_groups {
                        let Some(top_umi) = umi_group.first() else {
                            continue;
                        };
                        let group_count: u64 = umi_group
                            .iter()
                            .map(|u| counts.get(u).copied().unwrap_or(0))
                            .sum();

                        for umi in umi_group {
                            let Some(entry) = bundle.get_mut(umi) else {
                                continue;
                            };
                            let row = match mapping_out {
                                Some(_) => Some(GroupRow {
                                    gene,
                                    umi: utilities::decode_utf8(umi)?,
                                    umi_count: counts.get(umi).copied().unwrap_or(0),
                                    final_umi: utilities::decode_utf8(top_umi)?,
                                    final_umi_count: group_count,
                                    unique_id,
                                }),
                                None => None,
                            };
                            let contig = |id: Option<usize>| {
                                id.and_then(|id| reference_names.get(id))
                                    .map_or("", String::as_str)
                            };
                            match std::mem::take(&mut entry.reads) {
                                BundleReads::Many(reads) => {
                                    for mut read in reads {
                                        if let (Some(f), Some(row)) = (mapping_out.as_mut(), &row) {
                                            let position = get_read_position(
                                                &read,
                                                bundle_options.soft_clip_threshold,
                                            )?
                                            .pos;
                                            row.write(
                                                f,
                                                read.name().unwrap_or_default(),
                                                contig(read.reference_id()),
                                                position,
                                            )?;
                                        }
                                        if let Some(w) = writer.as_mut() {
                                            let tag_error = |e: std::io::Error| {
                                                PortError::alignment(e.to_string())
                                            };
                                            read.set_tag(
                                                *b"UG",
                                                TagValue::Int(
                                                    i64::try_from(unique_id).unwrap_or(i64::MAX),
                                                ),
                                            )
                                            .map_err(tag_error)?;
                                            read.set_tag(umi_group_tag, TagValue::String(top_umi))
                                                .map_err(tag_error)?;
                                            w.write(read)?;
                                        }
                                        n_output += 1;
                                    }
                                }
                                BundleReads::Summaries(reads) => {
                                    for read in reads {
                                        if let (Some(f), Some(row)) = (mapping_out.as_mut(), &row) {
                                            row.write(
                                                f,
                                                &read.name,
                                                contig(read.reference_id),
                                                read.position,
                                            )?;
                                        }
                                        n_output += 1;
                                    }
                                }
                                BundleReads::One(_) | BundleReads::CountOnly => {}
                            }
                        }
                        unique_id += 1;
                    }
                }
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
            &mut clusterer,
            &mut writer,
            &mut mapping_out,
        )
    })?;
    process(
        &mut bundler.finish(),
        &mut clusterer,
        &mut writer,
        &mut mapping_out,
    )?;

    if let Some(w) = writer.take() {
        w.finish()?;
    }
    if let Some(f) = mapping_out.take() {
        f.finish()?;
    }

    common::log_read_events(bundler.read_events());
    utilities::info(&format!(
        "Number of reads out: {n_output}, Number of groups: {unique_id}"
    ));
    common::log_cluster_stats(
        clusterer.positions(),
        clusterer.total_umis_per_position(),
        clusterer.max_umis_per_position(),
    );

    options::stop(&io, &run_id, started.elapsed().as_secs())?;
    Ok(0)
}

struct GroupRow<'a> {
    gene: Option<&'a str>,
    umi: String,
    umi_count: u64,
    final_umi: String,
    final_umi_count: u64,
    unique_id: u64,
}

impl GroupRow<'_> {
    fn write(
        &self,
        out: &mut utilities::OutputFile,
        name: &[u8],
        contig: &str,
        position: i64,
    ) -> Result<()> {
        writeln!(
            out,
            "{}\t{contig}\t{position}\t{}\t{}\t{}\t{}\t{}\t{}",
            String::from_utf8_lossy(name),
            self.gene.unwrap_or(contig),
            self.umi,
            self.umi_count,
            self.final_umi,
            self.final_umi_count,
            self.unique_id,
        )
        .map_err(common::io_error)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn about_names_the_command() {
        assert!(ABOUT.starts_with("group - "));
    }

    #[test]
    fn the_tsv_header_matches_pythons() {
        let columns: Vec<&str> = GROUP_TSV_HEADER.split('\t').collect();
        assert_eq!(
            columns,
            vec![
                "read_id",
                "contig",
                "position",
                "gene",
                "umi",
                "umi_count",
                "final_umi",
                "final_umi_count",
                "unique_id"
            ]
        );
    }

    #[test]
    fn an_output_file_requires_output_bam() {
        let args: Vec<String> = ["-S", "/tmp/out.bam", "-I", "/dev/null"]
            .iter()
            .map(|s| (*s).to_string())
            .collect();
        let err = main(&args, &["group".to_string()]).unwrap_err();
        assert!(matches!(err, PortError::Assertion(ref m)
            if m.contains("--output-bam")));
    }
}
