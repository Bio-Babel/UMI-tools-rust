//! `umi_tools count` — count molecules per gene from an alignment
//! (Python source: `umi_tools/count.py`).
//!
//! Counts distinct molecules rather than reads: each gene's reads are grouped
//! by UMI, near-identical UMIs are merged as sequencing errors, and the number
//! of surviving groups is the count. Intended for protocols that fragment
//! after amplification (10x, inDrop, Drop-seq, SCRB-seq, CEL-seq2), where the
//! mapping coordinate is no longer informative — which is why `--per-gene` is
//! hard-coded on, as the Python docstring notes.

use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufWriter, Write};

use umi_tools_core::bundles::{BundleFlags, BundleItem};
use umi_tools_core::network::UmiClusterer;
use umi_tools_core::options::validate_sam_options;
use umi_tools_core::utilities;
use umi_tools_core::Result;

use crate::common;
use crate::documentation;
use crate::options::{self, OptionGroups};

/// The subcommand's short description.
pub const ABOUT: &str = "count - Count reads per gene from BAM using UMIs and mapping coordinates";

/// Runs `umi_tools count` (Python source: `count.main`).
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Value`] for an invalid option combination, and
/// [`umi_tools_core::PortError::Alignment`] or [`umi_tools_core::PortError::Io`] for input problems.
pub fn main(args: &[String], argv: &[String]) -> Result<i32> {
    if crate::wants_help(args) {
        // `optparse`'s layout, reproduced verbatim; see `crate::help`.
        print!("{}", crate::help::COUNT);
        return Ok(0);
    }

    let started = std::time::Instant::now();
    let groups = OptionGroups {
        group_dedup: false,
        dedup_count_sam: true,
        sbcram: true,
        exclude_output_format: true,
        ..OptionGroups::default()
    };

    let command = options::build_parser("count", ABOUT, groups)
        .after_help(format!(
            "{}{}",
            documentation::GENERIC_DOCSTRING_GDC,
            documentation::GENERIC_DOCSTRING_SBCRAM_INPUT
        ))
        .arg(
            clap::Arg::new("wide_format_cell_counts")
                .overrides_with("wide_format_cell_counts")
                .long("wide-format-cell-counts")
                .action(clap::ArgAction::SetTrue)
                .help("output the cell counts in a wide format (rows=genes, columns=cells)"),
        );

    let matches = match command
        .try_get_matches_from(std::iter::once("count".to_string()).chain(args.iter().cloned()))
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
    // `count` hard-codes per-gene counting, as its docstring says.
    bundle_options.per_gene = true;

    let mut validation = options::sam_validation_input(&matches, argv)?;
    validation.per_gene = true;
    for warning in validate_sam_options(&mut validation, false)? {
        utilities::warn(&warning);
    }

    let dedup_options = options::dedup_options(&matches)?;
    let per_cell = bundle_options.per_cell;

    let mut opened = common::open_input(&io)?;
    let records = opened.read_all()?;
    let mut bundler = common::build_bundler(
        &matches,
        &bundle_options,
        BundleFlags {
            only_count_reads: true,
            ..BundleFlags::default()
        },
    )?;

    let mut clusterer = UmiClusterer::new(dedup_options.method);
    // gene -> cell -> count. `None` is the no-cell-barcode case.
    let mut gene_counts: BTreeMap<String, BTreeMap<Option<String>, usize>> = BTreeMap::new();
    let mut cells: BTreeSet<Option<String>> = BTreeSet::new();
    let mut n_input = 0u64;
    let mut n_output = 0u64;
    let mut reported = 0u64;

    let handle = |items: Vec<BundleItem>,
                  clusterer: &mut UmiClusterer,
                  gene_counts: &mut BTreeMap<String, BTreeMap<Option<String>, usize>>,
                  cells: &mut BTreeSet<Option<String>>,
                  n_input: &mut u64,
                  n_output: &mut u64,
                  reported: &mut u64|
     -> Result<()> {
        for item in items {
            let BundleItem::Bundle { bundle, key } = item else {
                continue;
            };
            let counts: umi_tools_core::network::CountMap = bundle
                .iter()
                .map(|(umi, entry)| (umi.clone(), entry.count))
                .collect();
            *n_input += counts.values().sum::<u64>();
            while *n_input >= *reported + 1_000_000 {
                *reported += 1_000_000;
                utilities::info(&format!("Parsed {reported} input reads"));
            }

            let order = umi_tools_core::bundles::insertion_order(&bundle);
            let groups = clusterer.cluster_ordered(&order, &counts, dedup_options.threshold)?;
            let gene = key.gene.clone().unwrap_or_default();
            let cell = key
                .cell
                .as_ref()
                .map(|c| utilities::decode_utf8(c))
                .transpose()?;
            cells.insert(cell.clone());
            *gene_counts
                .entry(gene)
                .or_default()
                .entry(cell)
                .or_insert(0) += groups.len();
            *n_output += groups.len() as u64;
        }
        Ok(())
    };

    for record in umi_tools_core::alignment::fetch(&records, false) {
        let reference_name = opened.reference_name(record);
        let items = bundler.push(record.clone(), reference_name.as_deref())?;
        handle(
            items,
            &mut clusterer,
            &mut gene_counts,
            &mut cells,
            &mut n_input,
            &mut n_output,
            &mut reported,
        )?;
    }
    handle(
        bundler.finish(),
        &mut clusterer,
        &mut gene_counts,
        &mut cells,
        &mut n_input,
        &mut n_output,
        &mut reported,
    )?;

    let mut out = BufWriter::new(std::io::stdout());
    if per_cell {
        let cell_names: Vec<String> = cells
            .iter()
            .map(|c| c.clone().unwrap_or_default())
            .collect();
        if matches.get_flag("wide_format_cell_counts") {
            writeln!(out, "gene\t{}", cell_names.join("\t")).map_err(common::io_error)?;
            for (gene, by_cell) in &gene_counts {
                let row: Vec<String> = cells
                    .iter()
                    .map(|c| by_cell.get(c).copied().unwrap_or(0).to_string())
                    .collect();
                writeln!(out, "{gene}\t{}", row.join("\t")).map_err(common::io_error)?;
            }
        } else {
            writeln!(out, "gene\tcell\tcount").map_err(common::io_error)?;
            for (gene, by_cell) in &gene_counts {
                for (cell, count) in by_cell {
                    writeln!(out, "{gene}\t{}\t{count}", cell.clone().unwrap_or_default())
                        .map_err(common::io_error)?;
                }
            }
        }
    } else {
        writeln!(out, "gene\tcount").map_err(common::io_error)?;
        for (gene, by_cell) in &gene_counts {
            let total: usize = by_cell.values().sum();
            writeln!(out, "{gene}\t{total}").map_err(common::io_error)?;
        }
    }
    out.flush().map_err(common::io_error)?;

    for (event, count) in bundler.read_events() {
        utilities::info(&format!("{event}: {count}"));
    }
    utilities::info(&format!(
        "Number of (post deduplication) reads counted: {n_output}"
    ));

    options::stop(&io, &run_id, started.elapsed().as_secs())?;
    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use umi_tools_core::PortError;

    #[test]
    fn about_names_the_command() {
        assert!(ABOUT.starts_with("count - "));
    }

    #[test]
    fn per_gene_requires_a_gene_source() {
        // `count` forces --per-gene on, so it must be told where genes come from.
        let args: Vec<String> = vec!["-I".to_string(), "/dev/null".to_string()];
        let err = main(&args, &["count".to_string()]).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.contains("must supply --per-contig or --gene-tag")));
    }
}
