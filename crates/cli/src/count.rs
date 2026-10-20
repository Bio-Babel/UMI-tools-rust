//! `umi_tools count` — count molecules per gene from an alignment
//! (Python source: `umi_tools/count.py`).
//!
//! Counts distinct molecules rather than reads: each gene's reads are grouped
//! by UMI, near-identical UMIs are merged as sequencing errors, and the number
//! of surviving groups is the count. Intended for protocols that fragment
//! after amplification (10x, inDrop, Drop-seq, SCRB-seq, CEL-seq2), where the
//! mapping coordinate is no longer informative — which is why `--per-gene` is
//! hard-coded on, as the Python docstring notes.

use std::collections::BTreeMap;
use std::io::Write;

use indexmap::IndexSet;

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
    let reference_names = common::reference_names(&opened.header);
    let mut bundler = common::build_bundler(
        &matches,
        &bundle_options,
        BundleFlags {
            only_count_reads: true,
            ..BundleFlags::default()
        },
    )?;

    let input = common::input(
        &io,
        validation.gene_transcript_map.as_deref(),
        &opened.header,
        false,
    )?;
    if let common::Input::Meta(map) = &input {
        bundler.with_metacontigs(map.clone());
    }

    let mut clusterer = UmiClusterer::new(dedup_options.method);
    let mut out = common::open_text_output(io.stdout.as_deref())?;
    // Without --per-cell Python writes a row per bundle, in bundle order.
    if !per_cell {
        writeln!(out, "gene\tcount").map_err(common::io_error)?;
    }
    // With --per-cell it keeps the last count per (gene, cell), sorted at the end.
    let mut cells: IndexSet<String> = IndexSet::new();
    let mut gene_counts: BTreeMap<String, Vec<(u32, usize)>> = BTreeMap::new();
    let mut n_input = 0u64;
    let mut n_output = 0u64;
    let mut reported = 0u64;

    let mut handle = |items: &mut dyn Iterator<Item = BundleItem>,
                      clusterer: &mut UmiClusterer,
                      out: &mut common::TextOutput|
     -> Result<()> {
        for item in items {
            let BundleItem::Bundle { bundle, key } = item else {
                continue;
            };
            let counts: umi_tools_core::network::CountMap = bundle
                .iter()
                .map(|(umi, entry)| (umi.clone(), entry.count))
                .collect();
            n_input += counts.values().sum::<u64>();
            while n_input >= reported + 1_000_000 {
                reported += 1_000_000;
                utilities::info(&format!("Parsed {reported} input reads"));
            }

            let order = umi_tools_core::bundles::insertion_order(&bundle);
            let groups = clusterer.cluster_ordered(&order, &counts, dedup_options.threshold)?;
            let gene = key.gene.unwrap_or_default();
            n_output += groups.len() as u64;
            if !per_cell {
                writeln!(out, "{gene}\t{}", groups.len()).map_err(common::io_error)?;
                continue;
            }
            let cell = match &key.cell {
                Some(c) => utilities::decode_utf8(c)?,
                None => String::new(),
            };
            let (cell, _) = cells.insert_full(cell);
            let cell = u32::try_from(cell)
                .map_err(|_| umi_tools_core::PortError::value("too many cell barcodes"))?;
            gene_counts
                .entry(gene)
                .or_default()
                .push((cell, groups.len()));
        }
        Ok(())
    };

    common::for_each_record(&mut opened, &input, |record| {
        let reference_name = record
            .reference_id()
            .and_then(|id| reference_names.get(id))
            .map(String::as_str);
        handle(
            &mut bundler.push(record, reference_name)?,
            &mut clusterer,
            &mut out,
        )
    })?;
    handle(&mut bundler.finish(), &mut clusterer, &mut out)?;

    if per_cell {
        let mut sorted_cells: Vec<u32> = (0..cells.len() as u32).collect();
        sorted_cells.sort_by(|&a, &b| cells[a as usize].cmp(&cells[b as usize]));
        let mut rank = vec![0u32; cells.len()];
        for (i, &cell) in sorted_cells.iter().enumerate() {
            rank[cell as usize] = i as u32;
        }
        // A cell seen twice keeps its last count, as Python's dict assignment does.
        let by_cell = |entries: &mut Vec<(u32, usize)>| {
            entries.sort_by_key(|&(cell, _)| rank[cell as usize]);
            entries.dedup_by(|later, kept| {
                let same = later.0 == kept.0;
                if same {
                    kept.1 = later.1;
                }
                same
            });
        };
        if matches.get_flag("wide_format_cell_counts") {
            let names: Vec<&str> = sorted_cells
                .iter()
                .map(|&c| cells[c as usize].as_str())
                .collect();
            writeln!(out, "gene\t{}", names.join("\t")).map_err(common::io_error)?;
            for (gene, entries) in &mut gene_counts {
                by_cell(entries);
                write!(out, "{gene}").map_err(common::io_error)?;
                let mut entries = entries.iter().peekable();
                for &cell in &sorted_cells {
                    let count = match entries.peek() {
                        Some(&&(c, n)) if c == cell => {
                            entries.next();
                            n
                        }
                        _ => 0,
                    };
                    write!(out, "\t{count}").map_err(common::io_error)?;
                }
                writeln!(out).map_err(common::io_error)?;
            }
        } else {
            writeln!(out, "gene\tcell\tcount").map_err(common::io_error)?;
            for (gene, entries) in &mut gene_counts {
                by_cell(entries);
                for &(cell, count) in entries.iter() {
                    writeln!(out, "{gene}\t{}\t{count}", cells[cell as usize])
                        .map_err(common::io_error)?;
                }
            }
        }
    }
    out.finish()?;

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
    fn counts_follow_pythons_rows_and_go_to_the_output_file() {
        use std::fmt::Write as _;
        let dir = tempfile::tempdir().unwrap();
        let sam = dir.path().join("in.sam");
        let mut text =
            String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:1000\n@SQ\tSN:c2\tLN:1000\n");
        // c1: G1 has UMIs AAAA and CCCC in cell X; c2: G1 has GGGG in X; G0 on c2.
        for (i, (contig, pos, umi, gene)) in [
            ("c1", 10, "AAAA", "G1"),
            ("c1", 20, "CCCC", "G1"),
            ("c2", 10, "GGGG", "G1"),
            ("c2", 30, "TTTT", "G0"),
        ]
        .iter()
        .enumerate()
        {
            writeln!(
                text,
                "r{i}\t0\t{contig}\t{pos}\t30\t4M\t*\t0\t0\tACGT\tIIII\tUB:Z:{umi}\tCB:Z:X\tGX:Z:{gene}"
            )
            .unwrap();
        }
        std::fs::write(&sam, text).unwrap();

        let run = |extra: &[&str]| -> String {
            let out = dir.path().join("counts.tsv");
            let mut args: Vec<String> = [
                "-I",
                sam.to_str().unwrap(),
                "--in-sam",
                "-S",
                out.to_str().unwrap(),
                "--extract-umi-method=tag",
                "--umi-tag=UB",
                "--cell-tag=CB",
                "--gene-tag=GX",
                "--log",
                dir.path().join("log").to_str().unwrap(),
            ]
            .iter()
            .map(|s| (*s).to_string())
            .collect();
            args.extend(extra.iter().map(|s| (*s).to_string()));
            assert_eq!(main(&args, &["count".to_string()]).unwrap(), 0);
            std::fs::read_to_string(out).unwrap()
        };
        assert_eq!(
            run(&["--per-cell"]),
            "gene\tcell\tcount\nG0\tX\t1\nG1\tX\t1\n"
        );
        assert_eq!(run(&[]), "gene\tcount\nG1\t2\nG0\t1\nG1\t1\n");
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
