//! `umi_tools count_tab` — count molecules per gene from a flat file
//! (Python source: `umi_tools/count_tab.py`).
//!
//! The same counting as [`crate::count`], but reading a two-column
//! `read_id<TAB>gene` table instead of an alignment — the shape
//! `featureCounts -R CORE` produces. The input must be sorted by gene, because
//! the count is a streaming group-by.

use std::io::{BufWriter, Write};

use umi_tools_core::network::UmiClusterer;
use umi_tools_core::sam::{GeneCountTab, StringBarcodeGetter};
use umi_tools_core::utilities;
use umi_tools_core::Result;

use crate::common;
use crate::documentation;
use crate::options::{self, OptionGroups};

/// The subcommand's short description.
pub const ABOUT: &str = "count_tab - Count reads per gene from flatfile using UMIs";

/// Runs `umi_tools count_tab` (Python source: `count_tab.main`).
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Value`] for an invalid option, [`umi_tools_core::PortError::Assertion`]
/// for a malformed input line, and [`umi_tools_core::PortError::Io`] for input problems.
pub fn main(args: &[String], argv: &[String]) -> Result<i32> {
    if crate::wants_help(args) {
        // `optparse`'s layout, reproduced verbatim; see `crate::help`.
        print!("{}", crate::help::COUNT_TAB);
        return Ok(0);
    }

    let started = std::time::Instant::now();
    let groups = OptionGroups {
        group_dedup: false,
        sam: false,
        sbcram: true,
        exclude_output_format: true,
        ..OptionGroups::default()
    };

    let command = options::build_parser("count_tab", ABOUT, groups)
        .after_help(format!(
            "{}{}",
            documentation::GENERIC_DOCSTRING_GDC,
            documentation::GENERIC_DOCSTRING_SBCRAM_INPUT
        ))
        .arg(
            clap::Arg::new("bc_sep")
                .overrides_with("bc_sep")
                .long("barcode-separator")
                .num_args(1)
                .help("separator between read id and UMI and (optionally) the cell barcode"),
        )
        .arg(
            clap::Arg::new("per_cell")
                .overrides_with("per_cell")
                .long("per-cell")
                .action(clap::ArgAction::SetTrue)
                .help("read name includes cell barcode as well as UMI: read[sep]UMI[sep]CB"),
        );

    let matches = match command
        .try_get_matches_from(std::iter::once("count_tab".to_string()).chain(args.iter().cloned()))
    {
        Ok(m) => m,
        Err(e) => {
            e.print().ok();
            return Ok(i32::from(e.use_stderr()));
        }
    };

    let io = options::io_options(&matches)?;
    let run_id = options::start(&io, argv)?;
    let dedup_options = options::dedup_options(&matches)?;

    let per_cell = matches.get_flag("per_cell");
    let sep = matches
        .get_one::<String>("bc_sep")
        .cloned()
        .unwrap_or_else(|| "_".to_string());
    let getter = if per_cell {
        StringBarcodeGetter::CellAndUmi
    } else {
        StringBarcodeGetter::UmiOnly
    };

    let input: Box<dyn std::io::BufRead> = match io.stdin.as_deref() {
        Some(path) => Box::new(utilities::open_file(path)?),
        None => Box::new(std::io::BufReader::new(std::io::stdin())),
    };

    let mut out = BufWriter::new(std::io::stdout());
    if per_cell {
        writeln!(out, "cell\tgene\tcount").map_err(common::io_error)?;
    } else {
        writeln!(out, "gene\tcount").map_err(common::io_error)?;
    }

    let mut clusterer = UmiClusterer::new(dedup_options.method);
    let mut n_input = 0u64;
    let mut n_output = 0u64;

    for group in GeneCountTab::new(input, getter, &sep) {
        let (gene, per_cell_counts) = group?;
        // Cells in first-seen order, and UMIs likewise, because `count_tab`
        // writes rows as it goes and the clustering order decides the count.
        for entry in per_cell_counts.cells() {
            n_input += entry.counts.values().sum::<u64>();
            let groups =
                clusterer.cluster_ordered(&entry.umis, &entry.counts, dedup_options.threshold)?;
            let gene_count = groups.len();

            if per_cell {
                let cell = entry
                    .cell
                    .as_ref()
                    .map(|c| utilities::decode_utf8(c))
                    .transpose()?
                    .unwrap_or_default();
                writeln!(out, "{cell}\t{gene}\t{gene_count}").map_err(common::io_error)?;
            } else {
                writeln!(out, "{gene}\t{gene_count}").map_err(common::io_error)?;
                // Python only accumulates nOutput on the non-per-cell path.
                n_output += gene_count as u64;
            }
        }
    }
    out.flush().map_err(common::io_error)?;

    let _ = n_input;
    utilities::info(&format!("Number of reads counted: {n_output}"));
    options::stop(&io, &run_id, started.elapsed().as_secs())?;
    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use umi_tools_core::PortError;

    #[test]
    fn about_names_the_command() {
        assert!(ABOUT.starts_with("count_tab - "));
    }

    #[test]
    fn counts_molecules_per_gene() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.tsv");
        // geneA: AAA seen twice and AAT once (one error away), CCC distinct.
        std::fs::write(
            &input,
            "r1_AAA\tgeneA\nr2_AAA\tgeneA\nr3_AAT\tgeneA\nr4_CCC\tgeneA\nr5_GGG\tgeneB\n",
        )
        .unwrap();

        let args: Vec<String> = vec!["-I".to_string(), input.to_string_lossy().into_owned()];
        assert_eq!(main(&args, &["count_tab".to_string()]).unwrap(), 0);
    }

    #[test]
    fn a_malformed_line_is_an_assertion_error() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.tsv");
        std::fs::write(&input, "r1_AAA\tgeneA\textra\n").unwrap();
        let args: Vec<String> = vec!["-I".to_string(), input.to_string_lossy().into_owned()];
        let err = main(&args, &["count_tab".to_string()]).unwrap_err();
        assert!(matches!(err, PortError::Assertion(_)));
    }
}
