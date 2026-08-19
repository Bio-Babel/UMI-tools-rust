//! Deduplicate a real BAM, the way `umi_tools dedup` does.
//!
//! Step 5 of the package's quick-start guide, as a library program: read an
//! alignment, group its reads by position and UMI, merge UMIs that are
//! sequencing errors of one another, and report one read per original
//! molecule.
//!
//! This is the whole pipeline the crate exists for, and it shows how the three
//! pieces fit: [`GetBundles`] turns a stream of records into per-position
//! bundles, [`ReadDeduplicator`] decides which UMIs within a bundle came from
//! the same molecule, and the caller writes out the survivors.
//!
//! Run with:
//!
//! ```text
//! cargo run --example deduplicate_bam -- path/to/alignment.bam
//! ```
//!
//! With no argument it uses the small BAM checked in under
//! `crates/core/tests/fixtures/bam/`.

use std::path::{Path, PathBuf};

use umi_tools_core::alignment::open_input_alignments;
use umi_tools_core::bundles::{insertion_order, BundleFlags, BundleItem, GetBundles};
use umi_tools_core::network::{ClusterMethod, CountMap, ReadDeduplicator};
use umi_tools_core::options::BundleOptions;
use umi_tools_core::sam::BarcodeGetter;
use umi_tools_core::utilities::AlignmentFormat;
use umi_tools_core::PortError;

fn default_input() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/bam/chr19_1mb.bam")
}

fn main() -> Result<(), PortError> {
    let input = std::env::args()
        .nth(1)
        .map_or_else(default_input, PathBuf::from);
    println!("Reading {}", input.display());

    let mut opened = open_input_alignments(&input, AlignmentFormat::Bam, None)?;
    let records = opened.read_all()?;
    println!("{} alignment records\n", records.len());

    // The UMI is appended to the read name by `umi_tools extract`, separated
    // by an underscore — the layout `extract_barcodes.rs` produces.
    let getter = BarcodeGetter::ReadId {
        cell_barcode: false,
        sep: "_".to_string(),
    };

    // `--random-seed` makes the tie-break between equally good reads
    // reproducible; without it the choice among tied reads is arbitrary.
    let options = BundleOptions {
        random_seed: Some(123_456_789),
        ..BundleOptions::default()
    };

    let mut bundler = GetBundles::new(&options, BundleFlags::default(), getter);
    let mut deduplicator = ReadDeduplicator::new(ClusterMethod::Directional);

    let mut reads_in = 0u64;
    let mut reads_out = 0u64;
    let mut examples_shown = 0;

    // Bundles are emitted as the scan passes them, so memory stays bounded
    // however large the input is.
    let handle = |items: Vec<BundleItem>,
                  deduplicator: &mut ReadDeduplicator,
                  reads_in: &mut u64,
                  reads_out: &mut u64,
                  examples_shown: &mut usize|
     -> Result<(), PortError> {
        for item in items {
            let BundleItem::Bundle { bundle, .. } = item else {
                continue;
            };
            let counts: CountMap = bundle
                .iter()
                .map(|(umi, entry)| (umi.clone(), entry.count))
                .collect();
            *reads_in += counts.values().sum::<u64>();

            // The UMIs in the order they were first seen; the clustering order
            // decides which UMI leads each group.
            let order = insertion_order(&bundle);
            let result = deduplicator.dedup_ordered(&order, &counts, 1)?;
            *reads_out += result.final_umis.len() as u64;

            // Show the first few positions where deduplication actually merged
            // something — those are the interesting ones.
            if *examples_shown < 3 && result.final_umis.len() < counts.len() {
                *examples_shown += 1;
                println!(
                    "position with {} distinct UMIs -> {} molecules",
                    counts.len(),
                    result.final_umis.len()
                );
                for (umi, count) in result.final_umis.iter().zip(&result.umi_counts) {
                    println!(
                        "    {} accounts for {count} read(s)",
                        String::from_utf8_lossy(umi)
                    );
                }
                println!();
            }
        }
        Ok(())
    };

    for record in umi_tools_core::alignment::fetch(&records, false) {
        let reference_name = opened.reference_name(record);
        let items = bundler.push(record.clone(), reference_name.as_deref())?;
        handle(
            items,
            &mut deduplicator,
            &mut reads_in,
            &mut reads_out,
            &mut examples_shown,
        )?;
    }
    handle(
        bundler.finish(),
        &mut deduplicator,
        &mut reads_in,
        &mut reads_out,
        &mut examples_shown,
    )?;

    let clusterer = deduplicator.clusterer();
    println!("Summary");
    println!("  reads in            {reads_in}");
    println!("  molecules out       {reads_out}");
    if reads_in > 0 {
        let rate = 100.0 * (1.0 - reads_out as f64 / reads_in as f64);
        println!("  duplication rate    {rate:.1}%");
    }
    println!("  positions           {}", clusterer.positions());
    if clusterer.positions() > 0 {
        println!(
            "  mean UMIs/position  {:.2}",
            clusterer.total_umis_per_position() as f64 / clusterer.positions() as f64
        );
        println!(
            "  max UMIs/position   {}",
            clusterer.max_umis_per_position()
        );
    }

    println!("\nRead events");
    for (event, count) in bundler.read_events() {
        println!("  {event}: {count}");
    }

    Ok(())
}
