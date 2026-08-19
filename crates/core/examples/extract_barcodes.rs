//! Move UMIs and cell barcodes out of read sequences and into read names.
//!
//! This is step 3 of the package's quick-start guide and step 2 of its
//! single-cell tutorial: before alignment, the barcode has to leave the
//! sequence (or the aligner will try to align it) and travel somewhere the
//! downstream tools can still find it — the read name.
//!
//! Both ways of saying where the barcode is are shown: a positional mask, and
//! a regular expression for the variable-length case.
//!
//! Run with:
//!
//! ```text
//! cargo run --example extract_barcodes
//! ```

use umi_tools_core::extract::{ExtractConfig, ExtractFilterAndUpdate, ExtractMethod};
use umi_tools_core::umi_methods::QualityEncoding;
use umi_tools_core::{PortError, Record};

fn show(label: &str, before: &Record, after: Option<&Record>) {
    println!("{label}");
    println!("  in   {}", before.seq);
    match after {
        Some(read) => {
            println!("  out  {}", read.seq);
            println!("  name {}", read.identifier);
        }
        None => println!("  out  (read discarded)"),
    }
    println!();
}

fn main() -> Result<(), PortError> {
    // --- The positional mask -------------------------------------------
    //
    // `NNNXXXXNN`: the first three bases are UMI, the next four stay on the
    // read, the last two are more UMI. Everything past the mask is untouched.
    let mut string_extractor = ExtractFilterAndUpdate::new(ExtractConfig {
        pattern: Some("NNNXXXXNN".to_string()),
        ..ExtractConfig::new()
    })?;

    let original = Record::new("read1 some description", "ACGTTTTGGCCCCC", "IIIIIIIIIIIIII");
    let mut read = original.clone();
    let kept = string_extractor.process(&mut read, None)?;
    show("string mask NNNXXXXNN", &original, kept.then_some(&read));

    // --- A cell barcode as well ----------------------------------------
    //
    // `C` marks cell-barcode bases. The barcodes join the read name as
    // `<id>_<cell>_<umi>`, which is what the alignment-side parsers expect.
    let mut cell_extractor = ExtractFilterAndUpdate::new(ExtractConfig {
        pattern: Some("CCCCCCNNNNNNNNNN".to_string()),
        extract_cell: true,
        ..ExtractConfig::new()
    })?;

    let original = Record::new(
        "read2",
        "GATCGAACGTACGTACGTAAAACCCC",
        "IIIIIIIIIIIIIIIIIIIIIIIIII",
    );
    let mut read = original.clone();
    let kept = cell_extractor.process(&mut read, None)?;
    show(
        "string mask CCCCCCNNNNNNNNNN",
        &original,
        kept.then_some(&read),
    );

    // --- A regular expression ------------------------------------------
    //
    // Needed when the barcode is not at a fixed offset. Group *names* carry
    // the meaning: `cell_*` is cell barcode, `umi_*` is UMI, `discard_*` is
    // dropped entirely. Here a fixed linker sits between two halves of the
    // cell barcode, which is the inDrop protocol's layout.
    let mut regex_extractor = ExtractFilterAndUpdate::new(ExtractConfig {
        method: ExtractMethod::Regex,
        pattern: Some(
            "(?P<cell_1>.{4})(?P<discard_2>GGGG)(?P<cell_3>.{4})(?P<umi_1>.{3}).*".to_string(),
        ),
        extract_cell: true,
        ..ExtractConfig::new()
    })?;

    let original = Record::new("read3", "AAAAGGGGCCCCTTTACGT", "IIIIIIIIIIIIIIIIIII");
    let mut read = original.clone();
    let kept = regex_extractor.process(&mut read, None)?;
    show(
        "regex with a linker to discard",
        &original,
        kept.then_some(&read),
    );

    // --- Quality filtering ---------------------------------------------
    //
    // A UMI base called with low confidence makes the whole UMI unreliable, so
    // a read can be dropped outright or the offending bases masked to `N`.
    let mut masking_extractor = ExtractFilterAndUpdate::new(ExtractConfig {
        pattern: Some("NNN".to_string()),
        quality_encoding: Some(QualityEncoding::Phred33),
        quality_filter_mask: Some(30),
        ..ExtractConfig::new()
    })?;

    let original = Record::new("read4", "ACGTTTT", "I!IIIII");
    let mut read = original.clone();
    let kept = masking_extractor.process(&mut read, None)?;
    show(
        "quality masking below Q30",
        &original,
        kept.then_some(&read),
    );

    println!("Event tallies from the masking extractor:");
    for (event, count) in masking_extractor.read_counts() {
        println!("  {event}: {count}");
    }

    Ok(())
}
