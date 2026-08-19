//! Pick the real cell barcodes out of an observed count distribution.
//!
//! Step 1 of the package's single-cell tutorial. A droplet run produces far
//! more distinct barcodes than there were cells — most are sequencing errors
//! or ambient RNA — so the count distribution has a characteristic "knee"
//! separating real cells from background, and this finds it.
//!
//! Both estimators are shown on the same data, along with the error-correction
//! map that folds near-miss barcodes into their whitelisted neighbours.
//!
//! Run with:
//!
//! ```text
//! cargo run --example cell_whitelist
//! ```

use umi_tools_core::whitelist::{
    get_cell_whitelist, knee_estimate_density, knee_estimate_distance, CellCounts, KneeMethod,
};
use umi_tools_core::PortError;

/// A synthetic distribution with the shape a real run has: a handful of true
/// cells at high counts, then a long tail of background.
fn observed_counts() -> CellCounts {
    let mut counts = CellCounts::new();
    for (barcode, count) in [
        ("AAAACCCC", 1000u64),
        ("AAAACCCA", 800), // one substitution from the first
        ("AAAACCTT", 600),
        ("AAAATTTT", 400),
        ("TTTTTTTT", 200),
        // The background tail.
        ("CCCCCCCC", 10),
        ("CCCCCCCG", 8),
        ("CCCCCCGG", 6),
        ("CCCCCGGG", 4),
        ("GGGGGGGG", 2),
    ] {
        counts.set(barcode.as_bytes(), count);
    }
    counts
}

fn render(barcodes: Option<&std::collections::BTreeSet<Vec<u8>>>) -> String {
    match barcodes {
        None => "(no threshold accepted)".to_string(),
        Some(set) => set
            .iter()
            .map(|b| String::from_utf8_lossy(b).into_owned())
            .collect::<Vec<_>>()
            .join(", "),
    }
}

fn main() -> Result<(), PortError> {
    let counts = observed_counts();

    println!("Observed barcodes, most abundant first:");
    for (barcode, count) in umi_tools_core::whitelist::most_common(&counts) {
        println!("  {:>10} x {count}", String::from_utf8_lossy(&barcode));
    }
    println!();

    // --- The distance method (the default) ------------------------------
    //
    // Draws a line from the first to the last point of the cumulative-count
    // curve and takes the point furthest from it — the elbow.
    let distance = knee_estimate_distance(&counts, None)?;
    println!("distance knee");
    println!("  cut at rank        {:?}", distance.knee_index);
    println!(
        "  selected           {}",
        render(distance.barcodes.as_ref())
    );
    println!(
        "  cumulative curve   {:?}...",
        &distance.cumulative[..distance.cumulative.len().min(5)]
    );
    println!();

    // --- The density method ---------------------------------------------
    //
    // Fits a kernel density estimate to the log counts and cuts at a local
    // minimum — the trough between the "real cells" and "background" modes.
    let density = knee_estimate_density(&counts, None, None)?;
    println!("density knee");
    println!("  selected           {}", render(density.barcodes.as_ref()));
    println!(
        "  candidate cuts     {:?} (barcodes each would accept)",
        density.local_mins_counts
    );
    println!("  chosen candidate   {:?}", density.selected_local_min);
    println!();

    // --- Error correction ------------------------------------------------
    //
    // Barcodes just off the whitelist are usually sequencing errors of a
    // whitelisted one. A barcode within threshold of *more than one* is
    // ambiguous and is dropped rather than guessed at.
    let (result, corrections) = get_cell_whitelist(&counts, KneeMethod::Distance, None, None, 1)?;
    println!("with error correction at edit distance 1");
    println!("  whitelist          {}", render(result.barcodes.as_ref()));
    match corrections {
        None => println!("  corrections        (none requested)"),
        Some(map) if map.is_empty() => println!("  corrections        none found"),
        Some(map) => {
            for (parent, children) in &map {
                let rendered: Vec<String> = children
                    .iter()
                    .map(|c| String::from_utf8_lossy(c).into_owned())
                    .collect();
                println!(
                    "  {} <- {}",
                    String::from_utf8_lossy(parent),
                    rendered.join(", ")
                );
            }
        }
    }

    Ok(())
}
