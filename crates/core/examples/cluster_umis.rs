//! Group the UMIs observed at one alignment position, by all five methods.
//!
//! This is the crate's core question in miniature: given the UMIs seen at a
//! position and how often each occurred, which of them are sequencing errors
//! of which others? The five methods answer it with increasing
//! error-awareness, and the differences below are the reason `--method`
//! exists.
//!
//! The input is the worked example from the upstream repository's own
//! `step1_unit_test.py`, expressed against the 1.1.6 API. (That file calls
//! `clusterer(input_data.keys(), input_data, threshold=1)`, a three-argument
//! signature that 1.1.6's `UMIClusterer.__call__` no longer has — it is stale,
//! and was not used as an oracle for this port.)
//!
//! Run with:
//!
//! ```text
//! cargo run --example cluster_umis
//! ```

use std::collections::BTreeMap;

use umi_tools_core::{ClusterMethod, PortError, UmiClusterer};

fn main() -> Result<(), PortError> {
    // Six UMIs at one position. ACGT dominates; several others are one
    // substitution away from it and are almost certainly PCR/sequencing
    // errors of it rather than distinct molecules.
    let counts: BTreeMap<Vec<u8>, u64> = [
        (b"ACGT".to_vec(), 456),
        (b"AAAT".to_vec(), 90),
        (b"ACAT".to_vec(), 72),
        (b"TCGT".to_vec(), 2),
        (b"CCGT".to_vec(), 2),
        (b"ACAG".to_vec(), 1),
    ]
    .into_iter()
    .collect();

    println!("Observed at this position:");
    for (umi, count) in &counts {
        println!("  {:>6} x {count}", String::from_utf8_lossy(umi));
    }
    println!();

    for method in ClusterMethod::all() {
        let mut clusterer = UmiClusterer::new(method);
        let groups = clusterer.cluster(&counts, 1)?;

        println!("{:<12} {} group(s)", method.to_string(), groups.len());
        for group in &groups {
            let rendered: Vec<String> = group
                .iter()
                .map(|u| String::from_utf8_lossy(u).into_owned())
                .collect();
            // The first member is the inferred parent — the UMI the whole
            // group's reads are attributed to.
            if rendered.len() == 1 {
                println!("    {} (alone)", rendered[0]);
            } else {
                println!("    {} absorbs {}", rendered[0], rendered[1..].join(", "));
            }
        }
        println!();
    }

    println!("Reading the output:");
    println!("  unique      no error correction: every distinct UMI is a molecule.");
    println!("  percentile  drops UMIs below 1% of the median count.");
    println!("  cluster     one group per connected component of the ");
    println!("              'within one substitution' graph — greedy, and prone");
    println!("              to chaining unrelated UMIs together.");
    println!("  adjacency   covers each component with the fewest hub UMIs.");
    println!("  directional the default: connects a -> b only when");
    println!("              count(a) >= 2*count(b) - 1, encoding that an error");
    println!("              is seen at most about half as often as its parent.");

    Ok(())
}
