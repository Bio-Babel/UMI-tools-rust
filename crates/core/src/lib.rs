//! Pure-Rust core of the UMI-tools port.
//!
//! [UMI-tools](https://github.com/CGATOxford/UMI-tools) identifies PCR
//! duplicates in sequencing data by grouping reads that share a mapping
//! position and a Unique Molecular Identifier (UMI), correcting for sequencing
//! errors in the UMI itself. This crate is a port of the Python package at
//! version 1.1.6, module for module, with the Python source recorded in each
//! module's own documentation.
//!
//! # Layout
//!
//! | Rust module | Python source |
//! |---|---|
//! | [`distance`] | `umi_tools/_dedup_umi.pyx` |
//! | [`network`] | `umi_tools/network.py` |
//! | [`umi_methods`] | `umi_tools/umi_methods.py` |
//! | [`extract`] | `umi_tools/extract_methods.py` |
//! | [`whitelist`] | `umi_tools/whitelist_methods.py` |
//! | [`sam`], [`bundles`] | `umi_tools/sam_methods.py` |
//! | [`utilities`], [`alignment`], [`options`] | `umi_tools/Utilities.py` |
//! | [`numeric`] | the NumPy/SciPy subset the package uses |
//! | [`stats`] | `umi_tools/dedup.py`'s pandas aggregation |
//! | [`bktree`] | the `pybktree` dependency |
//! | [`fuzzy`] | the `regex` module's fuzzy-matching subset |
//!
//! This crate has **no `PyO3` dependency**. Python bindings live in
//! `crates/pyo3-bindings`, which depends on this crate and never the other way
//! round, so dropping them later is a directory deletion with no changes here.
//!
//! # Examples
//!
//! Group a position's UMIs with the directional method, the default UMI-tools
//! uses:
//!
//! ```
//! use std::collections::BTreeMap;
//! use umi_tools_core::{ClusterMethod, UmiClusterer};
//!
//! let mut counts: BTreeMap<Vec<u8>, u64> = BTreeMap::new();
//! counts.insert(b"ACGT".to_vec(), 10);
//! counts.insert(b"ACGA".to_vec(), 2);   // one error away from ACGT
//! counts.insert(b"TTTT".to_vec(), 5);
//!
//! let mut clusterer = UmiClusterer::new(ClusterMethod::Directional);
//! let groups = clusterer.cluster(&counts, 1)?;
//!
//! assert_eq!(groups.len(), 2);
//! assert_eq!(groups[0][0], b"ACGT".to_vec());  // the parent UMI leads its group
//! # Ok::<(), umi_tools_core::PortError>(())
//! ```
#![warn(missing_docs)]
#![warn(clippy::pedantic)]
// The port runs under `clippy::pedantic`. The allowances below are deliberate,
// and each is here because the lint conflicts with reproducing the Python
// original rather than because the code is careless.
//
// `float_cmp`: exact float comparison is the point in several places — a knee
// threshold Python computed by the same expression, a `linspace` endpoint
// Python assigns rather than accumulates, a covariance the differential test
// asserts is bit-identical to SciPy's. Rounding those comparisons would hide
// exactly the divergences the tests exist to catch.
#![allow(clippy::float_cmp)]
// The numeric-cast family: every one mirrors an arithmetic step in the Python
// source, where ints are unbounded and int/float conversion is implicit. The
// values are read counts and barcode lengths — nowhere near the magnitudes at
// which the precision or sign warnings could bite — and restructuring the
// arithmetic to avoid the cast would change the result.
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_possible_wrap)]
#![allow(clippy::cast_sign_loss)]
// `struct_excessive_bools`: the option structs mirror `optparse` flag sets one
// field per flag. Grouping them into sub-structs would obscure the
// correspondence with the command line, which is what makes them checkable
// against the Python source.
#![allow(clippy::struct_excessive_bools)]
// `too_many_lines`: a few functions — the subcommand `main`s, the bundler's
// per-read dispatch — follow a long Python function step for step. Splitting
// them would make that correspondence harder to verify, which is the property
// that matters most here.
#![allow(clippy::too_many_lines)]
#![allow(clippy::must_use_candidate)]
#![allow(clippy::missing_errors_doc)]
#![allow(clippy::missing_panics_doc)]
#![allow(clippy::module_name_repetitions)]
#![allow(clippy::similar_names)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::needless_pass_by_value)]

pub mod alignment;
pub mod bktree;
pub mod bundles;
pub mod distance;
pub mod error;
pub mod extract;
pub mod fuzzy;
pub mod mersenne;
pub mod network;
pub mod numeric;
pub mod options;
pub mod pattern;
pub mod sam;
pub mod stats;
pub mod umi_methods;
pub mod utilities;
pub mod whitelist;

pub use bundles::{Bundle, BundleEntry, BundleItem, BundleReads, GetBundles};
pub use distance::{edit_distance, edit_distance_f64, within_threshold};
pub use error::{PortError, Result};
pub use extract::{ExtractConfig, ExtractFilterAndUpdate};
pub use network::{CellClusterer, ClusterMethod, ReadDeduplicator, UmiClusterer};
pub use umi_methods::{FastqIterator, JoinedFastqIterator, Record};
pub use whitelist::{knee_estimate_density, knee_estimate_distance, KneeResult};

/// The Python package version this crate was ported from.
///
/// The crate's own version tracks it: `umi_tools-core` 1.1.6 corresponds to
/// `umi_tools` 1.1.6. See the README's version-mapping note.
pub const VERSION: &str = "1.1.6";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_matches_python_package() {
        assert_eq!(VERSION, "1.1.6");
    }
}
