//! Hamming distance between barcodes (Python source: `umi_tools/_dedup_umi.pyx`).
//!
//! The Python original is a 19-line Cython function compiled into
//! `umi_tools._dedup_umi`:
//!
//! ```text
//! cpdef int edit_distance(a, b):
//!     cdef char * aa = a
//!     cdef char * bb = b
//!     ...
//!     if la != lb:
//!         return np.inf
//!     for k from 0 <= k < la:
//!         if aa[k] != bb[k]: c += 1
//!     return c
//! ```
//!
//! Despite the name it is a Hamming distance, not a Levenshtein distance:
//! there is no insertion or deletion handling, and unequal lengths are
//! rejected outright rather than aligned.
//!
//! # Deviation from the shipped Python build
//!
//! The `.pyx` above returns `np.inf` for unequal lengths, but the
//! **checked-in generated `_dedup_umi.c`** — which is what `pyproject.toml`
//! actually compiles — emits `np.Inf`, removed in `NumPy` 2.0. So the shipped
//! extension raises `AttributeError` on unequal-length input rather than
//! returning infinity. This module implements the `.pyx` intent, because that
//! is the gold standard; the `AttributeError` is an artefact of a stale
//! generated file. See `port_reports/umi_tools/05_design.md` deviation D1.

/// Hamming distance between two equal-length byte strings.
///
/// Returns `None` when the inputs differ in length, mirroring the Python
/// original's infinite-distance sentinel: an infinite distance is never
/// `<= threshold`, so every call site treats it as "not neighbours".
///
/// Barcodes are compared as raw bytes, not as text. The Python original
/// coerces its arguments to `char *`, so it too compares bytes and rejects
/// `str` input with `TypeError: expected bytes, str found`.
///
/// # Examples
///
/// ```
/// use umi_tools_core::distance::edit_distance;
/// assert_eq!(edit_distance(b"ACGT", b"ACGA"), Some(1));
/// assert_eq!(edit_distance(b"ACGT", b"ACG"), None);
/// ```
#[must_use]
pub fn edit_distance(a: &[u8], b: &[u8]) -> Option<u32> {
    if a.len() != b.len() {
        return None;
    }
    Some(
        a.iter()
            .zip(b.iter())
            .filter(|(x, y)| x != y)
            .count()
            .try_into()
            .unwrap_or(u32::MAX),
    )
}

/// Hamming distance as an `f64`, with unequal lengths yielding
/// [`f64::INFINITY`].
///
/// This is the form the Python original's *value* takes, and it is what
/// [`crate::umi_methods::get_average_umi_distance`] needs, since that function
/// averages the raw distances rather than comparing them to a threshold.
#[must_use]
pub fn edit_distance_f64(a: &[u8], b: &[u8]) -> f64 {
    match edit_distance(a, b) {
        Some(d) => f64::from(d),
        None => f64::INFINITY,
    }
}

/// Whether two barcodes are within `threshold` substitutions of each other.
///
/// Every call site in the Python original spells this as
/// `edit_distance(a, b) <= threshold`. Capturing it once means no caller has
/// to re-derive that an infinite distance is never within threshold.
///
/// # Examples
///
/// ```
/// use umi_tools_core::distance::within_threshold;
/// assert!(within_threshold(b"ACGT", b"ACGA", 1));
/// assert!(!within_threshold(b"ACGT", b"ACG", 99));
/// ```
#[must_use]
pub fn within_threshold(a: &[u8], b: &[u8], threshold: u32) -> bool {
    matches!(edit_distance(a, b), Some(d) if d <= threshold)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_is_zero() {
        assert_eq!(edit_distance(b"ACGT", b"ACGT"), Some(0));
    }

    #[test]
    fn single_substitution() {
        assert_eq!(edit_distance(b"ACGT", b"ACGA"), Some(1));
    }

    #[test]
    fn all_positions_differ() {
        assert_eq!(edit_distance(b"ACGT", b"TGCA"), Some(4));
    }

    #[test]
    fn empty_inputs_are_distance_zero() {
        assert_eq!(edit_distance(b"", b""), Some(0));
    }

    #[test]
    fn unequal_lengths_are_none() {
        assert_eq!(edit_distance(b"ACGT", b"ACG"), None);
        assert_eq!(edit_distance(b"", b"A"), None);
        assert_eq!(edit_distance(b"AAAA", b"AAAAA"), None);
    }

    #[test]
    fn non_ascii_bytes_compare_by_value() {
        assert_eq!(edit_distance(&[0x00, 0x01], &[0x00, 0x02]), Some(1));
        assert_eq!(edit_distance(&[0xff], &[0xff]), Some(0));
    }

    #[test]
    fn f64_form_uses_infinity_for_unequal_lengths() {
        assert_eq!(edit_distance_f64(b"ACGT", b"ACGA"), 1.0);
        assert!(edit_distance_f64(b"ACGT", b"ACG").is_infinite());
    }

    #[test]
    fn within_threshold_rejects_unequal_lengths_at_any_threshold() {
        assert!(within_threshold(b"ACGT", b"ACGT", 0));
        assert!(within_threshold(b"ACGT", b"ACGA", 1));
        assert!(!within_threshold(b"ACGT", b"ACAA", 1));
        assert!(!within_threshold(b"ACGT", b"ACG", u32::MAX));
    }

    /// Reference values dumped from the live Python oracle by
    /// `port_reports/umi_tools/data_exports/export_distance_numeric.py`.
    #[test]
    fn matches_python_reference_values() {
        let raw = include_str!("../tests/fixtures/reference_distance_numeric.json");
        let doc: serde_json::Value = serde_json::from_str(raw).expect("fixture parses");
        let cases = doc["edit_distance"]
            .as_array()
            .expect("edit_distance array");
        assert!(!cases.is_empty(), "fixture must carry cases");

        for case in cases {
            let a: Vec<u8> = case["a"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| u8::try_from(v.as_u64().unwrap()).unwrap())
                .collect();
            let b: Vec<u8> = case["b"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| u8::try_from(v.as_u64().unwrap()).unwrap())
                .collect();
            let expected = case["distance"].as_u64().map(|d| d as u32);
            assert_eq!(
                edit_distance(&a, &b),
                expected,
                "edit_distance({a:?}, {b:?})"
            );
        }
    }
}
