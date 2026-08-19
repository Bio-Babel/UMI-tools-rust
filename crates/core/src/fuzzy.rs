//! Approximate barcode matching (Python dependency: the `PyPI` `regex` module's
//! fuzzy-matching extensions).
//!
//! Three call sites in `umi_tools` compile a barcode into a fuzzy regular
//! expression and match another barcode against it:
//!
//! | Python | Pattern | Meaning |
//! |---|---|---|
//! | `whitelist_methods.checkError` | `(<barcode>){e<=n}` | up to `n` edits of any kind |
//! | `whitelist_methods.errorDetectAboveThreshold` | `(<barcode>){s<=n}` | up to `n` substitutions only |
//! | `network.CellClusterer._get_adj_list_directional` | `(<barcode>){e<=1}` and `(<reversed barcode>){i<=1}` | see [`cell_barcode_matches`] |
//!
//! The Rust `regex` crate has no approximate-matching extension and no
//! crates.io crate implements this syntax, so the three operations are
//! implemented directly. That is both smaller and easier to validate than a
//! general fuzzy-regex engine, and it makes the two semantics that are easy to
//! get wrong explicit:
//!
//! 1. **These are prefix matches, not whole-string matches.** Python uses
//!    `.match()`, which anchors at the start of the subject but not at the
//!    end. `(ACGT){e<=1}` matches `"ACGTXYZ"` with zero errors, because the
//!    pattern only has to consume a prefix. Verified against the live `regex`
//!    module rather than inferred from its documentation.
//! 2. **The three error classes are genuinely different.** `{e<=n}` allows
//!    substitutions, insertions and deletions; `{s<=n}` allows only
//!    substitutions, so the subject must be at least as long as the pattern
//!    and the comparison is a Hamming distance over the pattern's length;
//!    `{i<=n}` allows only insertions, so the pattern must appear as a
//!    subsequence of a prefix at most `n` bases longer than itself.
//!
//! # Input-form narrowing
//!
//! Python interpolates the barcode into a regular expression, so a barcode
//! containing regex metacharacters would be interpreted rather than matched
//! literally. Every barcode reaching these functions comes from a sequencer or
//! a whitelist file and consists of `ACGTN`, so this module treats the pattern
//! as a literal byte string. Recorded as an input-form narrowing in
//! `port_reports/umi_tools/05_design.md`.

/// Whether `pattern` matches a prefix of `text` within `max_errors`
/// substitutions, insertions and deletions
/// (Python: `regex.compile("(<pattern>){e<=n}").match(text)`).
///
/// The whole pattern must be consumed; the text need not be. The comparison is
/// therefore "can `pattern` be edited into some prefix of `text` in at most
/// `max_errors` steps".
///
/// # Examples
///
/// ```
/// use umi_tools_core::fuzzy::matches_within_errors;
/// assert!(matches_within_errors(b"ACGT", b"ACGA", 1));   // substitution
/// assert!(matches_within_errors(b"ACGT", b"CGT", 1));    // deletion
/// assert!(matches_within_errors(b"ACGT", b"TACGT", 1));  // insertion
/// assert!(matches_within_errors(b"ACGT", b"ACGTXYZ", 0)); // prefix, not whole
/// assert!(!matches_within_errors(b"ACGT", b"AAAA", 1));
/// ```
#[must_use]
pub fn matches_within_errors(pattern: &[u8], text: &[u8], max_errors: u32) -> bool {
    let m = pattern.len();
    let n = text.len();
    let cap = max_errors as usize;

    // row[j] is the cost of turning pattern[..i] into text[..j]. The answer is
    // whether any j leaves the final row within budget, which is what makes
    // this a prefix rather than a whole-string comparison.
    let mut prev: Vec<usize> = (0..=n).collect();
    let mut curr = vec![0usize; n + 1];

    for i in 1..=m {
        curr[0] = i;
        let mut row_min = curr[0];
        for j in 1..=n {
            let substitution = prev[j - 1] + usize::from(pattern[i - 1] != text[j - 1]);
            let deletion = prev[j] + 1;
            let insertion = curr[j - 1] + 1;
            curr[j] = substitution.min(deletion).min(insertion);
            row_min = row_min.min(curr[j]);
        }
        // Costs never decrease as `i` grows, so once a whole row exceeds the
        // budget no later row can come back under it.
        if row_min > cap {
            return false;
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    prev.iter().any(|&c| c <= cap)
}

/// Whether `pattern` matches a prefix of `text` within `max_subs`
/// substitutions and nothing else
/// (Python: `regex.compile("(<pattern>){s<=n}").match(text)`).
///
/// Because only substitutions are allowed, the matched region is exactly
/// `pattern.len()` bases long, so `text` shorter than `pattern` never matches
/// however large `max_subs` is.
///
/// # Examples
///
/// ```
/// use umi_tools_core::fuzzy::matches_within_substitutions;
/// assert!(matches_within_substitutions(b"ACGT", b"ACGA", 1));
/// assert!(matches_within_substitutions(b"ACG", b"ACGT", 0));   // prefix match
/// assert!(!matches_within_substitutions(b"ACGT", b"ACG", 9));  // too short
/// assert!(!matches_within_substitutions(b"ACGT", b"TACGT", 1));
/// ```
#[must_use]
pub fn matches_within_substitutions(pattern: &[u8], text: &[u8], max_subs: u32) -> bool {
    if text.len() < pattern.len() {
        return false;
    }
    let mismatches = pattern
        .iter()
        .zip(text.iter())
        .filter(|(p, t)| p != t)
        .count();
    mismatches <= max_subs as usize
}

/// Whether `pattern` matches a prefix of `text` using at most `max_inserts`
/// insertions and nothing else
/// (Python: `regex.compile("(<pattern>){i<=n}").match(text)`).
///
/// Equivalent to: some prefix of `text`, at most `max_inserts` bases longer
/// than `pattern`, contains `pattern` as a subsequence.
///
/// # Examples
///
/// ```
/// use umi_tools_core::fuzzy::matches_within_insertions;
/// assert!(matches_within_insertions(b"ACGT", b"TACGT", 1));
/// assert!(matches_within_insertions(b"ACGT", b"ACGTT", 0));  // prefix match
/// assert!(!matches_within_insertions(b"ACGT", b"ACGA", 1));  // needs a substitution
/// assert!(!matches_within_insertions(b"ACGT", b"ACG", 1));   // too short
/// ```
#[must_use]
pub fn matches_within_insertions(pattern: &[u8], text: &[u8], max_inserts: u32) -> bool {
    let limit = pattern.len().saturating_add(max_inserts as usize);
    if text.len() < pattern.len() {
        return false;
    }
    if pattern.is_empty() {
        // Zero characters of the pattern are already consumed, so an empty
        // pattern matches the empty prefix of any text.
        return true;
    }
    // Greedily consume the pattern; the greedy subsequence match uses the
    // fewest text characters possible, so if it finishes within `limit` no
    // other alignment does better.
    let mut pi = 0;
    for (ti, &b) in text.iter().enumerate() {
        if ti >= limit {
            break;
        }
        if pi < pattern.len() && pattern[pi] == b {
            pi += 1;
            if pi == pattern.len() {
                return ti < limit;
            }
        }
    }
    false
}

/// Whether cell barcode `b` should be treated as an error of cell barcode `a`
/// (Python source: the fuzzy branch of
/// `network.CellClusterer._get_adj_list_directional`).
///
/// Python builds two patterns from `a` and accepts a match from either:
///
/// ```text
/// comp_regex_err = regex.compile("(%s){e<=1}" % str(umi1))
/// comp_regex_del = regex.compile("(%s){i<=1}" % str(umi1)[::-1])
/// ...
/// if comp_regex_err.match(str(umi2)) or comp_regex_del.match(str(umi2)):
/// ```
///
/// The second pattern is built from the **reversed** barcode but matched
/// against the un-reversed subject, with a comment explaining it is meant to
/// catch deletions the first pattern misses (`UMI1 = "ATCG", UMI2 = "ATTC"`).
/// Reversing only one side means it almost never contributes a match — it
/// requires the reversed barcode to prefix the forward one, which happens only
/// for palindromes. It is reproduced exactly rather than repaired, because
/// "repairing" it would connect barcode pairs the Python original leaves
/// unconnected and so change which cells are merged.
///
/// # Examples
///
/// ```
/// use umi_tools_core::fuzzy::cell_barcode_matches;
/// assert!(cell_barcode_matches(b"ACGT", b"ACGA", 1));
/// assert!(!cell_barcode_matches(b"ACGT", b"ATTC", 1));
/// ```
#[must_use]
pub fn cell_barcode_matches(a: &[u8], b: &[u8], max_errors: u32) -> bool {
    if matches_within_errors(a, b, max_errors) {
        return true;
    }
    let reversed: Vec<u8> = a.iter().rev().copied().collect();
    matches_within_insertions(&reversed, b, max_errors)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Cases below were confirmed against the live `PyPI` `regex` module; see
    /// the module documentation for the exact patterns used.
    #[test]
    fn edit_matching_matches_python_regex() {
        let cases: &[(&[u8], &[u8], bool)] = &[
            (b"ACGT", b"ACGT", true),
            (b"ACGT", b"ACGA", true),
            (b"ACGT", b"ACG", true),
            (b"ACGT", b"ACGTT", true),
            (b"ACGT", b"ATTC", false),
            (b"ACG", b"ACGT", true),
            (b"ACGT", b"AAAA", false),
            (b"ACGT", b"TACGT", true),
            (b"ACGT", b"CGT", true),
            (b"ACGT", b"", false),
            (b"AAAA", b"AAAAAA", true),
            (b"ACGT", b"ACGTXYZ", true),
        ];
        for (pattern, text, expected) in cases {
            assert_eq!(
                matches_within_errors(pattern, text, 1),
                *expected,
                "(e<=1) pattern={} text={}",
                String::from_utf8_lossy(pattern),
                String::from_utf8_lossy(text)
            );
        }
    }

    #[test]
    fn substitution_matching_matches_python_regex() {
        let cases: &[(&[u8], &[u8], bool)] = &[
            (b"ACGT", b"ACGT", true),
            (b"ACGT", b"ACGA", true),
            (b"ACGT", b"ACG", false),
            (b"ACGT", b"ACGTT", true),
            (b"ACGT", b"ATTC", false),
            (b"ACG", b"ACGT", true),
            (b"ACGT", b"AAAA", false),
            (b"ACGT", b"TACGT", false),
            (b"ACGT", b"CGT", false),
            (b"ACGT", b"", false),
            (b"AAAA", b"AAAAAA", true),
            (b"ACGT", b"ACGTXYZ", true),
        ];
        for (pattern, text, expected) in cases {
            assert_eq!(
                matches_within_substitutions(pattern, text, 1),
                *expected,
                "(s<=1) pattern={} text={}",
                String::from_utf8_lossy(pattern),
                String::from_utf8_lossy(text)
            );
        }
    }

    #[test]
    fn insertion_matching_matches_python_regex() {
        let cases: &[(&[u8], &[u8], bool)] = &[
            (b"ACGT", b"ACGT", true),
            (b"ACGT", b"ACGA", false),
            (b"ACGT", b"ACG", false),
            (b"ACGT", b"ACGTT", true),
            (b"ACGT", b"ATTC", false),
            (b"ACG", b"ACGT", true),
            (b"ACGT", b"AAAA", false),
            (b"ACGT", b"TACGT", true),
            (b"ACGT", b"CGT", false),
            (b"ACGT", b"", false),
            (b"AAAA", b"AAAAAA", true),
            (b"ACGT", b"ACGTXYZ", true),
        ];
        for (pattern, text, expected) in cases {
            assert_eq!(
                matches_within_insertions(pattern, text, 1),
                *expected,
                "(i<=1) pattern={} text={}",
                String::from_utf8_lossy(pattern),
                String::from_utf8_lossy(text)
            );
        }
    }

    #[test]
    fn zero_error_budget_is_a_plain_prefix_match() {
        assert!(matches_within_errors(b"ACGT", b"ACGTXX", 0));
        assert!(!matches_within_errors(b"ACGT", b"ACGA", 0));
        assert!(matches_within_substitutions(b"ACGT", b"ACGTXX", 0));
        assert!(!matches_within_substitutions(b"ACGT", b"ACGA", 0));
    }

    #[test]
    fn empty_pattern_matches_anything() {
        assert!(matches_within_errors(b"", b"ACGT", 0));
        assert!(matches_within_errors(b"", b"", 0));
        assert!(matches_within_substitutions(b"", b"", 0));
        assert!(matches_within_insertions(b"", b"ACGT", 0));
    }

    #[test]
    fn larger_budgets_admit_more_edits() {
        assert!(!matches_within_errors(b"ACGT", b"AAAA", 1));
        assert!(!matches_within_errors(b"ACGT", b"AAAA", 2));
        assert!(matches_within_errors(b"ACGT", b"AAAA", 3));
        assert!(matches_within_substitutions(b"ACGT", b"AAAA", 3));
        assert!(!matches_within_substitutions(b"ACGT", b"AAAA", 2));
    }

    #[test]
    fn cell_matching_accepts_either_pattern() {
        assert!(cell_barcode_matches(b"ACGT", b"ACGA", 1));
        assert!(cell_barcode_matches(b"ACGT", b"ACG", 1));
        assert!(cell_barcode_matches(b"ACGT", b"ACGTT", 1));
        assert!(!cell_barcode_matches(b"ACGT", b"ATTC", 1));
        // The reversed-pattern branch only fires for a palindrome.
        assert!(cell_barcode_matches(b"AAAA", b"AAAAAA", 1));
    }

    #[test]
    fn non_ascii_bytes_are_compared_by_value() {
        assert!(matches_within_errors(&[0x00, 0xff], &[0x00, 0xff], 0));
        assert!(matches_within_errors(&[0x00, 0xff], &[0x00, 0xfe], 1));
        assert!(!matches_within_errors(&[0x00, 0xff], &[0x01, 0xfe], 1));
    }
}
