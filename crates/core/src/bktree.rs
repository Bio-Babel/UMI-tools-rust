//! A BK-tree over a discrete metric (Python dependency: `pybktree`).
//!
//! `whitelist_methods.getErrorCorrectMapping` builds a BK-tree over the
//! whitelist barcodes and queries it once per observed barcode, which turns an
//! all-against-all comparison into a bounded-radius search.
//!
//! Reimplemented rather than adopting the `bk-tree` crate because the caller
//! branches on `len(candidates)`, so which candidates come back — not merely
//! how many — is observable, and a local implementation keeps that under this
//! crate's control and under its differential tests.
//!
//! # The metric
//!
//! `pybktree` requires a metric satisfying the triangle inequality, which is
//! what makes the radius pruning sound. Python passes a *locally redefined*
//! Hamming distance rather than the Cython one, with the comment "Unexpected
//! results with cythonise hamming distance so redefine in python here", and
//! that local version returns `np.inf` for unequal lengths. Infinity is not a
//! metric value the triangle inequality handles, so unequal-length barcodes
//! are simply never within any finite radius — which is the behaviour
//! [`hamming_or_infinite`] reproduces.

use std::collections::BTreeMap;

/// The distance used by [`BkTree`]: Hamming distance, or `None` (conceptually
/// infinite) for unequal lengths.
///
/// Mirrors the `hamming_distance` closure defined inside
/// `whitelist_methods.getErrorCorrectMapping`, which returns `np.inf` rather
/// than attempting to align barcodes of different lengths.
#[must_use]
pub fn hamming_or_infinite(a: &[u8], b: &[u8]) -> Option<u32> {
    crate::distance::edit_distance(a, b)
}

/// A node of the tree: an item plus its children keyed by distance from it.
#[derive(Debug, Clone)]
struct Node {
    item: Vec<u8>,
    children: BTreeMap<u32, Node>,
}

/// A Burkhard-Keller tree over byte strings under [`hamming_or_infinite`]
/// (Python dependency: `pybktree.BKTree`).
///
/// # Examples
///
/// ```
/// use umi_tools_core::bktree::BkTree;
///
/// let tree = BkTree::from_items([b"AAAA".to_vec(), b"CCCC".to_vec()]);
/// let hits = tree.find(b"AAAT", 1);
/// assert_eq!(hits, vec![(1, b"AAAA".to_vec())]);
/// assert!(tree.find(b"GGGG", 1).is_empty());
/// ```
#[derive(Debug, Clone, Default)]
pub struct BkTree {
    root: Option<Node>,
    len: usize,
}

impl BkTree {
    /// Creates an empty tree.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Builds a tree from an iterator of items.
    ///
    /// Items whose distance to the tree is infinite (a different length from
    /// every node on their insertion path) still get inserted, as a child of
    /// the deepest reachable node under a sentinel key, so they remain
    /// findable by an exact-length query later. Python's `pybktree` would
    /// raise on a non-integer distance; this cannot arise from
    /// `getErrorCorrectMapping`, whose whitelist barcodes are all one length,
    /// and the sentinel keeps the structure total rather than partial.
    pub fn from_items(items: impl IntoIterator<Item = Vec<u8>>) -> Self {
        let mut tree = Self::new();
        for item in items {
            tree.add(item);
        }
        tree
    }

    /// The number of items in the tree.
    #[must_use]
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether the tree holds no items.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Distance key used for an item that is not comparable to a node.
    ///
    /// `u32::MAX` cannot collide with a real Hamming distance, which is
    /// bounded by the barcode length.
    const INCOMPARABLE: u32 = u32::MAX;

    /// Inserts `item`.
    pub fn add(&mut self, item: Vec<u8>) {
        self.len += 1;
        let Some(root) = self.root.as_mut() else {
            self.root = Some(Node {
                item,
                children: BTreeMap::new(),
            });
            return;
        };

        let mut node = root;
        loop {
            let key = hamming_or_infinite(&node.item, &item).unwrap_or(Self::INCOMPARABLE);
            if key == 0 {
                // An exact duplicate. `pybktree` stores duplicates as
                // distance-0 children of each other, forming a chain; this
                // keeps that structure so `find` reports each copy.
                // An exact duplicate. `pybktree` stores duplicates as
                // distance-0 children of each other, forming a chain; this
                // keeps that structure so `find` reports each copy.
                match node.children.entry(0) {
                    std::collections::btree_map::Entry::Vacant(slot) => {
                        slot.insert(Node {
                            item,
                            children: BTreeMap::new(),
                        });
                        return;
                    }
                    std::collections::btree_map::Entry::Occupied(slot) => {
                        node = slot.into_mut();
                    }
                }
                continue;
            }
            if let std::collections::btree_map::Entry::Vacant(e) = node.children.entry(key) {
                e.insert(Node {
                    item,
                    children: BTreeMap::new(),
                });
                return;
            }
            node = node.children.get_mut(&key).expect("checked");
        }
    }

    /// Every item within `radius` of `query`, as `(distance, item)` pairs
    /// (Python: `BKTree.find(item, threshold)`).
    ///
    /// Results are sorted by distance and then by item, so the order is
    /// defined. `pybktree` returns them sorted by distance with ties in tree
    /// order, which is an insertion-order artefact; the caller in
    /// `getErrorCorrectMapping` only counts and indexes the results, so a
    /// defined tie order is strictly better than an incidental one. See
    /// `port_reports/umi_tools/05_design.md` § Determinism decisions.
    #[must_use]
    pub fn find(&self, query: &[u8], radius: u32) -> Vec<(u32, Vec<u8>)> {
        let mut out = Vec::new();
        if let Some(root) = self.root.as_ref() {
            Self::search(root, query, radius, &mut out);
        }
        out.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
        out
    }

    fn search(node: &Node, query: &[u8], radius: u32, out: &mut Vec<(u32, Vec<u8>)>) {
        let Some(distance) = hamming_or_infinite(&node.item, query) else {
            // Incomparable to this node: it cannot be a hit, and the triangle
            // inequality gives no bound on its children, so they are all
            // searched.
            for child in node.children.values() {
                Self::search(child, query, radius, out);
            }
            return;
        };

        if distance <= radius {
            out.push((distance, node.item.clone()));
        }
        // Triangle inequality: any item at distance d from this node can be at
        // most `radius` from the query only if `|d - distance| <= radius`.
        let lo = distance.saturating_sub(radius);
        let hi = distance.saturating_add(radius);
        for (&key, child) in &node.children {
            if key == Self::INCOMPARABLE || (key >= lo && key <= hi) {
                Self::search(child, query, radius, out);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bc(s: &str) -> Vec<u8> {
        s.as_bytes().to_vec()
    }

    #[test]
    fn empty_tree_finds_nothing() {
        let tree = BkTree::new();
        assert!(tree.is_empty());
        assert_eq!(tree.len(), 0);
        assert!(tree.find(b"AAAA", 1).is_empty());
    }

    #[test]
    fn finds_exact_and_near_matches() {
        let tree = BkTree::from_items([bc("AAAA"), bc("CCCC"), bc("AAAT")]);
        assert_eq!(tree.len(), 3);
        assert_eq!(tree.find(b"AAAA", 0), vec![(0, bc("AAAA"))]);
        assert_eq!(
            tree.find(b"AAAA", 1),
            vec![(0, bc("AAAA")), (1, bc("AAAT"))]
        );
        assert!(tree.find(b"GGGG", 1).is_empty());
    }

    #[test]
    fn radius_widens_the_result_set() {
        let tree = BkTree::from_items([bc("AAAA"), bc("AAAT"), bc("AATT"), bc("ATTT")]);
        assert_eq!(tree.find(b"AAAA", 1).len(), 2);
        assert_eq!(tree.find(b"AAAA", 2).len(), 3);
        assert_eq!(tree.find(b"AAAA", 3).len(), 4);
    }

    #[test]
    fn unequal_lengths_are_never_within_radius() {
        let tree = BkTree::from_items([bc("AAA"), bc("AAAAA")]);
        assert!(tree.find(b"AAAA", 4).is_empty());
        assert_eq!(tree.find(b"AAA", 0), vec![(0, bc("AAA"))]);
        assert_eq!(tree.find(b"AAAAA", 0), vec![(0, bc("AAAAA"))]);
    }

    #[test]
    fn duplicates_are_all_reported() {
        let tree = BkTree::from_items([bc("AAAA"), bc("AAAA")]);
        assert_eq!(tree.len(), 2);
        assert_eq!(tree.find(b"AAAA", 0).len(), 2);
    }

    #[test]
    fn results_are_sorted_by_distance_then_item() {
        let tree = BkTree::from_items([bc("ACGT"), bc("ACGA"), bc("ACGC"), bc("ACGG")]);
        let hits = tree.find(b"ACGT", 1);
        assert_eq!(hits[0], (0, bc("ACGT")));
        let tail: Vec<Vec<u8>> = hits[1..].iter().map(|(_, i)| i.clone()).collect();
        let mut sorted = tail.clone();
        sorted.sort();
        assert_eq!(tail, sorted);
    }

    /// The pruning must not lose anything relative to a linear scan.
    #[test]
    fn agrees_with_exhaustive_search() {
        let alphabet = b"ACGT";
        let mut items = Vec::new();
        for a in alphabet {
            for b in alphabet {
                for c in alphabet {
                    items.push(vec![*a, *b, *c]);
                }
            }
        }
        let tree = BkTree::from_items(items.clone());

        for query in &items {
            for radius in 0..=3 {
                let mut expected: Vec<(u32, Vec<u8>)> = items
                    .iter()
                    .filter_map(|i| {
                        hamming_or_infinite(i, query)
                            .filter(|d| *d <= radius)
                            .map(|d| (d, i.clone()))
                    })
                    .collect();
                expected.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
                assert_eq!(
                    tree.find(query, radius),
                    expected,
                    "query={:?} radius={radius}",
                    String::from_utf8_lossy(query)
                );
            }
        }
    }
}
