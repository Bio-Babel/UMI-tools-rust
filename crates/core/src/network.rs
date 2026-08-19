//! UMI clustering: the algorithmic core (Python source: `umi_tools/network.py`).
//!
//! Given the UMIs observed at one alignment position and how often each was
//! seen, decide which of them are sequencing errors of which others. Five
//! methods are available, in increasing order of error-awareness:
//!
//! | Method | Behaviour |
//! |---|---|
//! | [`ClusterMethod::Unique`] | every distinct UMI is its own group; no error correction |
//! | [`ClusterMethod::Percentile`] | discard UMIs seen fewer than 1% of the median count |
//! | [`ClusterMethod::Cluster`] | one group per connected component of the "within threshold" graph |
//! | [`ClusterMethod::Adjacency`] | cover each component with the fewest hub UMIs |
//! | [`ClusterMethod::Directional`] | connect `a -> b` only when `count(a) >= 2*count(b) - 1` |
//!
//! # How the Python original dispatches, and why this module does not
//!
//! `UMIClusterer.__init__` selects a method by rebinding three instance
//! attributes to different bound methods:
//!
//! ```text
//! if cluster_method == "adjacency":
//!     self.get_adj_list = self._get_adj_list_adjacency
//!     self.get_connected_components = self._get_connected_components_adjacency
//!     self.get_groups = self._group_adjacency
//! elif ...
//! ```
//!
//! An unrecognised name binds nothing, so construction succeeds and the first
//! call fails with `AttributeError`. The `null` variants also make the
//! intermediate types polymorphic: `_get_adj_list_null` returns `None`, and
//! `_get_connected_components_null` returns the flat UMI list rather than a
//! list of components.
//!
//! This module resolves all of that into [`ClusterMethod`], [`AdjList`] and
//! [`Clusters`], so the method is validated at construction and the two
//! intermediate shapes are distinguishable at compile time. See
//! `port_reports/umi_tools/05_design.md` deviation D2.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use crate::distance::within_threshold;
use crate::error::{PortError, Result};
use crate::numeric::median;
use crate::utilities;

/// A UMI or cell barcode: raw bytes, because that is what the BAM path
/// produces and what [`crate::distance::edit_distance`] compares.
pub type Umi = Vec<u8>;

/// UMI to observed read count, as the Python original's `counts` dictionary.
///
/// A `BTreeMap` rather than a `HashMap`: iteration order is observable in
/// several places (component ordering, group ordering), and Python's own order
/// there is an insertion-order artefact. Sorted order makes the Rust result
/// reproducible run to run. See `05_design.md` § Determinism decisions.
pub type CountMap = BTreeMap<Umi, u64>;

/// UMI to the UMIs it is connected to, as the Python original's `adj_list`.
pub type AdjMap = BTreeMap<Umi, Vec<Umi>>;

/// Which clustering method to apply (Python source: the `cluster_method`
/// string `UMIClusterer.__init__` branches on).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ClusterMethod {
    /// No error correction: every distinct UMI is its own group.
    Unique,
    /// Retain only UMIs whose count exceeds 1% of the median count.
    Percentile,
    /// One group per connected component of the within-threshold graph.
    Cluster,
    /// Cover each component with the fewest "hub" UMIs.
    Adjacency,
    /// Directional: connect `a -> b` only when `count(a) >= 2*count(b) - 1`.
    Directional,
}

impl ClusterMethod {
    /// The method's name as spelled on the command line and in the Python
    /// source.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            ClusterMethod::Unique => "unique",
            ClusterMethod::Percentile => "percentile",
            ClusterMethod::Cluster => "cluster",
            ClusterMethod::Adjacency => "adjacency",
            ClusterMethod::Directional => "directional",
        }
    }

    /// Every method, in the order `--method`'s help text lists them.
    #[must_use]
    pub fn all() -> [ClusterMethod; 5] {
        [
            ClusterMethod::Unique,
            ClusterMethod::Percentile,
            ClusterMethod::Cluster,
            ClusterMethod::Adjacency,
            ClusterMethod::Directional,
        ]
    }
}

impl std::str::FromStr for ClusterMethod {
    type Err = PortError;

    /// # Errors
    ///
    /// Returns [`PortError::UnknownClusterMethod`] for any other string. The
    /// Python original accepts anything here and fails later with
    /// `AttributeError` on the first call — deviation D2.
    fn from_str(s: &str) -> Result<Self> {
        match s {
            "unique" => Ok(ClusterMethod::Unique),
            "percentile" => Ok(ClusterMethod::Percentile),
            "cluster" => Ok(ClusterMethod::Cluster),
            "adjacency" => Ok(ClusterMethod::Adjacency),
            "directional" => Ok(ClusterMethod::Directional),
            other => Err(PortError::UnknownClusterMethod(other.to_string())),
        }
    }
}

impl std::fmt::Display for ClusterMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// The adjacency structure a method produces.
///
/// The `unique` and `percentile` methods build no graph at all; Python models
/// that by returning `None` from `_get_adj_list_null`, which makes the
/// parameter's type `dict | None`. This enum makes the absence explicit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AdjList {
    /// A real graph, from `_get_adj_list_adjacency` or `_get_adj_list_directional`.
    Graph(AdjMap),
    /// No graph, from `_get_adj_list_null`.
    None,
}

impl AdjList {
    /// The underlying graph, if there is one.
    #[must_use]
    pub fn graph(&self) -> Option<&AdjMap> {
        match self {
            AdjList::Graph(g) => Some(g),
            AdjList::None => None,
        }
    }
}

/// The clustering stage's output.
///
/// `_get_connected_components_adjacency` returns a list of components, while
/// `_get_connected_components_null` returns the *flat* UMI list unchanged.
/// Python leaves these indistinguishable; `_group_percentile` relies on
/// receiving the flat form and passes it straight to `_get_best_percentile`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Clusters {
    /// Connected components, each already sorted.
    Components(Vec<Vec<Umi>>),
    /// The flat UMI list, passed through unchanged.
    Flat(Vec<Umi>),
}

/// Every UMI reachable from `node`, including `node` itself
/// (Python source: `network.breadth_first_search`).
///
/// Python accumulates into a `set` and returns it; this returns a
/// [`BTreeSet`] so the traversal result is ordered, which every caller wants
/// anyway — `_get_connected_components_adjacency` immediately sorts it.
///
/// A UMI absent from `adj` contributes no neighbours, matching Python, where
/// `adj_list` is built with an entry for every UMI so the lookup cannot miss.
///
/// # Examples
///
/// ```
/// use std::collections::BTreeMap;
/// use umi_tools_core::network::breadth_first_search;
///
/// let mut adj = BTreeMap::new();
/// adj.insert(b"A".to_vec(), vec![b"B".to_vec()]);
/// adj.insert(b"B".to_vec(), vec![b"A".to_vec()]);
/// adj.insert(b"C".to_vec(), vec![]);
///
/// let component = breadth_first_search(b"A", &adj);
/// assert_eq!(component.len(), 2);
/// assert!(component.contains(&b"B".to_vec()));
/// ```
#[must_use]
pub fn breadth_first_search(node: &[u8], adj: &AdjMap) -> BTreeSet<Umi> {
    let mut searched: BTreeSet<Umi> = BTreeSet::new();
    let mut queue: VecDeque<Umi> = VecDeque::new();
    searched.insert(node.to_vec());
    queue.push_back(node.to_vec());

    while let Some(current) = queue.pop_front() {
        if let Some(neighbours) = adj.get(&current) {
            for next in neighbours {
                if searched.insert(next.clone()) {
                    queue.push_back(next.clone());
                }
            }
        }
    }
    searched
}

/// Recursion depth at which [`recursive_search`] gives up.
///
/// Python's `network.py` calls `sys.setrecursionlimit(10000)` at import, and
/// `breadth_first_search_recursive` catches the resulting `RecursionError`.
/// Rust cannot catch a stack overflow, so the depth is tracked explicitly and
/// the same limit applied.
pub const RECURSION_LIMIT: usize = 10_000;

/// Depth-first traversal from `node` (Python source: `network.recursive_search`).
///
/// The Python original stores its accumulator on the *function object*
/// (`recursive_search.component = set((node,))`), which is why calling it
/// directly without first setting that attribute raises `AttributeError`.
/// Rust has no function attributes, so the accumulator is threaded through a
/// private helper while this public wrapper keeps the Python signature.
///
/// Returns `None` if the traversal exceeds [`RECURSION_LIMIT`], mirroring the
/// `RecursionError` Python raises at its own recursion limit.
///
/// # Examples
///
/// ```
/// use std::collections::BTreeMap;
/// use umi_tools_core::network::recursive_search;
///
/// let mut adj = BTreeMap::new();
/// adj.insert(b"A".to_vec(), vec![b"B".to_vec()]);
/// adj.insert(b"B".to_vec(), vec![]);
/// assert_eq!(recursive_search(b"A", &adj).unwrap().len(), 2);
/// ```
#[must_use]
pub fn recursive_search(node: &[u8], adj: &AdjMap) -> Option<BTreeSet<Umi>> {
    let mut component: BTreeSet<Umi> = BTreeSet::new();
    component.insert(node.to_vec());
    if recursive_search_inner(node, adj, &mut component, 0) {
        Some(component)
    } else {
        None
    }
}

/// Depth-first traversal driven by an explicit stack.
///
/// Returns `false` once the traversal would exceed [`RECURSION_LIMIT`]
/// frames, which is the condition Python surfaces as `RecursionError`.
///
/// The Python original really does recurse, and really can blow its stack —
/// which is why `network.py` raises the interpreter's limit to 10 000 at
/// import and why `breadth_first_search_recursive` catches the error. Rust
/// cannot catch a stack overflow: it aborts the process. So the traversal
/// order and the depth accounting are reproduced exactly while the frames
/// live on the heap, making the limit a *reported* condition rather than a
/// crash. This is an implementation adaptation, not a behaviour change — the
/// same graphs succeed, the same graphs report the limit.
fn recursive_search_inner(
    start: &[u8],
    adj: &AdjMap,
    component: &mut BTreeSet<Umi>,
    start_depth: usize,
) -> bool {
    // (node, depth) pairs still to expand, deepest last — the order a
    // depth-first recursion would visit them in.
    let mut stack: Vec<(Umi, usize)> = vec![(start.to_vec(), start_depth)];

    while let Some((node, depth)) = stack.pop() {
        if depth >= RECURSION_LIMIT {
            return false;
        }
        let Some(neighbours) = adj.get(&node) else {
            continue;
        };
        // Python filters against the accumulator before recursing and inserts
        // each survivor as it goes. Inserting first is equivalent in result,
        // because the accumulator is a set, and it avoids descending twice
        // into a node reachable twice from one parent.
        for child in neighbours.iter().rev() {
            if component.insert(child.clone()) {
                stack.push((child.clone(), depth + 1));
            }
        }
    }
    true
}

/// Depth-first traversal that falls back to breadth-first on deep graphs
/// (Python source: `network.breadth_first_search_recursive`).
///
/// Python wraps `recursive_search` in `try/except RecursionError` and retries
/// with the iterative [`breadth_first_search`]; the same fallback happens here
/// when [`recursive_search`] hits [`RECURSION_LIMIT`].
///
/// This function is dead code upstream — the only call site in
/// `_get_connected_components_adjacency` is commented out, with a `TS: TO DO:
/// Work out why recursive function doesn't lead to same final output` note. It
/// is ported because it is part of the module's public surface.
#[must_use]
pub fn breadth_first_search_recursive(node: &[u8], adj: &AdjMap) -> BTreeSet<Umi> {
    if let Some(component) = recursive_search(node, adj) {
        component
    } else {
        utilities::info("Recursion Error: maximum recursion depth exceeded");
        breadth_first_search(node, adj)
    }
}

/// The members of `cluster` left after removing `nodes` and everything
/// adjacent to them (Python source: `network.remove_umis`).
///
/// Used by [`UmiClusterer::get_best_min_account`] to test whether a candidate
/// set of hub UMIs accounts for the whole cluster.
#[must_use]
pub fn remove_umis(adj: &AdjMap, cluster: &[Umi], nodes: &[Umi]) -> BTreeSet<Umi> {
    let mut to_remove: BTreeSet<&Umi> = BTreeSet::new();
    for node in nodes {
        if let Some(neighbours) = adj.get(node) {
            to_remove.extend(neighbours.iter());
        }
        to_remove.insert(node);
    }
    cluster
        .iter()
        .filter(|u| !to_remove.contains(*u))
        .cloned()
        .collect()
}

/// Split points dividing a UMI of `umi_length` bases into `idx_size`
/// approximately equal substrings (Python source: `network.get_substr_slices`).
///
/// Reproduces Python's `divmod` split exactly, including its degenerate
/// behaviour when `idx_size > umi_length`: the remainder-sized slices come
/// first, and the leftovers are empty slices at the end. For example
/// `get_substr_slices(1, 3)` yields `[(0, 1), (1, 1), (1, 1)]` — two empty
/// slices at the same offset, which then collapse into one key in
/// [`build_substr_idx`]'s map, exactly as they do in Python's `defaultdict`.
///
/// # Examples
///
/// ```
/// use umi_tools_core::network::get_substr_slices;
/// assert_eq!(get_substr_slices(6, 2), vec![(0, 3), (3, 6)]);
/// assert_eq!(get_substr_slices(5, 2), vec![(0, 3), (3, 5)]);
/// assert_eq!(get_substr_slices(1, 3), vec![(0, 1), (1, 1), (1, 1)]);
/// ```
#[must_use]
pub fn get_substr_slices(umi_length: usize, idx_size: usize) -> Vec<(usize, usize)> {
    if idx_size == 0 {
        return Vec::new();
    }
    let cs = umi_length / idx_size;
    let r = umi_length % idx_size;
    let mut offset = 0;
    let mut slices = Vec::with_capacity(idx_size);
    for i in 0..idx_size {
        let size = if i < r { cs + 1 } else { cs };
        slices.push((offset, offset + size));
        offset += size;
    }
    slices
}

/// A substring index: for each slice, the UMIs sharing each substring
/// (Python source: the `defaultdict(lambda: defaultdict(set))` built by
/// `network.build_substr_idx`).
///
/// The pigeonhole principle underpins it: two UMIs within `min_edit`
/// substitutions must agree exactly on at least one of `min_edit + 1`
/// disjoint slices, so candidate neighbours can be found without comparing
/// every pair.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SubstrIdx {
    /// Slice bounds to substring to the UMIs carrying it.
    entries: BTreeMap<(usize, usize), BTreeMap<Vec<u8>, BTreeSet<Umi>>>,
}

impl SubstrIdx {
    /// The slices this index is built over, in insertion order.
    #[must_use]
    pub fn slices(&self) -> Vec<(usize, usize)> {
        self.entries.keys().copied().collect()
    }

    /// The UMIs sharing `umi`'s substring at each indexed slice, unioned.
    #[must_use]
    pub fn candidates(&self, umi: &[u8]) -> BTreeSet<Umi> {
        let mut out = BTreeSet::new();
        for (&(start, end), map) in &self.entries {
            let sub = &umi[start.min(umi.len())..end.min(umi.len())];
            if let Some(members) = map.get(sub) {
                out.extend(members.iter().cloned());
            }
        }
        out
    }

    /// Whether the index holds no slices.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// Builds the substring index used to shortlist neighbour candidates
/// (Python source: `network.build_substr_idx`).
///
/// `min_edit` is the edit-distance threshold; the UMI is split into
/// `min_edit + 1` slices.
///
/// Duplicate slices — which `get_substr_slices` produces when
/// `idx_size > umi_length` — collapse to a single map entry here, exactly as
/// they do in Python's `substr_idx[idx]` dictionary keyed by the slice tuple.
#[must_use]
pub fn build_substr_idx(umis: &[Umi], umi_length: usize, min_edit: u32) -> SubstrIdx {
    let mut entries: BTreeMap<(usize, usize), BTreeMap<Vec<u8>, BTreeSet<Umi>>> = BTreeMap::new();
    let slices = get_substr_slices(umi_length, min_edit as usize + 1);
    for (start, end) in slices {
        let per_slice = entries.entry((start, end)).or_default();
        for umi in umis {
            let sub = umi[start.min(umi.len())..end.min(umi.len())].to_vec();
            per_slice.entry(sub).or_default().insert(umi.clone());
        }
    }
    SubstrIdx { entries }
}

/// Candidate neighbour pairs from the substring index
/// (Python source: `network.iter_nearest_neighbours`).
///
/// Python's generator enumerates `umis` with `enumerate(umis, 1)` and calls
/// `neighbours.difference_update(umis[:i])`, which removes every UMI at or
/// before the current position — including the current UMI itself. Each
/// unordered pair is therefore yielded exactly once, in the order the outer
/// loop reaches its earlier member.
///
/// Python yields lazily from a `set`, so the order *within* one UMI's
/// neighbour block is an artefact. This returns an eager `Vec` in sorted
/// order; every consumer treats the stream as an unordered pair set, and the
/// differential tests compare it as one. See `05_design.md` § Determinism
/// decisions.
#[must_use]
pub fn iter_nearest_neighbours(umis: &[Umi], idx: &SubstrIdx) -> Vec<(Umi, Umi)> {
    let mut out = Vec::new();
    for (i, umi) in umis.iter().enumerate() {
        let seen: BTreeSet<&Umi> = umis[..=i].iter().collect();
        let mut neighbours = idx.candidates(umi);
        neighbours.retain(|n| !seen.contains(n));
        for nbr in neighbours {
            out.push((umi.clone(), nbr));
        }
    }
    out
}

/// Above this many UMIs, the substring index replaces exhaustive pairwise
/// comparison (Python source: the `if len(umis) > 25:` branch in both
/// `_get_adj_list_adjacency` and `_get_adj_list_directional`).
///
/// The two paths are not exactly equivalent: the substring index is an
/// approximation that can miss a pair whose only agreement falls outside the
/// indexed slices when `umi_length < min_edit + 1`. This threshold is
/// mirrored rather than tuned, so the port reproduces the Python result
/// including where the approximation bites.
pub const SUBSTR_INDEX_THRESHOLD: usize = 25;

fn candidate_pairs(umis: &[Umi], threshold: u32) -> Vec<(Umi, Umi)> {
    if umis.len() > SUBSTR_INDEX_THRESHOLD {
        let umi_length = umis[0].len();
        let idx = build_substr_idx(umis, umi_length, threshold);
        iter_nearest_neighbours(umis, &idx)
    } else {
        let mut pairs = Vec::new();
        for i in 0..umis.len() {
            for j in (i + 1)..umis.len() {
                pairs.push((umis[i].clone(), umis[j].clone()));
            }
        }
        pairs
    }
}

/// Clusters the UMIs observed at one alignment position
/// (Python source: `network.UMIClusterer`).
///
/// Construct with a [`ClusterMethod`], then call [`UmiClusterer::cluster`]
/// once per position. The clusterer accumulates per-run statistics across
/// calls — which is why `cluster` takes `&mut self` where Python's `__call__`
/// looks pure — and `dedup`/`group` report them at the end of a run.
///
/// # Examples
///
/// ```
/// use std::collections::BTreeMap;
/// use umi_tools_core::network::{ClusterMethod, UmiClusterer};
///
/// let counts: BTreeMap<Vec<u8>, u64> =
///     [(b"AAAA".to_vec(), 100), (b"AAAT".to_vec(), 2)].into_iter().collect();
///
/// let mut clusterer = UmiClusterer::new(ClusterMethod::Directional);
/// let groups = clusterer.cluster(&counts, 1)?;
///
/// // AAAT is one substitution from AAAA and 100 >= 2*2 - 1, so it is absorbed.
/// assert_eq!(groups, vec![vec![b"AAAA".to_vec(), b"AAAT".to_vec()]]);
/// assert_eq!(clusterer.positions(), 1);
/// # Ok::<(), umi_tools_core::PortError>(())
/// ```
#[derive(Debug, Clone)]
pub struct UmiClusterer {
    method: ClusterMethod,
    max_umis_per_position: usize,
    total_umis_per_position: u64,
    positions: u64,
}

impl UmiClusterer {
    /// Creates a clusterer for `method`.
    #[must_use]
    pub fn new(method: ClusterMethod) -> Self {
        Self {
            method,
            max_umis_per_position: 0,
            total_umis_per_position: 0,
            positions: 0,
        }
    }

    /// Creates a clusterer from a method name.
    ///
    /// # Errors
    ///
    /// Returns [`PortError::UnknownClusterMethod`] for an unrecognised name.
    /// Python accepts any string here and fails on first use — deviation D2 in
    /// `05_design.md`.
    pub fn from_name(name: &str) -> Result<Self> {
        Ok(Self::new(name.parse()?))
    }

    /// The method this clusterer applies.
    #[must_use]
    pub fn method(&self) -> ClusterMethod {
        self.method
    }

    /// How many positions have been clustered so far
    /// (Python source: `UMIClusterer.positions`).
    #[must_use]
    pub fn positions(&self) -> u64 {
        self.positions
    }

    /// Total distinct UMIs summed over every position clustered so far
    /// (Python source: `UMIClusterer.total_umis_per_position`).
    #[must_use]
    pub fn total_umis_per_position(&self) -> u64 {
        self.total_umis_per_position
    }

    /// The largest number of distinct UMIs seen at any one position
    /// (Python source: `UMIClusterer.max_umis_per_position`).
    #[must_use]
    pub fn max_umis_per_position(&self) -> usize {
        self.max_umis_per_position
    }

    /// Builds the adjacency structure for this method
    /// (Python source: the `get_adj_list` attribute).
    #[must_use]
    pub fn get_adj_list(&self, umis: &[Umi], counts: &CountMap, threshold: u32) -> AdjList {
        match self.method {
            ClusterMethod::Unique | ClusterMethod::Percentile => AdjList::None,
            ClusterMethod::Cluster | ClusterMethod::Adjacency => {
                AdjList::Graph(Self::adj_list_adjacency(umis, threshold))
            }
            ClusterMethod::Directional => {
                AdjList::Graph(Self::adj_list_directional(umis, counts, threshold))
            }
        }
    }

    /// Connects every pair of UMIs within `threshold` substitutions
    /// (Python source: `UMIClusterer._get_adj_list_adjacency`).
    ///
    /// The edge is symmetric: both UMIs gain the other as a neighbour.
    #[must_use]
    pub fn adj_list_adjacency(umis: &[Umi], threshold: u32) -> AdjMap {
        let mut adj: AdjMap = umis.iter().map(|u| (u.clone(), Vec::new())).collect();
        for (a, b) in candidate_pairs(umis, threshold) {
            if within_threshold(&a, &b, threshold) {
                adj.entry(a.clone()).or_default().push(b.clone());
                adj.entry(b).or_default().push(a);
            }
        }
        adj
    }

    /// Connects `a -> b` when `a` is within `threshold` substitutions of `b`
    /// **and** `count(a) >= 2*count(b) - 1`
    /// (Python source: `UMIClusterer._get_adj_list_directional`).
    ///
    /// The count condition encodes the model that a PCR/sequencing error is
    /// seen at most about half as often as the true UMI it derives from, so
    /// the edge is directional: an error points at its parent, not vice versa.
    /// Both directions are tested for every pair, so two UMIs with near-equal
    /// counts connect both ways and end up in one component.
    ///
    /// A missing count is treated as zero, matching Python's `counts[umi]` on
    /// a `dict` built from the same UMI list — where a genuine miss would be a
    /// `KeyError` that cannot occur, since the caller derives `umis` from
    /// `counts`.
    #[must_use]
    pub fn adj_list_directional(umis: &[Umi], counts: &CountMap, threshold: u32) -> AdjMap {
        let mut adj: AdjMap = umis.iter().map(|u| (u.clone(), Vec::new())).collect();
        let count_of = |u: &Umi| -> i128 { i128::from(counts.get(u).copied().unwrap_or(0)) };

        for (a, b) in candidate_pairs(umis, threshold) {
            if !within_threshold(&a, &b, threshold) {
                continue;
            }
            let (ca, cb) = (count_of(&a), count_of(&b));
            // `(count * 2) - 1` is computed in i128 so a count near u64::MAX
            // cannot overflow the way it silently would in a narrower type.
            // Python's ints are unbounded, so this matches its arithmetic.
            if ca >= cb * 2 - 1 {
                adj.entry(a.clone()).or_default().push(b.clone());
            }
            if cb >= ca * 2 - 1 {
                adj.entry(b).or_default().push(a);
            }
        }
        adj
    }

    /// Groups UMIs into connected components, or passes them through
    /// (Python source: the `get_connected_components` attribute).
    #[must_use]
    pub fn get_connected_components(
        &self,
        umis: &[Umi],
        adj: &AdjList,
        counts: &CountMap,
    ) -> Clusters {
        match adj {
            AdjList::None => Clusters::Flat(umis.to_vec()),
            AdjList::Graph(graph) => {
                Clusters::Components(Self::connected_components_adjacency(umis, graph, counts))
            }
        }
    }

    /// Finds the connected components of `graph`
    /// (Python source: `UMIClusterer._get_connected_components_adjacency`).
    ///
    /// Components are discovered by walking the graph's UMIs from most to
    /// least abundant, so the most abundant unvisited UMI seeds each
    /// component; each component is then sorted before being returned. Python
    /// sorts with `sorted(graph, key=lambda x: counts[x], reverse=True)`,
    /// which is stable — ties keep the graph's own iteration order — so this
    /// uses `sort_by` rather than `sort_unstable_by` to match.
    ///
    /// `umis` supplies the order ties are broken in. Python sorts the
    /// adjacency **dictionary** — `sorted(graph, key=lambda x: counts[x],
    /// reverse=True)` — and that dictionary was built as `{umi: [] for umi in
    /// umis}`, so its iteration order is the caller's UMI order. Since
    /// `sorted` is stable, equal-count UMIs seed components in that order, and
    /// the components come out in that order too. That decides the order
    /// `dedup` writes its reads, so it has to be threaded through rather than
    /// recovered from the map's own (sorted) order.
    ///
    /// A UMI in `umis` but absent from `graph` is skipped; in practice the two
    /// always hold the same keys.
    #[must_use]
    pub fn connected_components_adjacency(
        umis: &[Umi],
        graph: &AdjMap,
        counts: &CountMap,
    ) -> Vec<Vec<Umi>> {
        let mut nodes: Vec<&Umi> = umis.iter().filter(|u| graph.contains_key(*u)).collect();
        nodes.sort_by(|a, b| {
            counts
                .get(*b)
                .copied()
                .unwrap_or(0)
                .cmp(&counts.get(*a).copied().unwrap_or(0))
        });

        let mut found: BTreeSet<Umi> = BTreeSet::new();
        let mut components: Vec<Vec<Umi>> = Vec::new();
        for node in nodes {
            if found.contains(node) {
                continue;
            }
            let component: Vec<Umi> = breadth_first_search(node, graph).into_iter().collect();
            found.extend(component.iter().cloned());
            components.push(component);
        }
        components
    }

    /// The fewest highest-count UMIs whose neighbourhoods cover `cluster`
    /// (Python source: `UMIClusterer._get_best_min_account`).
    ///
    /// Greedily takes UMIs in descending count order until removing them and
    /// their neighbours leaves nothing. Returns `None` when no prefix covers
    /// the cluster, matching Python, whose loop falls off the end and
    /// implicitly returns `None`. That cannot happen for a genuinely connected
    /// component — taking all but one node necessarily covers the last, since
    /// it must be adjacent to one of them — so `None` signals a caller that
    /// passed a disconnected "cluster".
    #[must_use]
    pub fn get_best_min_account(
        cluster: &[Umi],
        adj: &AdjMap,
        counts: &CountMap,
    ) -> Option<Vec<Umi>> {
        if cluster.len() == 1 {
            return Some(cluster.to_vec());
        }
        let mut sorted: Vec<Umi> = cluster.to_vec();
        sorted.sort_by(|a, b| {
            counts
                .get(b)
                .copied()
                .unwrap_or(0)
                .cmp(&counts.get(a).copied().unwrap_or(0))
        });

        for i in 0..sorted.len().saturating_sub(1) {
            let candidate = &sorted[..=i];
            if remove_umis(adj, cluster, candidate).is_empty() {
                return Some(candidate.to_vec());
            }
        }
        None
    }

    /// UMIs whose count exceeds 1% of the median count
    /// (Python source: `UMIClusterer._get_best_percentile`).
    ///
    /// The median is taken over **every** count in `counts`, not only those
    /// in `cluster` — which is the same thing for the percentile method, whose
    /// "cluster" is the full UMI list, but is worth stating because the
    /// parameter name suggests otherwise.
    ///
    /// A single-element cluster is returned unfiltered, so a lone UMI is never
    /// discarded however small its count.
    #[must_use]
    pub fn get_best_percentile(cluster: &[Umi], counts: &CountMap) -> Vec<Umi> {
        if cluster.len() == 1 {
            return cluster.to_vec();
        }
        let all: Vec<f64> = counts.values().map(|&c| c as f64).collect();
        let threshold = median(&all) / 100.0;
        cluster
            .iter()
            .filter(|u| counts.get(*u).copied().unwrap_or(0) as f64 > threshold)
            .cloned()
            .collect()
    }

    /// Turns clusters into output groups, parent UMI first
    /// (Python source: the `get_groups` attribute).
    ///
    /// # Errors
    ///
    /// Returns [`PortError::Value`] if the clustering stage produced a shape
    /// this method cannot consume — which cannot happen through
    /// [`UmiClusterer::cluster`], since it pairs the two stages itself.
    pub fn get_groups(
        &self,
        clusters: &Clusters,
        adj: &AdjList,
        counts: &CountMap,
    ) -> Result<Vec<Vec<Umi>>> {
        match (self.method, clusters) {
            (ClusterMethod::Unique, Clusters::Flat(umis)) => Ok(Self::group_unique(umis)),
            (ClusterMethod::Percentile, Clusters::Flat(umis)) => {
                Ok(Self::group_percentile(umis, counts))
            }
            (ClusterMethod::Cluster, Clusters::Components(components)) => {
                Ok(Self::group_cluster(components, counts))
            }
            (ClusterMethod::Directional, Clusters::Components(components)) => {
                Ok(Self::group_directional(components, counts))
            }
            (ClusterMethod::Adjacency, Clusters::Components(components)) => {
                let graph = adj.graph().ok_or_else(|| {
                    PortError::value("adjacency grouping requires an adjacency list")
                })?;
                Ok(Self::group_adjacency(components, graph, counts))
            }
            (method, _) => Err(PortError::value(format!(
                "cluster shape does not match method {method}"
            ))),
        }
    }

    /// One group per UMI (Python source: `UMIClusterer._group_unique`).
    ///
    /// Python special-cases a single UMI as `groups = [clusters]` rather than
    /// the comprehension, which produces the same `[[umi]]` either way.
    #[must_use]
    pub fn group_unique(umis: &[Umi]) -> Vec<Vec<Umi>> {
        umis.iter().map(|u| vec![u.clone()]).collect()
    }

    /// One group per retained UMI (Python source: `UMIClusterer._group_percentile`).
    ///
    /// Grouping is not really meaningful for the percentile method; the Python
    /// original says as much in its docstring and returns the surviving UMIs in
    /// the same shape the other methods use.
    #[must_use]
    pub fn group_percentile(umis: &[Umi], counts: &CountMap) -> Vec<Vec<Umi>> {
        Self::get_best_percentile(umis, counts)
            .into_iter()
            .map(|u| vec![u])
            .collect()
    }

    /// Each component sorted by descending count
    /// (Python source: `UMIClusterer._group_cluster`).
    #[must_use]
    pub fn group_cluster(components: &[Vec<Umi>], counts: &CountMap) -> Vec<Vec<Umi>> {
        components
            .iter()
            .map(|component| {
                let mut sorted = component.clone();
                sorted.sort_by(|a, b| {
                    counts
                        .get(b)
                        .copied()
                        .unwrap_or(0)
                        .cmp(&counts.get(a).copied().unwrap_or(0))
                });
                sorted
            })
            .collect()
    }

    /// Each component sorted by descending count, skipping already-emitted
    /// UMIs (Python source: `UMIClusterer._group_directional`).
    ///
    /// The `observed` set is carried across components in the Python original.
    /// Because components are disjoint by construction it never actually
    /// excludes anything; it is mirrored anyway so that a caller passing
    /// overlapping "components" gets Python's behaviour rather than a
    /// divergence.
    #[must_use]
    pub fn group_directional(components: &[Vec<Umi>], counts: &CountMap) -> Vec<Vec<Umi>> {
        let mut observed: BTreeSet<Umi> = BTreeSet::new();
        let mut groups = Vec::with_capacity(components.len());
        for component in components {
            if component.len() == 1 {
                groups.push(component.clone());
                observed.extend(component.iter().cloned());
                continue;
            }
            let mut sorted = component.clone();
            sorted.sort_by(|a, b| {
                counts
                    .get(b)
                    .copied()
                    .unwrap_or(0)
                    .cmp(&counts.get(a).copied().unwrap_or(0))
            });
            let mut group = Vec::with_capacity(sorted.len());
            for node in sorted {
                if observed.insert(node.clone()) {
                    group.push(node);
                }
            }
            groups.push(group);
        }
        groups
    }

    /// One group per hub UMI, each followed by its unclaimed neighbours
    /// (Python source: `UMIClusterer._group_adjacency`).
    ///
    /// Unlike the directional method, `observed` is reset for each component,
    /// so a UMI can appear in two groups from two different components — which
    /// cannot happen, since components are disjoint, but is mirrored exactly.
    ///
    /// The trailing members of each group come from `list(connected_nodes -
    /// observed)` in Python, whose order is a `set`-iteration artefact. They
    /// are emitted in sorted order here; the leading hub UMI, which is the
    /// group's reported parent, is positionally exact either way. See
    /// `05_design.md` § Determinism decisions.
    #[must_use]
    pub fn group_adjacency(
        components: &[Vec<Umi>],
        adj: &AdjMap,
        counts: &CountMap,
    ) -> Vec<Vec<Umi>> {
        let mut groups = Vec::new();
        for component in components {
            if component.len() == 1 {
                groups.push(component.clone());
                continue;
            }
            let mut observed: BTreeSet<Umi> = BTreeSet::new();
            let Some(lead_umis) = Self::get_best_min_account(component, adj, counts) else {
                // Unreachable for a connected component; mirrors Python's
                // implicit `None`, which would raise `TypeError` on the next
                // line. Emitting the component unchanged is the closest
                // non-panicking behaviour.
                groups.push(component.clone());
                continue;
            };
            observed.extend(lead_umis.iter().cloned());

            for lead in &lead_umis {
                let connected: BTreeSet<Umi> =
                    adj.get(lead).into_iter().flatten().cloned().collect();
                let mut group = vec![lead.clone()];
                group.extend(connected.difference(&observed).cloned());
                groups.push(group);
                observed.extend(connected);
            }
        }
        groups
    }

    /// Clusters one position's UMIs in a caller-supplied order
    /// (Python source: `UMIClusterer.__call__`, whose `umis = list(umis.keys())`
    /// takes the counts dictionary's **insertion** order).
    ///
    /// The order matters, and not only cosmetically. It decides which UMI
    /// seeds each connected component when counts tie, and therefore the
    /// order the groups come out in — which becomes the order `dedup` and
    /// `group` write their reads. The package's own test suite compares that
    /// output byte for byte, so reproducing Python's order is what makes those
    /// comparisons pass.
    ///
    /// [`cluster`](Self::cluster) uses the counts map's own (sorted) order,
    /// which is deterministic but not Python's; a caller that has the original
    /// observation order — [`crate::bundles::insertion_order`] provides it —
    /// should pass it here instead.
    ///
    /// # Errors
    ///
    /// As [`cluster`](Self::cluster).
    pub fn cluster_ordered(
        &mut self,
        umis: &[Umi],
        counts: &CountMap,
        threshold: u32,
    ) -> Result<Vec<Vec<Umi>>> {
        self.positions += 1;
        self.total_umis_per_position += umis.len() as u64;
        self.max_umis_per_position = self.max_umis_per_position.max(umis.len());

        let (Some(min_len), Some(max_len)) = (
            umis.iter().map(Vec::len).min(),
            umis.iter().map(Vec::len).max(),
        ) else {
            return Err(PortError::value("max() arg is an empty sequence"));
        };
        if min_len != max_len {
            return Err(PortError::assertion(format!(
                "not all umis are the same length(!):  {min_len} - {max_len}"
            )));
        }

        let adj = self.get_adj_list(umis, counts, threshold);
        let clusters = self.get_connected_components(umis, &adj, counts);
        self.get_groups(&clusters, &adj, counts)
    }

    /// Clusters one position's UMIs (Python source: `UMIClusterer.__call__`).
    ///
    /// Each returned inner vector is one group of UMIs judged to come from the
    /// same original molecule, with the parent UMI first.
    ///
    /// Statistics ([`positions`](Self::positions),
    /// [`total_umis_per_position`](Self::total_umis_per_position),
    /// [`max_umis_per_position`](Self::max_umis_per_position)) are updated on
    /// every call, before the length check — matching Python, which increments
    /// them ahead of its `assert`.
    ///
    /// # Errors
    ///
    /// - [`PortError::Value`] if `counts` is empty. Python reaches
    ///   `max(len_umis)` on an empty list and raises
    ///   `ValueError: max() arg is an empty sequence`.
    /// - [`PortError::Assertion`] if the UMIs are not all the same length,
    ///   matching Python's `assert max(len_umis) == min(len_umis)` and its
    ///   message.
    pub fn cluster(&mut self, counts: &CountMap, threshold: u32) -> Result<Vec<Vec<Umi>>> {
        let umis: Vec<Umi> = counts.keys().cloned().collect();
        self.cluster_ordered(&umis, counts, threshold)
    }
}

/// The result of deduplicating one bundle of reads
/// (Python source: the 3-tuple `ReadDeduplicator.__call__` returns).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DedupResult {
    /// The UMI representing each retained group, in group order.
    pub final_umis: Vec<Umi>,
    /// Total read count behind each retained group, positionally aligned with
    /// [`final_umis`](Self::final_umis).
    pub umi_counts: Vec<u64>,
}

/// Applies [`UmiClusterer`] to a bundle of reads, optionally filtering the
/// result against a UMI whitelist (Python source: `network.ReadDeduplicator`).
///
/// The Python original also returns the representative *reads*; here the
/// caller looks them up from its own bundle using
/// [`final_umis`](DedupResult::final_umis), so this type stays independent of
/// the alignment-record representation and of `crates/cli`'s I/O.
#[derive(Debug, Clone)]
pub struct ReadDeduplicator {
    clusterer: UmiClusterer,
    umi_whitelist: Option<BTreeSet<Umi>>,
    umi_whitelist_counts: BTreeMap<String, u64>,
}

/// The counter key Python uses for reads whose group leader is off-whitelist.
pub const NON_WHITELIST_UMI: &str = "Non-whitelist UMI";

impl ReadDeduplicator {
    /// Creates a deduplicator for `method`, with no UMI whitelist.
    #[must_use]
    pub fn new(method: ClusterMethod) -> Self {
        Self {
            clusterer: UmiClusterer::new(method),
            umi_whitelist: None,
            umi_whitelist_counts: BTreeMap::new(),
        }
    }

    /// Creates a deduplicator that discards any group whose leading UMI is not
    /// in `whitelist` (Python source: the `options.filter_umi` branch of
    /// `ReadDeduplicator.__init__`).
    #[must_use]
    pub fn with_umi_whitelist(method: ClusterMethod, whitelist: BTreeSet<Umi>) -> Self {
        utilities::info(&format!("Length of UMI whitelist: {}", whitelist.len()));
        Self {
            clusterer: UmiClusterer::new(method),
            umi_whitelist: Some(whitelist),
            umi_whitelist_counts: BTreeMap::new(),
        }
    }

    /// The underlying clusterer, whose statistics `dedup` reports at the end
    /// of a run (Python source: `ReadDeduplicator.UMIClusterer`).
    #[must_use]
    pub fn clusterer(&self) -> &UmiClusterer {
        &self.clusterer
    }

    /// Per-whitelist-UMI read counts, plus the [`NON_WHITELIST_UMI`] total
    /// (Python source: `ReadDeduplicator.umi_whitelist_counts`).
    #[must_use]
    pub fn umi_whitelist_counts(&self) -> &BTreeMap<String, u64> {
        &self.umi_whitelist_counts
    }

    /// Deduplicates one bundle (Python source: `ReadDeduplicator.__call__`).
    ///
    /// Without a whitelist, every group is retained and its count is the sum
    /// over its members. With one, a group is retained only if its *leading*
    /// UMI is whitelisted — the whole group is discarded otherwise, on the
    /// reasoning that the group's parent is what the reads are attributed to.
    ///
    /// # Errors
    ///
    /// Propagates [`UmiClusterer::cluster`]'s errors, and returns
    /// [`PortError::Value`] if a whitelisted UMI is not valid UTF-8 — Python
    /// calls `.decode()` there, which raises `UnicodeDecodeError`.
    pub fn dedup(&mut self, counts: &CountMap, threshold: u32) -> Result<DedupResult> {
        let umis: Vec<Umi> = counts.keys().cloned().collect();
        self.dedup_ordered(&umis, counts, threshold)
    }

    /// Deduplicates one bundle, clustering its UMIs in the given order.
    ///
    /// See [`UmiClusterer::cluster_ordered`] for why the order matters.
    ///
    /// # Errors
    ///
    /// As [`dedup`](Self::dedup).
    pub fn dedup_ordered(
        &mut self,
        umis: &[Umi],
        counts: &CountMap,
        threshold: u32,
    ) -> Result<DedupResult> {
        let clusters = self.clusterer.cluster_ordered(umis, counts, threshold)?;

        let sum_of = |cluster: &[Umi]| -> u64 {
            cluster
                .iter()
                .map(|u| counts.get(u).copied().unwrap_or(0))
                .sum()
        };

        let Some(whitelist) = self.umi_whitelist.as_ref() else {
            let mut final_umis = Vec::with_capacity(clusters.len());
            let mut umi_counts = Vec::with_capacity(clusters.len());
            for cluster in &clusters {
                if let Some(first) = cluster.first() {
                    final_umis.push(first.clone());
                    umi_counts.push(sum_of(cluster));
                }
            }
            return Ok(DedupResult {
                final_umis,
                umi_counts,
            });
        };

        let mut final_umis = Vec::new();
        let mut umi_counts = Vec::new();
        for cluster in &clusters {
            let Some(top) = cluster.first() else { continue };
            let cluster_count = sum_of(cluster);
            if whitelist.contains(top) {
                let key = utilities::decode_utf8(top)?;
                final_umis.push(top.clone());
                umi_counts.push(cluster_count);
                *self.umi_whitelist_counts.entry(key).or_insert(0) += cluster_count;
            } else {
                *self
                    .umi_whitelist_counts
                    .entry(NON_WHITELIST_UMI.to_string())
                    .or_insert(0) += cluster_count;
            }
        }
        Ok(DedupResult {
            final_umis,
            umi_counts,
        })
    }
}

/// Clusters cell barcodes, allowing indels (Python source:
/// `network.CellClusterer`).
///
/// Unlike [`UmiClusterer`], the fuzzy path tolerates a single insertion or
/// deletion as well as substitutions, because cell barcodes in droplet
/// protocols vary in length.
///
/// # This type's `call` method does not work, and cannot
///
/// `CellClusterer.__call__` is dead code in the Python original, carrying
/// three independent defects: it passes four arguments to a three-parameter
/// method, it then calls `self.get_groups`, which `__init__` never assigns,
/// and `_get_connected_components_adjacency` calls `.sort()` on the `set` that
/// `breadth_first_search` returns. Running it raises
/// `TypeError: CellClusterer._get_connected_components_adjacency() takes 3
/// positional arguments but 4 were given` — confirmed against the live
/// package, not inferred. Nothing in `umi_tools` calls it.
///
/// [`CellClusterer::call`] therefore returns
/// [`PortError::UpstreamDefect`] rather than inventing an algorithm the
/// original never ran. The methods that *do* work
/// ([`adj_list_directional`](Self::adj_list_directional),
/// [`connected_components_adjacency`](Self::connected_components_adjacency),
/// [`get_best_min_account`](Self::get_best_min_account)) are ported and
/// usable. See `05_design.md` deviation D3.
#[derive(Debug, Clone)]
pub struct CellClusterer {
    dir_threshold: f64,
    fuzzy_match: bool,
}

impl CellClusterer {
    /// Creates a cell-barcode clusterer.
    ///
    /// `dir_threshold` is the count ratio above which one barcode absorbs
    /// another in the fuzzy path (Python default: 10). `fuzzy_match` selects
    /// between the indel-tolerant regex path and the plain Hamming path.
    ///
    /// # Errors
    ///
    /// Returns [`PortError::Value`] for any method other than `"directional"`,
    /// matching Python's `raise ValueError("CellClusterer currently only
    /// supports the directional method")`.
    pub fn new(cluster_method: &str, dir_threshold: f64, fuzzy_match: bool) -> Result<Self> {
        if cluster_method != "directional" {
            return Err(PortError::value(
                "CellClusterer currently only supports the directional method",
            ));
        }
        Ok(Self {
            dir_threshold,
            fuzzy_match,
        })
    }

    /// The count ratio used by the fuzzy path.
    #[must_use]
    pub fn dir_threshold(&self) -> f64 {
        self.dir_threshold
    }

    /// Whether indel-tolerant matching is enabled.
    #[must_use]
    pub fn fuzzy_match(&self) -> bool {
        self.fuzzy_match
    }

    /// Builds the directional adjacency list for cell barcodes
    /// (Python source: `CellClusterer._get_adj_list_directional`).
    ///
    /// Two different rules, selected by `fuzzy_match`:
    ///
    /// - **Fuzzy** (Python's default): for every ordered pair, connect
    ///   `a -> b` when `count(a) >= count(b) * dir_threshold`, the lengths
    ///   differ by at most one, and `b` matches `a` within one error — tested
    ///   with two regexes, one allowing any single edit and one allowing a
    ///   single insertion into the *reversed* barcode, which is how the
    ///   original catches deletions that the first pattern misses (its comment
    ///   gives `UMI1 = "ATCG", UMI2 = "ATTC"` as the motivating case).
    /// - **Exact**: the same rule as [`UmiClusterer::adj_list_directional`]
    ///   with a hard-coded threshold of 1 and the `2n - 1` count condition.
    ///   Note that `dir_threshold` is ignored on this path in the original.
    #[must_use]
    pub fn adj_list_directional(&self, umis: &[Umi], counts: &CountMap) -> AdjMap {
        if !self.fuzzy_match {
            return UmiClusterer::adj_list_directional(umis, counts, 1);
        }

        let mut adj: AdjMap = umis.iter().map(|u| (u.clone(), Vec::new())).collect();
        for a in umis {
            let ca = counts.get(a).copied().unwrap_or(0) as f64;
            for b in umis {
                if a == b {
                    continue;
                }
                let cb = counts.get(b).copied().unwrap_or(0) as f64;
                if ca < cb * self.dir_threshold {
                    continue;
                }
                if a.len().abs_diff(b.len()) > 1 {
                    continue;
                }
                if crate::fuzzy::cell_barcode_matches(a, b, 1) {
                    adj.entry(a.clone()).or_default().push(b.clone());
                }
            }
        }
        adj
    }

    /// Connected components of a cell-barcode graph
    /// (Python source: `CellClusterer._get_connected_components_adjacency`).
    ///
    /// The Python original calls `.sort()` on the `set` that
    /// `breadth_first_search` returns, which raises `AttributeError`. This is
    /// the intended traversal with that defect fixed, since the fix is
    /// unambiguous: `UMIClusterer`'s otherwise identical method does
    /// `component = list(breadth_first_search(...))` before sorting.
    #[must_use]
    pub fn connected_components_adjacency(
        umis: &[Umi],
        graph: &AdjMap,
        counts: &CountMap,
    ) -> Vec<Vec<Umi>> {
        UmiClusterer::connected_components_adjacency(umis, graph, counts)
    }

    /// The fewest highest-count barcodes covering `cluster`
    /// (Python source: `CellClusterer._get_best_min_account`, a verbatim copy
    /// of `UMIClusterer`'s).
    #[must_use]
    pub fn get_best_min_account(
        cluster: &[Umi],
        adj: &AdjMap,
        counts: &CountMap,
    ) -> Option<Vec<Umi>> {
        UmiClusterer::get_best_min_account(cluster, adj, counts)
    }

    /// Not usable — see the type-level documentation.
    ///
    /// # Errors
    ///
    /// Always returns [`PortError::UpstreamDefect`]. The Python original
    /// always raises `TypeError` here.
    pub fn call(&self, _umis: &[Umi], _counts: &CountMap) -> Result<Vec<Vec<Umi>>> {
        Err(PortError::UpstreamDefect {
            symbol: "network.CellClusterer.__call__",
            detail: "the Python original always raises TypeError: \
                     CellClusterer._get_connected_components_adjacency() takes 3 \
                     positional arguments but 4 were given; it also calls \
                     self.get_groups, which __init__ never assigns. Use \
                     adj_list_directional and connected_components_adjacency \
                     directly.",
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts_of(pairs: &[(&[u8], u64)]) -> CountMap {
        pairs.iter().map(|(u, c)| (u.to_vec(), *c)).collect()
    }

    fn umi(s: &str) -> Umi {
        s.as_bytes().to_vec()
    }

    #[test]
    fn method_names_round_trip() {
        for m in ClusterMethod::all() {
            assert_eq!(m.as_str().parse::<ClusterMethod>().unwrap(), m);
            assert_eq!(m.to_string(), m.as_str());
        }
    }

    #[test]
    fn unknown_method_is_rejected_at_construction() {
        let err = UmiClusterer::from_name("bogus").unwrap_err();
        assert!(matches!(err, PortError::UnknownClusterMethod(ref s) if s == "bogus"));
    }

    #[test]
    fn substr_slices_match_python_divmod_split() {
        assert_eq!(get_substr_slices(6, 2), vec![(0, 3), (3, 6)]);
        assert_eq!(get_substr_slices(5, 2), vec![(0, 3), (3, 5)]);
        assert_eq!(get_substr_slices(7, 3), vec![(0, 3), (3, 5), (5, 7)]);
        assert_eq!(get_substr_slices(1, 1), vec![(0, 1)]);
        assert_eq!(get_substr_slices(1, 3), vec![(0, 1), (1, 1), (1, 1)]);
        assert!(get_substr_slices(4, 0).is_empty());
    }

    #[test]
    fn bfs_finds_whole_component_and_stops_at_boundary() {
        let mut adj: AdjMap = BTreeMap::new();
        adj.insert(umi("A"), vec![umi("B"), umi("C")]);
        adj.insert(umi("B"), vec![umi("A")]);
        adj.insert(umi("C"), vec![umi("A"), umi("D")]);
        adj.insert(umi("D"), vec![umi("C")]);
        adj.insert(umi("E"), vec![]);

        let component = breadth_first_search(b"A", &adj);
        assert_eq!(component.len(), 4);
        assert!(!component.contains(&umi("E")));
        assert_eq!(breadth_first_search(b"E", &adj).len(), 1);
    }

    #[test]
    fn bfs_on_unknown_node_returns_just_that_node() {
        let adj: AdjMap = BTreeMap::new();
        assert_eq!(breadth_first_search(b"Z", &adj).len(), 1);
    }

    #[test]
    fn recursive_search_agrees_with_bfs() {
        let mut adj: AdjMap = BTreeMap::new();
        adj.insert(umi("A"), vec![umi("B"), umi("C")]);
        adj.insert(umi("B"), vec![umi("A")]);
        adj.insert(umi("C"), vec![umi("A"), umi("D")]);
        adj.insert(umi("D"), vec![umi("C")]);
        assert_eq!(
            recursive_search(b"A", &adj).unwrap(),
            breadth_first_search(b"A", &adj)
        );
        assert_eq!(
            breadth_first_search_recursive(b"A", &adj),
            breadth_first_search(b"A", &adj)
        );
    }

    #[test]
    fn recursive_search_falls_back_beyond_the_depth_limit() {
        // A path longer than RECURSION_LIMIT forces the fallback; the BFS
        // result must still be the whole path.
        let n = RECURSION_LIMIT + 5;
        let mut adj: AdjMap = BTreeMap::new();
        for i in 0..n {
            let mut neighbours = Vec::new();
            if i + 1 < n {
                neighbours.push(format!("{:06}", i + 1).into_bytes());
            }
            adj.insert(format!("{i:06}").into_bytes(), neighbours);
        }
        assert!(recursive_search(b"000000", &adj).is_none());
        assert_eq!(breadth_first_search_recursive(b"000000", &adj).len(), n);
    }

    #[test]
    fn remove_umis_drops_nodes_and_their_neighbours() {
        let mut adj: AdjMap = BTreeMap::new();
        adj.insert(umi("A"), vec![umi("B")]);
        adj.insert(umi("B"), vec![umi("A")]);
        adj.insert(umi("C"), vec![]);
        let cluster = vec![umi("A"), umi("B"), umi("C")];
        assert_eq!(
            remove_umis(&adj, &cluster, &[umi("A")]),
            [umi("C")].into_iter().collect::<BTreeSet<_>>()
        );
        assert!(remove_umis(&adj, &cluster, &[umi("A"), umi("C")]).is_empty());
    }

    #[test]
    fn directional_absorbs_only_when_counts_justify_it() {
        let counts = counts_of(&[(b"AAAA", 100), (b"AAAT", 2)]);
        let adj = UmiClusterer::adj_list_directional(&[umi("AAAA"), umi("AAAT")], &counts, 1);
        assert_eq!(adj[&umi("AAAA")], vec![umi("AAAT")]);
        assert!(adj[&umi("AAAT")].is_empty());

        // Near-equal counts connect both ways: 5 >= 2*5 - 1 is false, but
        // 5 >= 2*3 - 1 is false too, while 5 >= 2*3-1=5 holds.
        let counts = counts_of(&[(b"AAAA", 5), (b"AAAT", 3)]);
        let adj = UmiClusterer::adj_list_directional(&[umi("AAAA"), umi("AAAT")], &counts, 1);
        assert_eq!(adj[&umi("AAAA")], vec![umi("AAAT")]);
        assert!(adj[&umi("AAAT")].is_empty());
    }

    #[test]
    fn equal_counts_connect_in_both_directions() {
        let counts = counts_of(&[(b"AAAA", 1), (b"AAAT", 1)]);
        let adj = UmiClusterer::adj_list_directional(&[umi("AAAA"), umi("AAAT")], &counts, 1);
        assert_eq!(adj[&umi("AAAA")], vec![umi("AAAT")]);
        assert_eq!(adj[&umi("AAAT")], vec![umi("AAAA")]);
    }

    #[test]
    fn adjacency_edges_are_symmetric() {
        let adj = UmiClusterer::adj_list_adjacency(&[umi("AAAA"), umi("AAAT"), umi("CCCC")], 1);
        assert_eq!(adj[&umi("AAAA")], vec![umi("AAAT")]);
        assert_eq!(adj[&umi("AAAT")], vec![umi("AAAA")]);
        assert!(adj[&umi("CCCC")].is_empty());
    }

    #[test]
    fn threshold_zero_connects_nothing_distinct() {
        let adj = UmiClusterer::adj_list_adjacency(&[umi("AAAA"), umi("AAAT")], 0);
        assert!(adj[&umi("AAAA")].is_empty());
        assert!(adj[&umi("AAAT")].is_empty());
    }

    #[test]
    fn single_umi_clusters_to_itself_under_every_method() {
        let counts = counts_of(&[(b"ACGT", 7)]);
        for method in ClusterMethod::all() {
            let mut c = UmiClusterer::new(method);
            assert_eq!(
                c.cluster(&counts, 1).unwrap(),
                vec![vec![umi("ACGT")]],
                "method {method}"
            );
        }
    }

    #[test]
    fn empty_input_is_a_value_error() {
        let mut c = UmiClusterer::new(ClusterMethod::Directional);
        let err = c.cluster(&CountMap::new(), 1).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m == "max() arg is an empty sequence"));
    }

    #[test]
    fn unequal_lengths_are_an_assertion_error() {
        let counts = counts_of(&[(b"AAA", 1), (b"AAAA", 2)]);
        let mut c = UmiClusterer::new(ClusterMethod::Directional);
        let err = c.cluster(&counts, 1).unwrap_err();
        assert!(matches!(err, PortError::Assertion(ref m)
            if m == "not all umis are the same length(!):  3 - 4"));
    }

    #[test]
    fn statistics_accumulate_across_positions() {
        let mut c = UmiClusterer::new(ClusterMethod::Unique);
        c.cluster(&counts_of(&[(b"AA", 1), (b"CC", 1)]), 1).unwrap();
        c.cluster(&counts_of(&[(b"AA", 1)]), 1).unwrap();
        assert_eq!(c.positions(), 2);
        assert_eq!(c.total_umis_per_position(), 3);
        assert_eq!(c.max_umis_per_position(), 2);
    }

    #[test]
    fn statistics_update_even_when_clustering_fails() {
        // Python increments the counters before its length assertion, so a
        // failing position still counts.
        let mut c = UmiClusterer::new(ClusterMethod::Directional);
        assert!(c
            .cluster(&counts_of(&[(b"AAA", 1), (b"AAAA", 1)]), 1)
            .is_err());
        assert_eq!(c.positions(), 1);
        assert_eq!(c.total_umis_per_position(), 2);
    }

    #[test]
    fn percentile_keeps_only_counts_above_one_percent_of_the_median() {
        // median([500, 1, 3]) == 3, threshold == 0.03, so all three survive.
        let counts = counts_of(&[(b"ACGT", 500), (b"AAAA", 1), (b"CCCC", 3)]);
        let mut c = UmiClusterer::new(ClusterMethod::Percentile);
        assert_eq!(c.cluster(&counts, 1).unwrap().len(), 3);

        // With a large median, a count of 1 falls below the 1% threshold.
        let counts = counts_of(&[(b"ACGT", 1000), (b"AAAA", 1), (b"CCCC", 900)]);
        let mut c = UmiClusterer::new(ClusterMethod::Percentile);
        let groups = c.cluster(&counts, 1).unwrap();
        assert_eq!(groups.len(), 2);
        assert!(!groups.iter().any(|g| g[0] == umi("AAAA")));
    }

    #[test]
    fn substring_index_and_pairwise_paths_agree_at_the_boundary() {
        // 26 UMIs takes the substring-index path, 25 the pairwise one. Over a
        // UMI set where every pair is comparable, both must find the same graph.
        let mut umis: Vec<Umi> = Vec::new();
        for a in b"ACGT" {
            for b in b"ACGT" {
                for c in b"ACGT" {
                    umis.push(vec![*a, *b, *c, b'A']);
                }
            }
        }
        umis.sort();
        assert!(umis.len() > SUBSTR_INDEX_THRESHOLD);

        let indexed = UmiClusterer::adj_list_adjacency(&umis, 1);
        let pairwise = {
            let mut adj: AdjMap = umis.iter().map(|u| (u.clone(), Vec::new())).collect();
            for i in 0..umis.len() {
                for j in (i + 1)..umis.len() {
                    if within_threshold(&umis[i], &umis[j], 1) {
                        adj.get_mut(&umis[i]).unwrap().push(umis[j].clone());
                        adj.get_mut(&umis[j]).unwrap().push(umis[i].clone());
                    }
                }
            }
            adj
        };

        for (k, v) in &indexed {
            let mut a: Vec<&Umi> = v.iter().collect();
            let mut b: Vec<&Umi> = pairwise[k].iter().collect();
            a.sort();
            b.sort();
            assert_eq!(a, b, "neighbours of {k:?}");
        }
    }

    #[test]
    fn dedup_without_whitelist_sums_group_counts() {
        let counts = counts_of(&[(b"AAAA", 100), (b"AAAT", 2), (b"CCCC", 7)]);
        let mut d = ReadDeduplicator::new(ClusterMethod::Directional);
        let out = d.dedup(&counts, 1).unwrap();
        assert_eq!(out.final_umis, vec![umi("AAAA"), umi("CCCC")]);
        assert_eq!(out.umi_counts, vec![102, 7]);
        assert_eq!(d.clusterer().positions(), 1);
    }

    #[test]
    fn dedup_with_whitelist_discards_whole_off_whitelist_groups() {
        let counts = counts_of(&[(b"AAAA", 100), (b"AAAT", 2), (b"CCCC", 7)]);
        let whitelist: BTreeSet<Umi> = [umi("AAAA")].into_iter().collect();
        let mut d = ReadDeduplicator::with_umi_whitelist(ClusterMethod::Directional, whitelist);
        let out = d.dedup(&counts, 1).unwrap();
        assert_eq!(out.final_umis, vec![umi("AAAA")]);
        assert_eq!(out.umi_counts, vec![102]);
        assert_eq!(d.umi_whitelist_counts()[NON_WHITELIST_UMI], 7);
        assert_eq!(d.umi_whitelist_counts()["AAAA"], 102);
    }

    #[test]
    fn cell_clusterer_rejects_non_directional_methods() {
        let err = CellClusterer::new("cluster", 10.0, true).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.contains("only supports the directional method")));
    }

    #[test]
    fn cell_clusterer_call_reports_the_upstream_defect() {
        let c = CellClusterer::new("directional", 10.0, false).unwrap();
        let err = c.call(&[], &CountMap::new()).unwrap_err();
        assert!(matches!(err, PortError::UpstreamDefect { symbol, .. }
            if symbol == "network.CellClusterer.__call__"));
    }

    #[test]
    fn cell_clusterer_exact_path_matches_umi_directional() {
        let counts = counts_of(&[(b"ACGT", 100), (b"ACGA", 2), (b"TTTT", 50)]);
        let umis: Vec<Umi> = counts.keys().cloned().collect();
        let c = CellClusterer::new("directional", 10.0, false).unwrap();
        assert_eq!(
            c.adj_list_directional(&umis, &counts),
            UmiClusterer::adj_list_directional(&umis, &counts, 1)
        );
    }
}
