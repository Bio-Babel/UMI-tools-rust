//! Python bindings exposing `umi_tools_core` as the `umi_tools_rs` package.
//!
//! Deliberately thin: it marshals arguments and errors at the boundary and
//! delegates every decision to `umi_tools_core`. Deleting this crate would
//! remove the Python API and change nothing else.
//!
//! # Why `umi_tools_rs` and not `umi_tools._native`
//!
//! The conventional layout installs the extension as a submodule of the
//! package being ported. That does not work here: the original `umi_tools` is
//! `pip install -e .`-ed in the same environment, and `maturin develop` of a
//! mixed-layout `umi_tools` package silently shadows it — `import umi_tools`
//! then resolves to the shim and the pure-Python modules the differential
//! harness needs are gone.
//!
//! A distinct top-level name keeps both importable in one interpreter, which
//! is what same-process differential testing actually requires. See
//! `port_reports/umi_tools/05_design.md` deviation D9.
//!
//! # What is exposed
//!
//! The functions the differential harness compares, plus the clustering API a
//! Python user would plausibly want to call. Anything typed in terms of
//! `noodles` records is deliberately absent: marshalling those into `pysam`
//! objects would cost more than it buys, and the alignment path is Rust-only.

// The bindings marshal Python values by value at the boundary.
#![allow(clippy::needless_pass_by_value)]

use std::collections::{BTreeMap, BTreeSet};

use pyo3::exceptions::{
    PyAssertionError, PyAttributeError, PyIOError, PyNotImplementedError, PyRuntimeError,
    PyTypeError, PyValueError,
};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyList};

use umi_tools_core::extract::{ExtractConfig, ExtractFilterAndUpdate, ExtractMethod};
use umi_tools_core::network::{CountMap, UmiClusterer};
use umi_tools_core::umi_methods::{QualityEncoding, Record};
use umi_tools_core::whitelist::{CellCounts, ErrorCorrectionMode, KneeMethod, ResolutionMethod};
use umi_tools_core::PortError;

/// Reconstructs the Python exception type a caller of the *original* package
/// would have caught for the equivalent failure.
///
/// Never falls through to a generic exception: preserving the exception
/// *type*, not just its message, is what lets an existing `except ValueError:`
/// keep working against the port. The mapping is tabulated in
/// `port_reports/umi_tools/05_design.md` § Error model.
///
/// This is a free function rather than `impl From<PortError> for PyErr`
/// because, with the core and bindings in separate crates, both types are
/// foreign here and that impl would violate the orphan rule.
fn to_py_err(err: PortError) -> PyErr {
    let message = err.to_string();
    match err {
        PortError::Value(_) | PortError::Format(_) => PyValueError::new_err(message),
        PortError::Assertion(_) => PyAssertionError::new_err(message),
        PortError::NotImplemented(_) => PyNotImplementedError::new_err(message),
        PortError::Io { .. } => PyIOError::new_err(message),
        PortError::UnknownClusterMethod(_) => PyAttributeError::new_err(message),
        PortError::UpstreamDefect { .. } => PyTypeError::new_err(message),
        PortError::Alignment(_) => PyRuntimeError::new_err(message),
    }
}

/// Hamming distance between two barcodes
/// (Python source: `umi_tools._dedup_umi.edit_distance`).
///
/// Returns `None` for unequal-length inputs, where the Cython original's
/// `.pyx` source returns `numpy.inf`. The extension as built against NumPy 2
/// raises `AttributeError` instead, because its checked-in generated `.c`
/// still references the removed `np.Inf`; this implements the `.pyx` intent.
#[pyfunction]
fn edit_distance(a: &[u8], b: &[u8]) -> Option<u32> {
    umi_tools_core::edit_distance(a, b)
}

/// Whether two barcodes are within `threshold` substitutions.
#[pyfunction]
fn within_threshold(a: &[u8], b: &[u8], threshold: u32) -> bool {
    umi_tools_core::within_threshold(a, b, threshold)
}

/// The mean pairwise Hamming distance between a set of UMIs
/// (Python source: `umi_methods.get_average_umi_distance`).
///
/// Returns `-1.0` for fewer than two UMIs, as the original does.
#[pyfunction]
fn get_average_umi_distance(umis: Vec<Vec<u8>>) -> f64 {
    umi_tools_core::umi_methods::get_average_umi_distance(&umis)
}

fn count_map_from(counts: &Bound<'_, PyDict>) -> PyResult<CountMap> {
    let mut out = CountMap::new();
    for (key, value) in counts.iter() {
        // The Python API accepts both `bytes` (what the BAM path produces)
        // and `str` (what the flat-file path produces), as the original does.
        let umi: Vec<u8> = match key.extract::<Vec<u8>>() {
            Ok(b) => b,
            Err(_) => key.extract::<String>()?.into_bytes(),
        };
        out.insert(umi, value.extract::<u64>()?);
    }
    Ok(out)
}

fn groups_to_py(py: Python<'_>, groups: &[Vec<Vec<u8>>]) -> PyResult<Py<PyList>> {
    let outer = PyList::empty(py);
    for group in groups {
        let inner = PyList::empty(py);
        for umi in group {
            inner.append(PyBytes::new(py, umi))?;
        }
        outer.append(inner)?;
    }
    Ok(outer.unbind())
}

/// Clusters the UMIs observed at one position
/// (Python source: `network.UMIClusterer`).
///
/// The counterpart of Python's functor: construct with a method name, then
/// call it once per position. Statistics accumulate across calls, as they do
/// upstream.
///
/// # Examples
///
/// ```python
/// import umi_tools_rs
/// clusterer = umi_tools_rs.UMIClusterer("directional")
/// clusterer.cluster({b"ACGT": 10, b"ACGA": 2}, 1)
/// # [[b'ACGT', b'ACGA']]
/// ```
#[pyclass(name = "UMIClusterer")]
struct PyUmiClusterer {
    inner: UmiClusterer,
}

#[pymethods]
impl PyUmiClusterer {
    /// Creates a clusterer for the named method.
    ///
    /// Raises `AttributeError` for an unknown name. The Python original
    /// accepts any string and fails on first use with the same exception
    /// type; raising at construction is more informative and cannot change a
    /// working program's behaviour.
    #[new]
    #[pyo3(signature = (cluster_method = "directional"))]
    fn new(cluster_method: &str) -> PyResult<Self> {
        Ok(Self {
            inner: UmiClusterer::from_name(cluster_method).map_err(to_py_err)?,
        })
    }

    /// Clusters one position's UMIs, returning a list of groups.
    ///
    /// `umis` maps a UMI to how often it was seen; the parent UMI leads each
    /// returned group.
    #[pyo3(signature = (umis, threshold = 1))]
    fn cluster(
        &mut self,
        py: Python<'_>,
        umis: &Bound<'_, PyDict>,
        threshold: u32,
    ) -> PyResult<Py<PyList>> {
        let counts = count_map_from(umis)?;
        let groups = self.inner.cluster(&counts, threshold).map_err(to_py_err)?;
        groups_to_py(py, &groups)
    }

    /// Alias for [`cluster`](Self::cluster), matching the Python original's
    /// functor call.
    #[pyo3(signature = (umis, threshold = 1))]
    fn __call__(
        &mut self,
        py: Python<'_>,
        umis: &Bound<'_, PyDict>,
        threshold: u32,
    ) -> PyResult<Py<PyList>> {
        self.cluster(py, umis, threshold)
    }

    /// The method this clusterer applies.
    #[getter]
    fn method(&self) -> &'static str {
        self.inner.method().as_str()
    }

    /// How many positions have been clustered so far.
    #[getter]
    fn positions(&self) -> u64 {
        self.inner.positions()
    }

    /// Total distinct UMIs summed over every position clustered so far.
    #[getter]
    fn total_umis_per_position(&self) -> u64 {
        self.inner.total_umis_per_position()
    }

    /// The largest number of distinct UMIs seen at any one position.
    #[getter]
    fn max_umis_per_position(&self) -> usize {
        self.inner.max_umis_per_position()
    }
}

/// The adjacency list a clustering method builds
/// (Python source: `UMIClusterer.get_adj_list`).
///
/// Returned as a dict of UMI to neighbour list, with each neighbour list
/// sorted so the result is order-independent; the Python original appends in
/// pair-iteration order, which is an artefact.
#[pyfunction]
#[pyo3(signature = (umis, threshold = 1, method = "directional"))]
fn get_adj_list(
    py: Python<'_>,
    umis: &Bound<'_, PyDict>,
    threshold: u32,
    method: &str,
) -> PyResult<Py<PyDict>> {
    let counts = count_map_from(umis)?;
    let keys: Vec<Vec<u8>> = counts.keys().cloned().collect();
    let clusterer = UmiClusterer::from_name(method).map_err(to_py_err)?;
    let adj = clusterer.get_adj_list(&keys, &counts, threshold);

    let out = PyDict::new(py);
    if let Some(graph) = adj.graph() {
        for (node, neighbours) in graph {
            let mut sorted = neighbours.clone();
            sorted.sort();
            let list = PyList::empty(py);
            for n in sorted {
                list.append(PyBytes::new(py, &n))?;
            }
            out.set_item(PyBytes::new(py, node), list)?;
        }
    }
    Ok(out.unbind())
}

/// Extracts barcodes from a read and rewrites it
/// (Python source: `extract_methods.ExtractFilterAndUpdate`).
///
/// Takes the read's three FASTQ fields and returns the rewritten triple, or
/// `None` when the read was filtered out — the same verdict the original
/// signals by returning `None`.
#[pyclass(name = "ExtractFilterAndUpdate")]
struct PyExtractFilterAndUpdate {
    inner: ExtractFilterAndUpdate,
}

#[pymethods]
impl PyExtractFilterAndUpdate {
    /// Creates an extractor.
    ///
    /// The parameters mirror the Python constructor's, minus the ones its own
    /// body never reads.
    #[new]
    #[pyo3(signature = (
        method = "string",
        pattern = None,
        pattern2 = None,
        prime3 = false,
        extract_cell = false,
        quality_encoding = None,
        quality_filter_threshold = None,
        quality_filter_mask = None,
        retain_umi = false,
        umi_separator = "_",
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        method: &str,
        pattern: Option<String>,
        pattern2: Option<String>,
        prime3: bool,
        extract_cell: bool,
        quality_encoding: Option<String>,
        quality_filter_threshold: Option<i32>,
        quality_filter_mask: Option<i32>,
        retain_umi: bool,
        umi_separator: &str,
    ) -> PyResult<Self> {
        let encoding = quality_encoding
            .map(|e| e.parse::<QualityEncoding>())
            .transpose()
            .map_err(to_py_err)?;
        let config = ExtractConfig {
            method: method.parse::<ExtractMethod>().map_err(to_py_err)?,
            pattern,
            pattern2,
            prime3,
            extract_cell,
            quality_encoding: encoding,
            quality_filter_threshold,
            quality_filter_mask,
            retain_umi,
            umi_separator: umi_separator.to_string(),
            ..ExtractConfig::new()
        };
        Ok(Self {
            inner: ExtractFilterAndUpdate::new(config).map_err(to_py_err)?,
        })
    }

    /// Processes one read, returning `(identifier, seq, quals)` or `None`.
    #[pyo3(signature = (identifier, seq, quals))]
    fn process(
        &mut self,
        identifier: &str,
        seq: &str,
        quals: &str,
    ) -> PyResult<Option<(String, String, String)>> {
        let mut read = Record::new(identifier, seq, quals);
        let kept = self.inner.process(&mut read, None).map_err(to_py_err)?;
        Ok(kept.then_some((read.identifier, read.seq, read.quals)))
    }

    /// Processes a read pair, returning both rewritten reads or `None`.
    #[pyo3(signature = (read1, read2))]
    #[allow(clippy::type_complexity)]
    fn process_pair(
        &mut self,
        read1: (String, String, String),
        read2: (String, String, String),
    ) -> PyResult<Option<((String, String, String), (String, String, String))>> {
        let mut r1 = Record::new(read1.0, read1.1, read1.2);
        let mut r2 = Record::new(read2.0, read2.1, read2.2);
        let kept = self
            .inner
            .process(&mut r1, Some(&mut r2))
            .map_err(to_py_err)?;
        Ok(kept.then_some({
            (
                (r1.identifier, r1.seq, r1.quals),
                (r2.identifier, r2.seq, r2.quals),
            )
        }))
    }

    /// The event tallies accumulated so far
    /// (Python source: `getReadCounts`).
    fn get_read_counts(&self, py: Python<'_>) -> PyResult<Py<PyDict>> {
        let out = PyDict::new(py);
        for (event, count) in self.inner.read_counts() {
            out.set_item(event, count)?;
        }
        Ok(out.unbind())
    }
}

/// Builds a counter from a Python dict, preserving its insertion order — which
/// is what `Counter.most_common` breaks ties on.
fn cell_counts_from(counts: &Bound<'_, PyDict>) -> PyResult<CellCounts> {
    let mut out = CellCounts::new();
    for (key, value) in counts.iter() {
        let barcode: Vec<u8> = match key.extract::<Vec<u8>>() {
            Ok(b) => b,
            Err(_) => key.extract::<String>()?.into_bytes(),
        };
        out.set(&barcode, value.extract::<u64>()?);
    }
    Ok(out)
}

/// Selects cell barcodes by one of the two knee methods
/// (Python source: `whitelist_methods.getCellWhitelist`).
///
/// Returns the selected barcodes as a sorted list of `str`, or `None` when no
/// threshold was accepted — which is what `whitelist` reports as "No local
/// minima was accepted".
#[pyfunction]
#[pyo3(signature = (
    cell_barcode_counts,
    knee_method = "distance",
    expect_cells = None,
    cell_number = None,
))]
fn get_cell_whitelist(
    cell_barcode_counts: &Bound<'_, PyDict>,
    knee_method: &str,
    expect_cells: Option<u64>,
    cell_number: Option<usize>,
) -> PyResult<Option<Vec<String>>> {
    let counts = cell_counts_from(cell_barcode_counts)?;
    let method: KneeMethod = knee_method.parse().map_err(to_py_err)?;
    let result = match method {
        KneeMethod::Distance => {
            umi_tools_core::knee_estimate_distance(&counts, cell_number).map_err(to_py_err)?
        }
        KneeMethod::Density => {
            umi_tools_core::knee_estimate_density(&counts, expect_cells, cell_number)
                .map_err(to_py_err)?
        }
    };
    Ok(result.barcodes.map(|set| {
        set.iter()
            .map(|b| String::from_utf8_lossy(b).into_owned())
            .collect()
    }))
}

/// Maps whitelisted barcodes to the observed barcodes correcting to them
/// (Python source: `whitelist_methods.getErrorCorrectMapping`).
#[pyfunction]
#[pyo3(signature = (cell_barcodes, whitelist, threshold = 1))]
fn get_error_correct_mapping(
    py: Python<'_>,
    cell_barcodes: Vec<String>,
    whitelist: Vec<String>,
    threshold: u32,
) -> PyResult<Py<PyDict>> {
    let barcodes: Vec<Vec<u8>> = cell_barcodes
        .iter()
        .map(|b| b.as_bytes().to_vec())
        .collect();
    let whitelist: BTreeSet<Vec<u8>> = whitelist.iter().map(|b| b.as_bytes().to_vec()).collect();
    let mapping =
        umi_tools_core::whitelist::get_error_correct_mapping(&barcodes, &whitelist, threshold);

    let out = PyDict::new(py);
    for (parent, children) in &mapping {
        let list = PyList::empty(py);
        for child in children {
            list.append(String::from_utf8_lossy(child).into_owned())?;
        }
        out.set_item(String::from_utf8_lossy(parent).into_owned(), list)?;
    }
    Ok(out.unbind())
}

/// Whitelisted barcodes that `barcode` could be an error of
/// (Python source: `whitelist_methods.checkError`).
#[pyfunction]
#[pyo3(signature = (barcode, whitelist, errors = 1))]
fn check_error(barcode: &str, whitelist: Vec<String>, errors: u32) -> Vec<String> {
    let whitelist: Vec<Vec<u8>> = whitelist.iter().map(|b| b.as_bytes().to_vec()).collect();
    umi_tools_core::whitelist::check_error(barcode.as_bytes(), &whitelist, errors)
        .iter()
        .map(|b| String::from_utf8_lossy(b).into_owned())
        .collect()
}

/// Reads a user-supplied whitelist file
/// (Python source: `whitelist_methods.getUserDefinedBarcodes`).
///
/// Returns `(whitelist, correction_map)`; the map is `None` when no error
/// correction was requested, and a value of `None` inside it marks a barcode
/// reachable from more than one whitelist entry.
#[pyfunction]
#[pyo3(signature = (
    whitelist_tsv,
    whitelist_tsv2 = None,
    get_error_correction = false,
    derive_error_correction = false,
    threshold = 1,
))]
fn get_user_defined_barcodes(
    py: Python<'_>,
    whitelist_tsv: &str,
    whitelist_tsv2: Option<&str>,
    get_error_correction: bool,
    derive_error_correction: bool,
    threshold: u32,
) -> PyResult<(Vec<String>, Option<Py<PyDict>>)> {
    // Python prefers `derive` when both flags are set.
    let mode = if derive_error_correction {
        ErrorCorrectionMode::Derive
    } else if get_error_correction {
        ErrorCorrectionMode::FromFile
    } else {
        ErrorCorrectionMode::None
    };
    let (whitelist, mapping) = umi_tools_core::whitelist::get_user_defined_barcodes(
        std::path::Path::new(whitelist_tsv),
        whitelist_tsv2.map(std::path::Path::new),
        mode,
        threshold,
    )
    .map_err(to_py_err)?;

    let names: Vec<String> = whitelist
        .iter()
        .map(|b| String::from_utf8_lossy(b).into_owned())
        .collect();
    let map = mapping
        .map(|m| -> PyResult<Py<PyDict>> {
            let out = PyDict::new(py);
            for (error, parent) in &m {
                out.set_item(
                    String::from_utf8_lossy(error).into_owned(),
                    parent
                        .as_ref()
                        .map(|p| String::from_utf8_lossy(p).into_owned()),
                )?;
            }
            Ok(out.unbind())
        })
        .transpose()?;
    Ok((names, map))
}

/// Removes or corrects whitelisted barcodes that look like errors of others
/// (Python source: `whitelist_methods.errorDetectAboveThreshold`).
///
/// Returns the surviving whitelist, sorted.
#[pyfunction]
#[pyo3(signature = (cell_barcode_counts, cell_whitelist, errors = 1, resolution_method = "discard"))]
fn error_detect_above_threshold(
    cell_barcode_counts: &Bound<'_, PyDict>,
    cell_whitelist: Vec<String>,
    errors: u32,
    resolution_method: &str,
) -> PyResult<Vec<String>> {
    let counts = cell_counts_from(cell_barcode_counts)?;
    // Most-abundant-first, the order the estimator returns and the order the
    // tie-break depends on.
    let selected: BTreeSet<Vec<u8>> = cell_whitelist
        .iter()
        .map(|b| b.as_bytes().to_vec())
        .collect();
    let whitelist: Vec<Vec<u8>> = umi_tools_core::whitelist::most_common(&counts)
        .into_iter()
        .filter(|(b, _)| selected.contains(b))
        .map(|(b, _)| b)
        .collect();
    let method: ResolutionMethod = resolution_method.parse().map_err(to_py_err)?;
    let (remaining, _) = umi_tools_core::whitelist::error_detect_above_threshold(
        &counts, &whitelist, None, errors, method,
    )
    .map_err(to_py_err)?;
    Ok(remaining
        .iter()
        .map(|b| String::from_utf8_lossy(b).into_owned())
        .collect())
}

/// Splits a read's sequence and qualities by base role
/// (Python source: `extract_methods.extractSeqAndQuals`).
///
/// Returns `(new_seq, new_quals, umi_quals, cell_quals)`, matching the
/// original's tuple.
#[pyfunction]
#[pyo3(signature = (seq, quals, umi_bases, cell_bases, discard_bases, retain_umi = false))]
fn extract_seq_and_quals(
    seq: &str,
    quals: &str,
    umi_bases: Vec<usize>,
    cell_bases: Vec<usize>,
    discard_bases: Vec<usize>,
    retain_umi: bool,
) -> (String, String, String, String) {
    let out = umi_tools_core::extract::extract_seq_and_quals(
        seq,
        quals,
        &umi_bases.into_iter().collect(),
        &cell_bases.into_iter().collect(),
        &discard_bases.into_iter().collect(),
        retain_umi,
    );
    (out.new_seq, out.new_quals, out.umi_quals, out.cell_quals)
}

/// Whether any UMI base falls below `threshold`
/// (Python source: `extract_methods.umi_below_threshold`).
#[pyfunction]
fn umi_below_threshold(
    umi_quals: &str,
    quality_encoding: &str,
    quality_filter_threshold: i32,
) -> PyResult<bool> {
    let encoding: QualityEncoding = quality_encoding.parse().map_err(to_py_err)?;
    Ok(umi_tools_core::extract::umi_below_threshold(
        umi_quals,
        encoding,
        quality_filter_threshold,
    ))
}

/// Replaces low-quality UMI bases with `N`
/// (Python source: `extract_methods.mask_umi`).
#[pyfunction]
fn mask_umi(
    umi: &str,
    umi_quals: &str,
    quality_encoding: &str,
    quality_filter_threshold: i32,
) -> PyResult<String> {
    let encoding: QualityEncoding = quality_encoding.parse().map_err(to_py_err)?;
    Ok(umi_tools_core::extract::mask_umi(
        umi,
        umi_quals,
        encoding,
        quality_filter_threshold,
    ))
}

/// Every quality encoding compatible with a quality string
/// (Python source: `umi_methods.Record.guessFormat`).
#[pyfunction]
fn guess_format(quals: &str) -> PyResult<Vec<String>> {
    Record::new("", "", quals)
        .guess_format()
        .map(|v| v.into_iter().map(str::to_string).collect())
        .map_err(to_py_err)
}

/// Aggregates per-position observations into per-UMI summaries
/// (Python source: `dedup.aggregateStatsDF`).
///
/// Returns a dict of UMI to `(median_counts, times_observed, total_counts)`.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn aggregate_stats(py: Python<'_>, observations: Vec<(Vec<u8>, u64)>) -> PyResult<Py<PyDict>> {
    let stats: BTreeMap<_, _> = umi_tools_core::stats::aggregate_stats(&observations);
    let out = PyDict::new(py);
    for (umi, s) in &stats {
        out.set_item(
            PyBytes::new(py, umi),
            (s.median_counts, s.times_observed, s.total_counts),
        )?;
    }
    Ok(out.unbind())
}

/// Matches an extraction pattern against a read, returning the group spans
/// (Python dependency: the `regex` module's `.match()` plus `.span(group)`).
///
/// Returns a dict of group name to `(start, end)`, or `None` when the pattern
/// does not match. This is what `core::extract` slices the read by, so it is
/// the natural unit for a differential comparison.
#[pyfunction]
fn match_pattern(py: Python<'_>, pattern: &str, subject: &str) -> PyResult<Option<Py<PyDict>>> {
    let compiled = umi_tools_core::pattern::Pattern::parse(pattern).map_err(to_py_err)?;
    let Some(hit) = compiled.match_at_start(subject) else {
        return Ok(None);
    };
    let out = PyDict::new(py);
    for (name, span) in &hit.groups {
        out.set_item(name, *span)?;
    }
    Ok(Some(out.unbind()))
}

/// The version of the Python package this crate was ported from.
#[pyfunction]
fn version() -> &'static str {
    umi_tools_core::VERSION
}

/// The `umi_tools_rs._native` extension module.
#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", umi_tools_core::VERSION)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;

    m.add_function(wrap_pyfunction!(edit_distance, m)?)?;
    m.add_function(wrap_pyfunction!(within_threshold, m)?)?;
    m.add_function(wrap_pyfunction!(get_average_umi_distance, m)?)?;
    m.add_function(wrap_pyfunction!(get_adj_list, m)?)?;

    m.add_function(wrap_pyfunction!(extract_seq_and_quals, m)?)?;
    m.add_function(wrap_pyfunction!(umi_below_threshold, m)?)?;
    m.add_function(wrap_pyfunction!(mask_umi, m)?)?;
    m.add_function(wrap_pyfunction!(guess_format, m)?)?;
    m.add_function(wrap_pyfunction!(match_pattern, m)?)?;

    m.add_function(wrap_pyfunction!(get_cell_whitelist, m)?)?;
    m.add_function(wrap_pyfunction!(get_error_correct_mapping, m)?)?;
    m.add_function(wrap_pyfunction!(check_error, m)?)?;
    m.add_function(wrap_pyfunction!(get_user_defined_barcodes, m)?)?;
    m.add_function(wrap_pyfunction!(error_detect_above_threshold, m)?)?;
    m.add_function(wrap_pyfunction!(aggregate_stats, m)?)?;

    m.add_class::<PyUmiClusterer>()?;
    m.add_class::<PyExtractFilterAndUpdate>()?;
    Ok(())
}
