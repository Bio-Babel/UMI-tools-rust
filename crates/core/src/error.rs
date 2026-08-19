//! Crate-wide error type, mapping the Python package's exception surface
//! (Python source: raises throughout `umi_tools/*.py`).
//!
//! One variant per distinct Python exception class the ported API can raise,
//! so `crates/pyo3-bindings` can reconstruct the exact Python exception type
//! rather than collapsing everything into `RuntimeError`. The mapping is
//! tabulated in `port_reports/umi_tools/05_design.md` § Error model.

use std::path::PathBuf;

use thiserror::Error;

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, PortError>;

/// Every recoverable failure the ported API surface can produce.
///
/// Panics are reserved for violated invariants internal to this crate; no
/// user input reaches one.
#[derive(Debug, Error)]
pub enum PortError {
    /// Equivalent of Python's `ValueError`.
    #[error("{0}")]
    Value(String),

    /// Equivalent of Python's `AssertionError` (an `assert` in the original).
    #[error("{0}")]
    Assertion(String),

    /// Equivalent of Python's `NotImplementedError`.
    #[error("{0}")]
    NotImplemented(String),

    /// Equivalent of Python's `OSError`/`IOError`.
    #[error("{path}: {source}")]
    Io {
        /// Path that could not be opened, read or written.
        path: PathBuf,
        /// Underlying operating-system error.
        source: std::io::Error,
    },

    /// Malformed input in a FASTQ, SAM or TSV stream.
    ///
    /// The Python original calls `Utilities.error()` for these, which logs at
    /// CRITICAL and calls `sys.exit(1)` from inside a generator. A library
    /// must not terminate the process, so this is a recoverable error and the
    /// CLI reproduces the original's message and exit code — deviation D6 in
    /// `05_design.md`.
    #[error("{0}")]
    Format(String),

    /// A cluster-method name that does not name one of the five methods.
    ///
    /// Python binds no dispatch attributes for an unknown name and fails later
    /// with `AttributeError` on first use; this crate rejects it in the
    /// constructor — deviation D2 in `05_design.md`.
    #[error("unknown cluster method: {0}")]
    UnknownClusterMethod(String),

    /// A code path that is unreachable in the Python original because of a
    /// defect there, and therefore has no behaviour to port.
    ///
    /// Deviation D3 in `05_design.md`.
    #[error("{symbol} is not usable: {detail}")]
    UpstreamDefect {
        /// Fully-qualified Python symbol, e.g. `network.CellClusterer.__call__`.
        symbol: &'static str,
        /// What goes wrong in the Python original.
        detail: &'static str,
    },

    /// A SAM/BAM/CRAM read, write or sort failure.
    #[error("alignment error: {0}")]
    Alignment(String),
}

impl PortError {
    /// Builds a [`PortError::Value`] from anything displayable.
    pub fn value(msg: impl Into<String>) -> Self {
        PortError::Value(msg.into())
    }

    /// Builds a [`PortError::Assertion`] from anything displayable.
    pub fn assertion(msg: impl Into<String>) -> Self {
        PortError::Assertion(msg.into())
    }

    /// Builds a [`PortError::Format`] from anything displayable.
    pub fn format(msg: impl Into<String>) -> Self {
        PortError::Format(msg.into())
    }

    /// Builds a [`PortError::Alignment`] from anything displayable.
    pub fn alignment(msg: impl Into<String>) -> Self {
        PortError::Alignment(msg.into())
    }

    /// Builds a [`PortError::Io`] tagged with the path that failed.
    pub fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        PortError::Io {
            path: path.into(),
            source,
        }
    }
}
