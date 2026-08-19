//! The top-level `umi_tools --help` text
//! (Python source: `umi_tools/umi_tools.py`'s module docstring, which `main`
//! prints verbatim).
//!
//! Reproduced from the live module's `__doc__` rather than retyped, so it is
//! exactly what the original prints — including the author and tag lines, and
//! the "6 tools" count that the list of seven immediately contradicts. The
//! docstring is user-facing output that the package's own test suite compares
//! byte for byte; correcting it would be a change, not a fix.

/// The module docstring `umi_tools --help` prints.
pub const MODULE_DOC: &str = r"
umi_tools.py - Tools for UMI analyses
=====================================

:Author: Tom Smith & Ian Sudbery, CGAT
:Tags: Genomics UMI

There are 6 tools:

  - whitelist
  - extract
  - group
  - dedup
  - count
  - count_tab
  - prepare_for_em

To get help on a specific tool, type:

    umi_tools <tool> --help

To use a specific tool, type::

    umi_tools <tool> [tool options] [tool arguments]
";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_docstring_lists_every_tool() {
        for tool in [
            "whitelist",
            "extract",
            "group",
            "dedup",
            "count",
            "count_tab",
            "prepare_for_em",
        ] {
            assert!(MODULE_DOC.contains(tool), "docstring omits {tool}");
        }
    }

    #[test]
    fn the_docstring_is_reproduced_verbatim() {
        // The upstream text says "6 tools" while listing seven. Reproducing
        // that is the point: this is compared byte for byte.
        assert!(MODULE_DOC.contains("There are 6 tools:"));
        assert!(MODULE_DOC.starts_with('\n'));
    }
}
