//! The `umi_tools` command-line interface
//! (Python source: `umi_tools/umi_tools.py` and the seven subcommand modules).
//!
//! Everything here is presentation and orchestration: argument parsing, file
//! plumbing, and printing. All the algorithms live in `umi_tools_core`, which
//! has no dependency on this crate and no notion of a command line — a split
//! the Python package cannot make, because every module there imports
//! `Utilities`, which owns the option parser and a pile of process globals.
//!
//! # The seven tools
//!
//! | Command | What it does |
//! |---|---|
//! | [`extract`] | move the UMI from the read sequence into the read name |
//! | [`whitelist`] | pick the real cell barcodes out of an observed distribution |
//! | [`dedup`] | keep one read per original molecule |
//! | [`group`] | tag reads with the group they were assigned to, without removing any |
//! | [`count`] | count molecules per gene from an alignment |
//! | [`count_tab`] | count molecules per gene from a flat file |
//! | [`prepare_for_em`] | make `dedup`/`group` output usable by RSEM and Salmon |
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

pub mod common;
pub mod count;
pub mod count_tab;
pub mod dedup;
pub mod documentation;
pub mod extract;
pub mod group;
pub mod help;
pub mod options;
pub mod prepare_for_em;
pub mod usage;
pub mod whitelist;

use umi_tools_core::{PortError, Result};

/// The documentation URL Python prints alongside the usage text.
pub const DOCS_URL: &str =
    "For full UMI-tools documentation, see: https://umi-tools.readthedocs.io/en/latest/\n";

/// Dispatches to a subcommand (Python source: `umi_tools.umi_tools.main`).
///
/// `argv` is the full argument vector including the program name, as
/// `std::env::args` yields it.
///
/// Returns the process exit code. Python's `main` returns 0 for `--help` and
/// `--version`, 1 for an unknown command, and otherwise delegates — so an
/// exit code rather than a `Result` is the faithful signature.
///
/// # Errors
///
/// Returns whatever the dispatched subcommand returns.
pub fn run(argv: &[String]) -> Result<i32> {
    let Some(command) = argv.get(1) else {
        print_usage();
        return Ok(0);
    };

    if command == "--help" || command == "-h" {
        print_usage();
        return Ok(0);
    }
    if command == "--version" || command == "-v" {
        println!("UMI-tools version: {}", umi_tools_core::VERSION);
        return Ok(0);
    }

    // Python substitutes the old name for the renamed command.
    let command = if command == "prepare-for-rsem" {
        "prepare-for-em"
    } else {
        command.as_str()
    };

    // `del sys.argv[0]` in Python: the subcommand sees its own name first.
    let sub_argv: Vec<String> = argv[1..].to_vec();
    let args: Vec<String> = sub_argv[1..].to_vec();

    match command {
        "extract" => extract::main(&args, &sub_argv),
        "whitelist" => whitelist::main(&args, &sub_argv),
        "group" => group::main(&args, &sub_argv),
        "dedup" => dedup::main(&args, &sub_argv),
        "count" => count::main(&args, &sub_argv),
        "count_tab" => count_tab::main(&args, &sub_argv),
        "prepare-for-em" | "prepare_for_em" => prepare_for_em::main(&args, &sub_argv),
        unknown => {
            println!("'{unknown}' is not a UMI-tools command. See 'umi_tools -h'.\n");
            print_usage();
            Ok(1)
        }
    }
}

/// Whether the argument list asks for help.
///
/// Intercepted before `clap` parses, so `--help` prints the original's own
/// text (see [`help`]) rather than `clap`'s rendering of the same options.
///
/// `--help-extended` is included because `Utilities.Start` declares it with
/// `action='help'`, making it a second spelling of `--help` rather than a
/// different, longer output.
#[must_use]
pub fn wants_help(args: &[String]) -> bool {
    args.iter()
        .any(|a| a == "--help" || a == "-h" || a == "--help-extended")
}

/// Prints the documentation pointer and the module docstring, as Python's
/// `main` does with two `print` calls — each of which appends its own newline,
/// and the docstring itself begins with one.
fn print_usage() {
    println!("{DOCS_URL}");
    println!("{}", usage::MODULE_DOC);
}

/// Reports a fatal error the way `Utilities.error` does, and returns the exit
/// code Python's `sys.exit(1)` produces.
///
/// Python's `error()` logs at ERROR, writes a fixed notice to stderr, and
/// exits 1 — from inside library code. Here the library returns the error and
/// the CLI does all three, so the observable behaviour is unchanged while the
/// library stays usable. See `05_design.md` deviation D6.
pub fn report_fatal(error: &PortError) -> i32 {
    eprint!("{}", umi_tools_core::utilities::error(&error.to_string()));
    1
}

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(args: &[&str]) -> Vec<String> {
        std::iter::once("umi_tools".to_string())
            .chain(args.iter().map(|a| (*a).to_string()))
            .collect()
    }

    #[test]
    fn no_arguments_prints_usage_and_succeeds() {
        assert_eq!(run(&argv(&[])).unwrap(), 0);
    }

    #[test]
    fn help_and_version_succeed() {
        assert_eq!(run(&argv(&["--help"])).unwrap(), 0);
        assert_eq!(run(&argv(&["-h"])).unwrap(), 0);
        assert_eq!(run(&argv(&["--version"])).unwrap(), 0);
        assert_eq!(run(&argv(&["-v"])).unwrap(), 0);
    }

    #[test]
    fn help_extended_is_a_second_spelling_of_help() {
        assert!(wants_help(&["--help-extended".to_string()]));
        assert!(wants_help(&["--help".to_string()]));
        assert!(wants_help(&["-h".to_string()]));
        assert!(!wants_help(&["--helpful".to_string()]));
    }

    #[test]
    fn an_unknown_command_exits_one() {
        assert_eq!(run(&argv(&["nonsense"])).unwrap(), 1);
    }

    #[test]
    fn usage_lists_every_tool() {
        for tool in [
            "whitelist",
            "extract",
            "group",
            "dedup",
            "count",
            "count_tab",
            "prepare_for_em",
        ] {
            assert!(
                usage::MODULE_DOC.contains(tool),
                "usage does not mention {tool}"
            );
        }
    }
}
