"""Export each subcommand's `--help` text, captured from the live CLI.

The package's own test suite compares `umi_tools <tool> --help` byte for byte
against checked-in references, and that text is produced by `optparse`'s
`IndentedHelpFormatter` as customised by `Utilities.BetterFormatter` — its
wrapping, its indent positions, its `%default` expansion. Re-deriving it in
Rust would mean hand-transcribing some sixty option help strings per
subcommand and reimplementing the formatter, with no way to notice a dropped
word.

Capturing it from the live command removes that risk: what lands in Rust is
exactly what the original prints. The text is fixed at export time, which is
right for a port of a specific version — the crate already pins itself to
umi_tools 1.1.6.

`COLUMNS=80` matches what the test harness sets, since `optparse` wraps to the
terminal width.
"""
import importlib
import io
import json
import os
import sys

OUT = os.path.dirname(os.path.abspath(__file__))

SUBCOMMANDS = [
    "extract", "whitelist", "group", "dedup", "count", "count_tab",
    "prepare_for_em",
]


def run_help(name):
    """Run `<name> --help` and capture what it prints before exiting."""
    module = importlib.import_module("umi_tools." + name)
    buf = io.StringIO()
    saved_out = sys.stdout
    sys.stdout = buf
    try:
        module.main([name, "--help"])
    except SystemExit:
        pass
    finally:
        sys.stdout = saved_out
    return buf.getvalue()


def main():
    os.environ["COLUMNS"] = "80"
    payload = {}
    for name in SUBCOMMANDS:
        try:
            payload[name] = run_help(name)
        except Exception as exc:  # noqa: BLE001
            payload[name] = None
            print("could not render %s: %s: %s" % (name, type(exc).__name__, exc),
                  file=sys.stderr)
    path = os.path.join(OUT, "reference_help.json")
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=1)
    ok = sum(1 for v in payload.values() if v)
    print("wrote %s (%d/%d subcommands captured)" % (path, ok, len(SUBCOMMANDS)))


if __name__ == "__main__":
    main()
