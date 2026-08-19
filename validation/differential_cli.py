"""Differential test: the `umi_tools` CLI, Python vs Rust, case by case.

Runs every case in the package's own `tests/tests.yaml` twice — once against
the original Python entry point, once against the Rust binary — and compares
the outputs the case declares. This is the end-to-end evidence: it exercises
the same code paths a user does, over the package's own fixtures, and compares
the same files the upstream suite compares.

It is deliberately *not* a copy of `tests/test_umi_tools.py`. That harness
compares against checked-in reference files; this one compares the two
implementations against **each other**, on the same machine in the same
session, so a difference cannot be blamed on a stale reference.

Both binaries are located explicitly rather than through `PATH`, so there is no
question which one ran.

Usage:

    python validation/differential_cli.py [--python-bin PATH] [--rust-bin PATH]
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _recorder import ResultRecorder  # noqa: E402

import yaml  # noqa: E402

SOURCE = Path(os.environ.get("UMI_TOOLS_SRC", "../umi_tools_src"))
TESTS = SOURCE / "tests"


def load_cases():
    with open(TESTS / "tests.yaml") as fh:
        return yaml.safe_load(fh)


def run_case(binary, options, stdin, workdir):
    """Run one case, returning (returncode, {output name: bytes})."""
    tmpdir = Path(tempfile.mkdtemp())
    options = re.sub(r"%TMP%|<TMP>", str(tmpdir), options)
    options = re.sub(r"%DIR%|<DIR>", str(TESTS), options)
    options = options.replace("\n", "")

    stdin_arg = "--stdin=%s" % (TESTS / stdin) if stdin else ""
    stdout_path = tmpdir / "stdout"
    statement = "%s %s %s > %s" % (binary, options, stdin_arg, stdout_path)

    env = os.environ.copy()
    env["COLUMNS"] = "80"
    proc = subprocess.run(["/bin/bash", "-c", statement], cwd=tmpdir, env=env,
                          capture_output=True, check=False)
    return proc.returncode, tmpdir, stdout_path


def read_output(path):
    """Read one output, decompressing `.gz` as the upstream harness does.

    A gzip stream embeds a modification time and an OS byte in its header, so
    two runs producing identical content still differ byte for byte. The
    package's own harness gunzips before comparing; so does this one.
    """
    if not path.exists():
        return None
    if path.suffix == ".gz":
        import gzip
        with gzip.open(path, "rb") as fh:
            return fh.read()
    return path.read_bytes()


def collect(outputs, tmpdir, stdout_path):
    """Read each declared output, keyed by its declared name."""
    collected = {}
    for output in outputs:
        if output == "stdout":
            path = stdout_path
        elif output.startswith("<DIR>/") or output.startswith("%DIR%/"):
            path = TESTS / output[6:]
        else:
            path = tmpdir / output
        collected[output] = read_output(path)
    return collected


def normalise(data, sort):
    """Prepare an output for comparison, as the upstream harness does.

    Two adjustments, both taken from `tests/test_umi_tools.py`:

    - `#`-prefixed and timestamped log lines are dropped. Both implementations
      stamp their own runtime, host and process id into the header, footer and
      log lines, so those differ by construction.
    - When the case declares `sort: True`, lines are sorted before comparing.
      The upstream harness does this for cases whose output order is not
      itself meaningful — `group --method=adjacency`, whose trailing group
      members come from a Python `set` difference, is the example.
    """
    if data is None:
        return None
    lines = []
    for line in data.split(b"\n"):
        if line.startswith(b"#") or line.startswith(b"@PG"):
            continue
        # Log lines from either implementation.
        if re.match(rb"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", line):
            continue
        lines.append(line)
    if sort:
        lines.sort()
    return b"\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python-bin",
                        default=os.environ.get("UMI_TOOLS_PY", "umi_tools"))
    parser.add_argument("--rust-bin",
                        default=str(Path(__file__).resolve().parents[1]
                                    / "target/release/umi_tools"))
    args = parser.parse_args()

    rec = ResultRecorder("cli")
    rec.log("setup", "python binary: %s" % args.python_bin)
    rec.log("setup", "rust binary:   %s" % args.rust_bin)

    # Cases blocked by a documented capability gap rather than a behaviour
    # difference. Listing them here keeps them visible in the results rather
    # than folded into a pass count — see `05_design.md` deviations D12/D13.
    KNOWN_GAPS = {
        "dedup_from_cram": "D12: reference given as a URL; htslib fetches it, this port reads local files",
        "dedup_to_cram": "D12 and D13: URL reference, and CRAM output of sequence-less records",
        "dedup_bam_to_cram": "D12 and D13: URL reference, and CRAM output of sequence-less records",
    }

    cases = load_cases()
    identical = differing = errored = blocked = 0
    differing_names = []
    blocked_names = []

    for name in sorted(cases):
        values = cases[name]
        options = values.get("options", "")
        stdin = values.get("stdin")
        outputs = values.get("outputs", ["stdout"])
        sort = bool(values.get("sort", False))

        py_code, py_tmp, py_stdout = run_case(args.python_bin, options, stdin, TESTS)
        rs_code, rs_tmp, rs_stdout = run_case(args.rust_bin, options, stdin, TESTS)

        if py_code != rs_code:
            if name in KNOWN_GAPS:
                blocked += 1
                blocked_names.append("%s — %s" % (name, KNOWN_GAPS[name]))
            else:
                errored += 1
                differing_names.append("%s (exit %d vs %d)"
                                       % (name, py_code, rs_code))
        else:
            py_files = collect(outputs, py_tmp, py_stdout)
            rs_files = collect(outputs, rs_tmp, rs_stdout)
            same = all(normalise(py_files[o], sort) == normalise(rs_files[o], sort)
                       for o in outputs)
            if same:
                identical += 1
            else:
                differing += 1
                differing_names.append(name)

        shutil.rmtree(py_tmp, ignore_errors=True)
        shutil.rmtree(rs_tmp, ignore_errors=True)

    total = identical + differing + errored + blocked
    comparable = total - blocked

    rec.record("cli_cases", "identical_outputs", value=identical,
               threshold=comparable, tier=1, passed=identical == comparable)
    rec.record("cli_cases", "differing_outputs", value=differing, threshold=0,
               tier=1, passed=differing == 0)
    rec.record("cli_cases", "unexpected_exit_code_differences", value=errored,
               threshold=0, tier=1, passed=errored == 0)
    # Reported so the gaps stay visible; passing means the count matches what
    # is documented, not that the gaps do not exist.
    rec.record("cli_cases", "blocked_by_documented_gap", value=blocked,
               threshold=len(KNOWN_GAPS), tier=1,
               passed=blocked == len(KNOWN_GAPS))

    rec.log("cli_cases", "%d/%d comparable cases produced identical output"
            % (identical, comparable))
    rec.log("cli_cases", "%d case(s) blocked by a documented gap" % blocked)
    for name in blocked_names:
        rec.log("cli_cases", "blocked: %s" % name)
    for name in differing_names:
        rec.log("cli_cases", "differs: %s" % name)
    rec.save()
    return 0


if __name__ == "__main__":
    sys.exit(main())
