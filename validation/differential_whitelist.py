"""Differential test: `umi_tools.whitelist_methods` vs `umi_tools_rs`.

The only genuinely float-sensitive part of the package: both knee estimators
run on floating-point curves, and the density one fits a Gaussian kernel
density estimate whose Python implementation is SciPy's. That makes this a
**Tier 2** comparison for the curves — a different numeric backend for one
step — while the *selected barcode sets* they produce must still match
exactly, because a knee threshold either includes a barcode or it does not.

See `differential_network.py`'s docstring for why the extension is imported as
`umi_tools_rs`.
"""

import collections
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _recorder import ResultRecorder  # noqa: E402

import numpy as np  # noqa: E402
import umi_tools.whitelist_methods as py_wl  # noqa: E402
import umi_tools_rs as rs  # noqa: E402


def distributions():
    """Barcode count distributions spanning the shapes the knee methods see."""
    cases = [
        ("bimodal", {"AAAA": 1000, "AAAT": 800, "AATT": 600, "ATTT": 400,
                     "TTTT": 200, "CCCC": 10, "CCCG": 8, "CCGG": 6,
                     "CGGG": 4, "GGGG": 2}),
    ]
    # Geometric distributions, which is roughly what a droplet run produces.
    for i, seed in enumerate((11, 22, 33, 44)):
        rng = np.random.RandomState(seed)
        bases = "ACGT"
        out = {}
        while len(out) < 120:
            bc = "".join(bases[j] for j in rng.randint(0, 4, 8))
            out[bc] = int(max(1, rng.geometric(1.0 / 400)))
        cases.append(("geometric_%d" % i, out))
    return cases


def compare_distance_knee(rec, cases):
    mismatches = 0
    for name, counts in cases:
        counter = collections.Counter(counts)
        try:
            py = sorted(py_wl.getKneeEstimateDistance(counter))
        except Exception as exc:  # noqa: BLE001
            rec.log("knee_distance[%s]" % name,
                    "python raised %s" % type(exc).__name__)
            continue
        got = rs.get_cell_whitelist(dict(counts), "distance")
        got = sorted(got) if got else []
        if py != got:
            mismatches += 1
            rec.log("knee_distance[%s]" % name,
                    "python selected %d, rust %d" % (len(py), len(got)))
    rec.record("knee_estimate_distance", "mismatched_distributions",
               value=mismatches, threshold=0, tier=1, passed=mismatches == 0)


def compare_density_knee(rec, cases):
    """Tier 2: the KDE is SciPy's on one side and this crate's on the other."""
    mismatches = 0
    compared = 0
    for name, counts in cases:
        counter = collections.Counter(counts)
        try:
            py = py_wl.getKneeEstimateDensity(counter)
        except Exception as exc:  # noqa: BLE001
            rec.log("knee_density[%s]" % name,
                    "python raised %s" % type(exc).__name__)
            continue
        try:
            got = rs.get_cell_whitelist(dict(counts), "density")
        except Exception as exc:  # noqa: BLE001
            rec.log("knee_density[%s]" % name,
                    "rust raised %s" % type(exc).__name__)
            mismatches += 1
            continue
        compared += 1
        py = sorted(py) if py else []
        got = sorted(got) if got else []
        if py != got:
            mismatches += 1
            rec.log("knee_density[%s]" % name,
                    "python selected %d, rust %d" % (len(py), len(got)))
    rec.record("knee_estimate_density", "mismatched_distributions",
               value=mismatches, threshold=0, tier=2, passed=mismatches == 0)
    rec.log("knee_estimate_density", "%d distributions compared" % compared)


def compare_error_correction(rec, cases):
    mismatches = 0
    for name, counts in cases:
        counter = collections.Counter(counts)
        try:
            whitelist = py_wl.getKneeEstimateDistance(counter)
        except Exception:  # noqa: BLE001
            continue
        py = {k: sorted(v) for k, v in py_wl.getErrorCorrectMapping(
            list(counts.keys()), set(whitelist), 1).items()}
        got = {k: sorted(v) for k, v in rs.get_error_correct_mapping(
            list(counts.keys()), list(whitelist), 1).items()}
        if py != got:
            mismatches += 1
            rec.log("error_correct_mapping[%s]" % name,
                    "python %d parents, rust %d" % (len(py), len(got)))
    rec.record("getErrorCorrectMapping", "mismatched_distributions",
               value=mismatches, threshold=0, tier=1, passed=mismatches == 0)


def compare_check_error(rec):
    cases = [
        ("AAAA", ["AAAT", "AATT", "AAAA", "AAAG"], 1),
        ("AAAA", ["CCCC", "GGGG"], 1),
        ("AAAA", ["AAA", "AAAAA"], 1),
        ("AAAA", ["AATT"], 2),
        ("ACGT", ["ACGA", "ACGC", "ACGG"], 1),
        ("ACGT", [], 1),
        ("ACGTACGT", ["ACGTACGA", "ACGTACGT", "TTTTTTTT"], 3),
    ]
    mismatches = 0
    for barcode, whitelist, errors in cases:
        py = py_wl.checkError(barcode, whitelist, errors)
        got = rs.check_error(barcode, whitelist, errors)
        if py != got:
            mismatches += 1
    rec.record("checkError", "mismatches", value=mismatches, threshold=0,
               tier=1, passed=mismatches == 0)


def compare_user_defined(rec):
    mismatches = 0
    with tempfile.TemporaryDirectory() as d:
        plain = Path(d) / "wl.tsv"
        plain.write_text("# comment\nAAAA\tAAAT,AATT\nCCCC\tCCCG\n")
        derive = Path(d) / "wl2.tsv"
        derive.write_text("AC\nGT\n")

        for path, kwargs, rs_kwargs in [
            (plain, {}, {}),
            (plain, dict(getErrorCorrection=True),
             dict(get_error_correction=True)),
            (derive, dict(deriveErrorCorrection=True, threshold=1),
             dict(derive_error_correction=True, threshold=1)),
        ]:
            py_wl_set, py_map = py_wl.getUserDefinedBarcodes(str(path), **kwargs)
            got_wl, got_map = rs.get_user_defined_barcodes(str(path), **rs_kwargs)
            if sorted(py_wl_set) != sorted(got_wl):
                mismatches += 1
            if (py_map or {}) != (got_map or {}):
                mismatches += 1
    rec.record("getUserDefinedBarcodes", "mismatches", value=mismatches,
               threshold=0, tier=1, passed=mismatches == 0)


def compare_error_detect(rec, cases):
    mismatches = 0
    for name, counts in cases:
        counter = collections.Counter(counts)
        try:
            whitelist = py_wl.getKneeEstimateDistance(counter)
        except Exception:  # noqa: BLE001
            continue
        for method in ("discard", "correct"):
            py_out, _ = py_wl.errorDetectAboveThreshold(
                counter, set(whitelist),
                collections.defaultdict(set), errors=1,
                resolution_method=method)
            got = rs.error_detect_above_threshold(
                dict(counts), list(whitelist), 1, method)
            if sorted(py_out) != sorted(got):
                mismatches += 1
                rec.log("errorDetectAboveThreshold[%s/%s]" % (name, method),
                        "python kept %d, rust %d" % (len(py_out), len(got)))
    rec.record("errorDetectAboveThreshold", "mismatched_distributions",
               value=mismatches, threshold=0, tier=1, passed=mismatches == 0)


def main():
    rec = ResultRecorder("whitelist")
    cases = distributions()
    rec.log("setup", "%d barcode count distributions" % len(cases))
    rec.log("setup", "scipy %s, numpy %s"
            % (__import__("scipy").__version__, np.__version__))

    compare_distance_knee(rec, cases)
    compare_density_knee(rec, cases)
    compare_error_correction(rec, cases)
    compare_check_error(rec)
    compare_user_defined(rec)
    compare_error_detect(rec, cases)
    rec.save()
    return 0


if __name__ == "__main__":
    sys.exit(main())
