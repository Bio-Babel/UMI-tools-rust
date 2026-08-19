"""Export reference values for core::whitelist, core::bktree and core::fuzzy.

Dumps both knee estimators' full intermediate chain (cumulative counts,
distance-to-line curve, chosen index) alongside the selected barcode set, so a
Rust divergence can be attributed to the geometry or to the selection rule.
"""
import collections
import json
import os
import tempfile

import numpy as np
import numpy.matlib as npm
import regex
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde

from umi_tools.whitelist_methods import (
    checkError, errorDetectAboveThreshold, getCellWhitelist,
    getErrorCorrectMapping, getKneeEstimateDensity, getKneeEstimateDistance,
    getUserDefinedBarcodes,
)

OUT = os.path.dirname(os.path.abspath(__file__))


def knee_distance_internals(values):
    """Mirror of getKneeDistance's internals, for stage-by-stage comparison."""
    nPoints = len(values)
    allCoord = np.vstack((range(nPoints), values)).T
    firstPoint = allCoord[0]
    lineVec = allCoord[-1] - allCoord[0]
    lineVecNorm = lineVec / np.sqrt(np.sum(lineVec ** 2))
    vecFromFirst = allCoord - firstPoint
    scalarProduct = np.sum(
        vecFromFirst * npm.repmat(lineVecNorm, nPoints, 1), axis=1)
    vecFromFirstParallel = np.outer(scalarProduct, lineVecNorm)
    vecToLine = vecFromFirst - vecFromFirstParallel
    distToLine = np.sqrt(np.sum(vecToLine ** 2, axis=1))
    return [float(x) for x in distToLine], int(np.argmax(distToLine))


DISTRIBUTIONS = [
    ("bimodal_small", {"AAAA": 1000, "AAAT": 800, "AATT": 600, "ATTT": 400,
                       "TTTT": 200, "CCCC": 10, "CCCG": 8, "CCGG": 6,
                       "CGGG": 4, "GGGG": 2}),
    ("steep", {"A" * 4: 10000, "AAAT": 5, "AATT": 4, "ATTT": 3, "TTTT": 2,
               "CCCC": 1}),
    ("flat", {b: 100 for b in ["AAAA", "CCCC", "GGGG", "TTTT", "ACGT",
                               "TGCA", "AACC", "GGTT"]}),
]


def geometric_distribution(n, seed):
    rng = np.random.RandomState(seed)
    bases = "ACGT"
    out = {}
    while len(out) < n:
        bc = "".join(bases[i] for i in rng.randint(0, 4, 8))
        out[bc] = int(max(1, rng.geometric(1.0 / 400)))
    return out


for i, seed in enumerate((11, 22, 33)):
    DISTRIBUTIONS.append(("geometric_%d" % i,
                          geometric_distribution(120, seed)))


def dump_knee():
    rows = []
    for name, counts in DISTRIBUTIONS:
        c = collections.Counter(counts)
        ordered = c.most_common()
        raw_counts = [x[1] for x in ordered]
        values = list(np.cumsum(raw_counts))
        dist, idx = knee_distance_internals(values)

        entry = {
            "name": name,
            "counts": dict(counts),
            "most_common": [[k, v] for k, v in ordered],
            "cumulative": [int(v) for v in values],
            "dist_to_line": dist,
            "argmax": idx,
        }
        try:
            entry["distance_whitelist"] = sorted(getKneeEstimateDistance(c))
        except Exception as exc:  # noqa: BLE001
            entry["distance_whitelist"] = None
            entry["distance_error"] = "%s: %s" % (type(exc).__name__, exc)

        for cell_number in (None, 3):
            key = "distance_whitelist_n%s" % cell_number
            try:
                wl = getKneeEstimateDistance(c, cell_number or False)
                entry[key] = sorted(wl) if wl else None
            except Exception as exc:  # noqa: BLE001
                entry[key] = None
                entry[key + "_error"] = "%s: %s" % (type(exc).__name__, exc)

        # Density knee: also dump the KDE evaluation so the Tier-2 comparison
        # can be scored on the curve, not only on the selected set.
        threshold = 0.001 * c.most_common(1)[0][1]
        sorted_counts = sorted(c.values(), reverse=True)
        thresh_counts = [x for x in sorted_counts if x > threshold]
        entry["density_log_counts"] = [float(x) for x in
                                       np.log10(thresh_counts)]
        if len(set(thresh_counts)) > 1:
            log_counts = np.log10(thresh_counts)
            density = gaussian_kde(log_counts, bw_method=0.1)
            xx = np.linspace(log_counts.min(), log_counts.max(), 10000)
            dens = density(xx)
            entry["density_xx_head"] = [float(x) for x in xx[:5]]
            entry["density_values_head"] = [float(x) for x in dens[:5]]
            entry["density_values_sample"] = [float(dens[i])
                                              for i in range(0, 10000, 500)]
            entry["density_sum"] = float(dens.sum())
            entry["density_covariance_factor"] = float(density.factor)
            entry["density_covariance"] = float(density.covariance[0][0])
            entry["local_mins"] = [int(x) for x in
                                   argrelextrema(dens, np.less)[0]]
            try:
                wl = getKneeEstimateDensity(c)
                entry["density_whitelist"] = sorted(wl) if wl else None
            except Exception as exc:  # noqa: BLE001
                entry["density_whitelist"] = None
                entry["density_error"] = "%s: %s" % (type(exc).__name__, exc)
        rows.append(entry)
    return rows


def dump_error_correction():
    rows = []
    counts = {"AAAA": 1000, "AAAT": 800, "AATT": 600, "ATTT": 400,
              "TTTT": 200, "CCCC": 10, "CCCG": 8}
    for whitelist in (["AAAA", "CCCC"], ["AAAA"], ["AAAA", "AAAT"]):
        for threshold in (1, 2):
            m = getErrorCorrectMapping(list(counts), set(whitelist), threshold)
            rows.append({
                "barcodes": list(counts), "whitelist": sorted(whitelist),
                "threshold": threshold,
                "mapping": {k: sorted(v) for k, v in m.items()},
            })
    return rows


def dump_check_error():
    rows = []
    cases = [
        ("AAAA", ["AAAT", "AATT", "AAAA", "AAAG"], 1),
        ("AAAA", ["CCCC", "GGGG"], 1),
        ("AAAA", ["AAA", "AAAAA"], 1),
        ("AAAA", ["AATT"], 2),
        ("ACGT", ["ACGA", "ACGC", "ACGG"], 1),
        ("ACGT", [], 1),
    ]
    for barcode, whitelist, errors in cases:
        rows.append({"barcode": barcode, "whitelist": whitelist,
                     "errors": errors,
                     "near_matches": checkError(barcode, whitelist, errors)})
    return rows


def dump_error_detect():
    rows = []
    c = collections.Counter({"AAAA": 1000, "AAAT": 800, "AATT": 600,
                             "ATTT": 400, "TTTT": 200})
    wl, t2f = getCellWhitelist(c, "distance", False, False, 1, None)
    for method in ("discard", "correct"):
        new_wl, new_map = errorDetectAboveThreshold(
            c, set(wl), collections.defaultdict(set, t2f or {}),
            errors=1, resolution_method=method)
        rows.append({
            "counts": dict(c), "whitelist_in": sorted(wl),
            "true_to_false_in": {k: sorted(v) for k, v in (t2f or {}).items()},
            "resolution_method": method,
            "whitelist_out": sorted(new_wl),
            "true_to_false_out": {k: sorted(v) for k, v in new_map.items()},
        })
    return rows


def dump_user_defined():
    rows = []
    with tempfile.TemporaryDirectory() as d:
        p1 = os.path.join(d, "wl1.tsv")
        with open(p1, "w") as fh:
            fh.write("# comment\nAAAA\tAAAT,AATT\nCCCC\tCCCG\n")
        p2 = os.path.join(d, "wl2.tsv")
        with open(p2, "w") as fh:
            fh.write("GG\nTT\n")

        wl, m = getUserDefinedBarcodes(p1)
        rows.append({"case": "plain", "whitelist": sorted(wl), "map": m})

        wl, m = getUserDefinedBarcodes(p1, getErrorCorrection=True)
        rows.append({"case": "get_error_correction", "whitelist": sorted(wl),
                     "map": m})

        wl, m = getUserDefinedBarcodes(p1, p2)
        rows.append({"case": "paired", "whitelist": sorted(wl), "map": m})

        p3 = os.path.join(d, "wl3.tsv")
        with open(p3, "w") as fh:
            fh.write("AC\nGT\n")
        wl, m = getUserDefinedBarcodes(p3, deriveErrorCorrection=True,
                                       threshold=1)
        rows.append({"case": "derive_threshold1", "whitelist": sorted(wl),
                     "map": {k: v for k, v in sorted(m.items())}})
    return rows


def dump_fuzzy():
    """The three fuzzy patterns umi_tools uses, over an exhaustive small set."""
    rows = []
    alphabet = "ACG"
    words = [""]
    for length in (1, 2, 3, 4):
        words += ["".join(w) for w in
                  __import__("itertools").product(alphabet, repeat=length)]
    words = words[:60]
    for pat in ("ACG", "AC", "ACGA"):
        for text in words:
            rows.append({
                "pattern": pat, "text": text,
                "e1": bool(regex.compile("(%s){e<=1}" % pat).match(text)),
                "s1": bool(regex.compile("(%s){s<=1}" % pat).match(text)),
                "i1": bool(regex.compile("(%s){i<=1}" % pat).match(text)),
                "e2": bool(regex.compile("(%s){e<=2}" % pat).match(text)),
            })
    return rows


def main():
    payload = {
        "knee": dump_knee(),
        "error_correct_mapping": dump_error_correction(),
        "check_error": dump_check_error(),
        "error_detect_above_threshold": dump_error_detect(),
        "user_defined_barcodes": dump_user_defined(),
        "fuzzy": dump_fuzzy(),
        "scipy_version": __import__("scipy").__version__,
        "numpy_version": np.__version__,
    }
    path = os.path.join(OUT, "reference_whitelist.json")
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=1)
    print("wrote", path)


if __name__ == "__main__":
    main()
