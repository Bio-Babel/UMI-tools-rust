"""Differential test: `umi_tools.network` (Python) vs `umi_tools_rs` (Rust).

Imports the original Python package and the `maturin develop`-installed Rust
extension into the **same interpreter** and calls both on the same inputs, so
the comparison involves no serialisation round-trip.

Both packages coexist because the extension is installed under its own
top-level name (`umi_tools_rs`) rather than as `umi_tools._native`; installing
it as a submodule of a package that is separately `pip install -e .`-ed
shadows the original, which is documented in
`references/differential-testing-and-validation.md` and recorded as deviation
D9 in `port_reports/umi_tools/05_design.md`.

Run after `maturin develop` from `crates/pyo3-bindings/`:

    python validation/differential_network.py
"""

import itertools
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _recorder import ResultRecorder  # noqa: E402

import umi_tools.network as py_network  # noqa: E402
import umi_tools_rs as rs  # noqa: E402
from umi_tools._dedup_umi import edit_distance as py_edit_distance  # noqa: E402

METHODS = ["unique", "percentile", "cluster", "adjacency", "directional"]


def random_cases(n_cases, seed):
    """Seeded random UMI count dictionaries, spanning both clustering paths.

    The `> 25 UMIs` case matters: it switches `_get_adj_list_*` from exhaustive
    pairwise comparison to the substring index, which is a different algorithm.
    """
    rng = random.Random(seed)
    cases = []
    for _ in range(n_cases):
        n = rng.randint(2, 60)
        length = rng.choice([4, 6, 8, 10])
        umis = set()
        while len(umis) < n:
            umis.add(bytes(rng.choice(b"ACGT") for _ in range(length)))
        # Sorted, so the Python dict's insertion order matches the Rust
        # BTreeMap's iteration order and the comparison is exact rather than
        # order-tolerant.
        cases.append({u: rng.randint(1, 500) for u in sorted(umis)})
    return cases


def compare_edit_distance(rec):
    pairs = [
        (b"ACGT", b"ACGT"), (b"ACGT", b"ACGA"), (b"ACGT", b"TTTT"),
        (b"", b""), (b"A", b"T"), (b"ACGTACGTACGT", b"ACGTACGTACGA"),
        (b"\x00\x01", b"\x00\x02"), (b"NNNN", b"ACGT"),
    ]
    mismatches = 0
    for a, b in pairs:
        # Equal lengths only: the shipped extension raises AttributeError on
        # unequal lengths under NumPy 2 (see 05_design.md deviation D1).
        if py_edit_distance(a, b) != rs.edit_distance(a, b):
            mismatches += 1
    rec.record("edit_distance", "mismatches", value=mismatches, threshold=0,
               tier=1, passed=mismatches == 0)

    # The unequal-length case has no live oracle, so it is asserted against the
    # `.pyx` source's documented intent instead.
    rec.record("edit_distance/unequal_lengths", "is_none",
               value=1 if rs.edit_distance(b"ACGT", b"ACG") is None else 0,
               threshold=1, tier=1)


def compare_clustering(rec, cases):
    for method in METHODS:
        mismatches = 0
        compared = 0
        for counts in cases:
            py = py_network.UMIClusterer(method)(dict(counts), 1)
            rs_clusterer = rs.UMIClusterer(method)
            got = rs_clusterer.cluster(dict(counts), 1)
            compared += 1
            if method == "adjacency":
                # The trailing members of an adjacency group come from a Python
                # `set` difference, whose order is an artefact; the leading
                # (parent) UMI is positionally exact either way.
                same = len(py) == len(got) and all(
                    p[0] == g[0] and sorted(p) == sorted(g)
                    for p, g in zip(py, got))
            else:
                same = py == got
            if not same:
                mismatches += 1
        rec.record("UMIClusterer[%s]" % method, "mismatched_positions",
                   value=mismatches, threshold=0, tier=1,
                   passed=mismatches == 0)
        rec.log("UMIClusterer[%s]" % method, "%d positions compared" % compared)


def compare_adjacency(rec, cases):
    """The adjacency list, before grouping, so a divergence is localised."""
    for method in ("adjacency", "directional"):
        mismatches = 0
        for counts in cases:
            clusterer = py_network.UMIClusterer(method)
            umis = list(counts.keys())
            py = {k: sorted(v) for k, v in
                  clusterer.get_adj_list(umis, counts, 1).items()}
            got = {k: sorted(v) for k, v in
                   rs.get_adj_list(dict(counts), 1, method).items()}
            if py != got:
                mismatches += 1
        rec.record("get_adj_list[%s]" % method, "mismatched_positions",
                   value=mismatches, threshold=0, tier=1,
                   passed=mismatches == 0)


def compare_average_distance(rec, cases):
    worst = 0.0
    for counts in cases:
        umis = list(counts.keys())
        py = py_network.umi_methods.get_average_umi_distance(umis) \
            if hasattr(py_network, "umi_methods") else None
        if py is None:
            import umi_tools.umi_methods as um
            py = um.get_average_umi_distance(umis)
        got = rs.get_average_umi_distance(umis)
        worst = max(worst, abs(py - got))
    rec.record("get_average_umi_distance", "max_abs_diff", value=worst,
               threshold=1e-9, tier=1, passed=worst <= 1e-9)


def main():
    rec = ResultRecorder("network")
    cases = random_cases(200, seed=20260818)
    rec.log("setup", "%d randomised UMI count dictionaries" % len(cases))
    rec.log("setup", "python umi_tools %s, rust umi_tools_rs %s"
            % (__import__("umi_tools").__version__, rs.version()))

    compare_edit_distance(rec)
    compare_adjacency(rec, cases)
    compare_clustering(rec, cases)
    compare_average_distance(rec, cases)
    rec.save()
    return 0


if __name__ == "__main__":
    sys.exit(main())
