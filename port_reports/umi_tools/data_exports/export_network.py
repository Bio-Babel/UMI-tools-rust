"""Export Tier-1 reference values for core::network.

Runs the *live* Python oracle over deterministic and randomised UMI count
dictionaries, dumping the full intermediate chain (adjacency list -> connected
components -> groups) rather than only the final result, so a Rust failure can
be localised to the stage that diverged.

Covers both the `<= 25 UMIs` pairwise path and the `> 25 UMIs`
substring-index path, which take completely different code in
`_get_adj_list_*`.
"""
import itertools
import json
import os
import random

import umi_tools.network as network

OUT = os.path.dirname(os.path.abspath(__file__))
METHODS = ["unique", "percentile", "cluster", "adjacency", "directional"]


def b(x):
    return x.decode() if isinstance(x, bytes) else x


def make_cases():
    """Deterministic cases first, then seeded random ones."""
    cases = []

    # The canonical example from the upstream repo's step1_unit_test.py,
    # re-expressed as bytes for the 1.1.6 API.
    cases.append({
        "name": "upstream_example",
        "counts": {b"ACGT": 456, b"AAAT": 90, b"ACAT": 72,
                   b"TCGT": 2, b"CCGT": 2, b"ACAG": 1},
        "threshold": 1,
    })
    cases.append({"name": "single_umi", "counts": {b"ACGT": 7}, "threshold": 1})
    cases.append({"name": "two_identical_counts",
                  "counts": {b"ACGT": 5, b"ACGA": 5}, "threshold": 1})
    cases.append({"name": "no_neighbours",
                  "counts": {b"AAAA": 3, b"CCCC": 3, b"GGGG": 3}, "threshold": 1})
    cases.append({"name": "chain",
                  "counts": {b"AAAA": 100, b"AAAT": 40, b"AATT": 10, b"ATTT": 2},
                  "threshold": 1})
    cases.append({"name": "threshold_two",
                  "counts": {b"AAAA": 100, b"AATT": 10, b"TTTT": 1},
                  "threshold": 2})
    cases.append({"name": "threshold_zero",
                  "counts": {b"AAAA": 5, b"AAAT": 5}, "threshold": 0})
    # Exactly at the >25 boundary that switches on the substring index.
    for n in (25, 26, 40):
        alphabet = b"ACGT"
        umis = [bytes(p) for p in itertools.islice(
            itertools.product(alphabet, repeat=4), n)]
        cases.append({
            "name": "boundary_%d" % n,
            "counts": {u: (i * 7) % 23 + 1 for i, u in enumerate(umis)},
            "threshold": 1,
        })

    rng = random.Random(20260818)
    for trial in range(60):
        n = rng.randint(2, 60)
        length = rng.choice([4, 6, 8, 10])
        umis = set()
        while len(umis) < n:
            umis.add(bytes(rng.choice(b"ACGT") for _ in range(length)))
        cases.append({
            "name": "random_%03d" % trial,
            "counts": {u: rng.randint(1, 500) for u in sorted(umis)},
            "threshold": rng.choice([1, 1, 1, 2]),
        })
    return cases


def dump_helpers():
    rows = {"get_substr_slices": [], "breadth_first_search": [],
            "remove_umis": [], "build_substr_idx": [],
            "iter_nearest_neighbours": []}

    for umi_length in range(1, 13):
        for idx_size in range(1, 5):
            rows["get_substr_slices"].append({
                "umi_length": umi_length, "idx_size": idx_size,
                "slices": [list(s) for s in
                           network.get_substr_slices(umi_length, idx_size)],
            })

    adj = {b"A": [b"B", b"C"], b"B": [b"A"], b"C": [b"A", b"D"],
           b"D": [b"C"], b"E": []}
    for node in sorted(adj):
        rows["breadth_first_search"].append({
            "node": b(node),
            "adj": {b(k): [b(x) for x in v] for k, v in adj.items()},
            "component": sorted(b(x) for x in
                                network.breadth_first_search(node, adj)),
        })

    for nodes in ([b"A"], [b"C"], [b"A", b"C"], [b"E"]):
        rows["remove_umis"].append({
            "nodes": [b(x) for x in nodes],
            "cluster": [b(x) for x in sorted(adj)],
            "remaining": sorted(b(x) for x in
                                network.remove_umis(adj, sorted(adj), nodes)),
        })

    umis = [b"AAAA", b"AAAT", b"AATT", b"ATTT", b"TTTT", b"CCCC"]
    for min_edit in (1, 2):
        idx = network.build_substr_idx(umis, 4, min_edit)
        rows["build_substr_idx"].append({
            "umis": [b(u) for u in umis], "umi_length": 4, "min_edit": min_edit,
            "index": {str(list(k)): {b(sub): sorted(b(x) for x in members)
                                     for sub, members in v.items()}
                      for k, v in idx.items()},
        })
        rows["iter_nearest_neighbours"].append({
            "umis": [b(u) for u in umis], "min_edit": min_edit,
            "pairs": sorted([b(a), b(x)] for a, x in
                            network.iter_nearest_neighbours(umis, idx)),
        })
    return rows


def dump_clustering(cases):
    rows = []
    for case in cases:
        counts = case["counts"]
        threshold = case["threshold"]
        entry = {"name": case["name"],
                 "counts": {b(k): v for k, v in counts.items()},
                 "threshold": threshold,
                 "methods": {}}
        # Rust's CountMap is a BTreeMap, so it iterates in sorted-key order
        # while Python's dict iterates in insertion order. Feeding Python a
        # dict rebuilt in sorted-key order makes the two orders identical, so
        # the comparison is exact rather than order-tolerant. See
        # 05_design.md § Determinism decisions.
        counts = {k: counts[k] for k in sorted(counts)}
        entry["counts"] = {b(k): v for k, v in counts.items()}
        for method in METHODS:
            c = network.UMIClusterer(method)
            umis = list(counts.keys())
            adj = c.get_adj_list(umis, counts, threshold)
            clusters = c.get_connected_components(umis, adj, counts)
            groups = [list(x) for x in c.get_groups(clusters, adj, counts)]
            final = c(counts, threshold)
            entry["methods"][method] = {
                "adj_list": (None if adj is None else
                             {b(k): sorted(b(x) for x in v)
                              for k, v in adj.items()}),
                "clusters": (None if clusters is None else
                             [[b(x) for x in cl] if isinstance(cl, list)
                              else b(cl) for cl in clusters]),
                "groups": [[b(x) for x in g] for g in groups],
                "final": [[b(x) for x in g] for g in final],
            }
        rows.append(entry)
    return rows


def dump_cell_clusterer():
    """CellClusterer's __call__ is broken upstream; export what works."""
    rows = []
    counts = {b"ACGT": 100, b"ACGA": 2, b"TTTT": 50}
    cc = network.CellClusterer("directional", dir_threshold=10, fuzzy_match=False)
    adj = cc._get_adj_list_directional(list(counts.keys()), counts)
    rows.append({"fuzzy": False,
                 "counts": {b(k): v for k, v in counts.items()},
                 "adj_list": {b(k): sorted(b(x) for x in v)
                              for k, v in adj.items()}})

    scounts = {"ACGT": 100, "ACGA": 2, "ATTC": 1, "ACG": 3, "ACGTT": 4}
    cc2 = network.CellClusterer("directional", dir_threshold=10, fuzzy_match=True)
    adj2 = cc2._get_adj_list_directional(list(scounts.keys()), scounts)
    rows.append({"fuzzy": True, "dir_threshold": 10, "counts": scounts,
                 "adj_list": {k: sorted(v) for k, v in adj2.items()}})

    # __call__ is expected to raise; record the exact failure.
    try:
        cc(list(counts.keys()), counts)
        call_error = None
    except Exception as exc:  # noqa: BLE001 - recording upstream behaviour
        call_error = "%s: %s" % (type(exc).__name__, exc)
    return {"cases": rows, "call_error": call_error}


def main():
    cases = make_cases()
    payload = {
        "helpers": dump_helpers(),
        "clustering": dump_clustering(cases),
        "cell_clusterer": dump_cell_clusterer(),
    }
    path = os.path.join(OUT, "reference_network.json")
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=1)
    print("wrote %s (%d clustering cases)" % (path, len(payload["clustering"])))


if __name__ == "__main__":
    main()
