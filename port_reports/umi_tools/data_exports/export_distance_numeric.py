"""Export Tier-1 reference values for core::distance and core::numeric.

Runs the *live* Python oracle (not a transcription of it) and dumps
{args, result} records as JSON for consumption by `cargo test` fixtures.

The `edit_distance` unequal-length case is exported as the `.pyx` gold
standard's intent (infinity) rather than the NumPy-2 `AttributeError` the
checked-in generated `.c` actually raises -- see 05_design.md deviation D1.
"""
import json
import os
import sys

import numpy as np
from umi_tools._dedup_umi import edit_distance

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)))

PAIRS = [
    (b"ACGT", b"ACGT"), (b"ACGT", b"ACGA"), (b"ACGT", b"TTTT"),
    (b"", b""), (b"A", b"A"), (b"A", b"T"),
    (b"ACGTACGTACGT", b"ACGTACGTACGA"),
    (b"NNNN", b"ACGT"), (b"\x00\x01", b"\x00\x02"),
    (b"ACGT", b"ACG"), (b"", b"A"), (b"AAAA", b"AAAAA"),
]


def dump_edit_distance():
    rows = []
    for a, b in PAIRS:
        if len(a) == len(b):
            d = edit_distance(a, b)
            rows.append({"a": list(a), "b": list(b), "distance": d})
        else:
            # gold standard = _dedup_umi.pyx, which returns np.inf here.
            rows.append({"a": list(a), "b": list(b), "distance": None})
    return rows


VECTORS = [
    [1.0], [1.0, 2.0], [1.0, 2.0, 3.0], [3.0, 1.0, 2.0],
    [1.0, 2.0, 3.0, 4.0], [-1.0, 0.0, 1.0],
    [1e-12, 1e12], [0.0, 0.0, 0.0],
    [5.0, 3.0, 8.0, 1.0, 9.0, 2.0, 7.0],
    [float(x) for x in range(1, 101)],
]


def dump_numeric():
    rows = []
    for v in VECTORS:
        a = np.array(v)
        rows.append({
            "input": v,
            "median": float(np.median(a)),
            "cumsum": [float(x) for x in np.cumsum(a)],
            "argmax": int(np.argmax(a)),
            "sum": float(np.sum(a)),
            "max": float(np.max(a)),
            "min": float(np.min(a)),
            "log10": [float(x) for x in np.log10(a)] if all(x > 0 for x in v) else None,
        })
    return rows


def dump_linspace():
    rows = []
    for lo, hi, n in [(0.0, 1.0, 5), (1.0, 3.0, 10000), (-2.0, 2.0, 3),
                      (0.5, 0.5, 4), (0.0, 1.0, 1), (0.0, 1.0, 2)]:
        xs = np.linspace(lo, hi, n)
        rows.append({"lo": lo, "hi": hi, "n": n,
                     "first": float(xs[0]), "last": float(xs[-1]),
                     "sum": float(xs.sum()),
                     "head": [float(x) for x in xs[:5]]})
    return rows


def dump_digitize_bincount():
    rows = []
    cases = [
        ([-1.0, 0.0, 0.5, 1.0, 2.0, 3.0], list(range(-1, 5))),
        ([0.0, 1.0, 1.5, 2.0], list(range(-1, 4))),
        ([-1.0], list(range(-1, 2))),
        ([2.5, 2.5, 2.5], list(range(-1, 5))),
    ]
    for vals, bins in cases:
        binned = np.digitize(vals, bins, right=True)
        rows.append({
            "values": vals,
            "bins": bins,
            "digitize_right": [int(x) for x in binned],
            "bincount_minlength_6": [int(x) for x in np.bincount(binned, minlength=6)],
            "bincount_default": [int(x) for x in np.bincount(binned)],
        })
    return rows


def main():
    payload = {
        "edit_distance": dump_edit_distance(),
        "numeric": dump_numeric(),
        "linspace": dump_linspace(),
        "digitize_bincount": dump_digitize_bincount(),
        "numpy_version": np.__version__,
        "python_version": sys.version,
    }
    path = os.path.join(OUT, "reference_distance_numeric.json")
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=1)
    print("wrote", path)


if __name__ == "__main__":
    main()
