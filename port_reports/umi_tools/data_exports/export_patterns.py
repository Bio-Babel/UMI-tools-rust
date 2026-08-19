"""Export reference match results for core::pattern.

Runs every --bc-pattern that appears in the package's test suite and docs
against real reads drawn from the checked-in FASTQ fixtures, recording the
per-named-group spans. Those spans, not merely whether the pattern matched,
are what core::extract depends on.
"""
import gzip
import itertools
import json
import os

import regex

OUT = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(os.environ.get("UMI_TOOLS_SRC", "../umi_tools_src"), "tests")

PATTERNS = [
    "^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})",
    ".*(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})$",
    "^(?P<cell_1>.{6})(?P<umi_1>.{10})",
    "^(?P<umi_1>.{10}).*",
    "(?P<cell_1>.{8,12})(?P<discard_2>GAGTGATTGCTTGTGACGCCTT)(?P<cell_3>.{8})(?P<umi_1>.{6})T{3}.*",
    "(?P<cell_1>.{8,12})(?P<discard_2>GAGTGATTGCTTGTGACGCCTT{s<=2})(?P<cell_3>.{8})(?P<umi_1>.{6})T{3}.*",
    "(?P<discard_1>.{8,12})(?P<discard_2>GAGTGATTGCTTGTGACGCCTT{s<=2})(?P<discard_3>.{8})(?P<umi_1>.{6})T{3}.*",
]


def reads(filename, limit):
    path = os.path.join(TESTS, filename)
    if not os.path.exists(path):
        return []
    opener = gzip.open if filename.endswith(".gz") else open
    out = []
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                out.append(line.rstrip("\n"))
                if len(out) >= limit:
                    break
    return out


def main():
    subjects = []
    for fn, n in [("slim.fastq.gz", 60), ("scrb_seq_fastq.1.gz", 60),
                  ("indrop.fastq.1.gz", 80), ("indrop.fastq.2.gz", 40)]:
        subjects += reads(fn, n)
    # Plus synthetic edge cases: empty, very short, all-N.
    subjects += ["", "A", "ACGT", "N" * 30, "ACGT" * 20]

    rows = []
    for pattern in PATTERNS:
        compiled = regex.compile(pattern)
        for subject in subjects:
            m = compiled.match(subject)
            rows.append({
                "pattern": pattern,
                "subject": subject,
                "matched": bool(m),
                "span": list(m.span()) if m else None,
                "groups": ({k: list(m.span(k)) for k in compiled.groupindex
                            if m.span(k) != (-1, -1)} if m else None),
            })

    path = os.path.join(OUT, "reference_patterns.json")
    with open(path, "w") as fh:
        json.dump({"cases": rows}, fh, indent=1)
    print("wrote %s (%d cases, %d subjects)" % (path, len(rows), len(subjects)))


if __name__ == "__main__":
    main()
