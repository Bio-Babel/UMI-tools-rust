"""Export reference extraction results for core::extract.

Runs the live ExtractFilterAndUpdate over real reads from the checked-in
FASTQ fixtures, for every --bc-pattern/--extract-method combination the test
suite exercises, and records the full rewritten record plus the read_counts
tally. That tally is what `umi_tools extract` prints at the end of a run, so
it is part of the observable output, not just an internal.
"""
import gzip
import json
import os

import regex

import umi_tools.extract_methods as em
from umi_tools.umi_methods import Record

OUT = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(os.environ.get("UMI_TOOLS_SRC", "../umi_tools_src"), "tests")


def reads(filename, limit):
    path = os.path.join(TESTS, filename)
    if not os.path.exists(path):
        return []
    opener = gzip.open if filename.endswith(".gz") else open
    out, buf = [], []
    with opener(path, "rt") as fh:
        for line in fh:
            buf.append(line.rstrip("\n"))
            if len(buf) == 4:
                out.append(Record(buf[0][1:], buf[1], buf[3]))
                buf = []
                if len(out) >= limit:
                    break
    return out


def rec(r):
    return {"identifier": r.identifier, "seq": r.seq, "quals": r.quals}


CONFIGS = [
    # (label, kwargs, read-source, paired-source)
    ("string_simple", dict(method="string", pattern="NNNXXXXNN"), "slim.fastq.gz", None),
    ("string_3prime", dict(method="string", pattern="NNNNN", prime3=True),
     "slim.fastq.gz", None),
    ("string_cell", dict(method="string", pattern="CCCCCCNNNNNNNNNN", extract_cell=True),
     "scrb_seq_fastq.1.gz", None),
    ("regex_5prime", dict(method="regex",
                          pattern=regex.compile("^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})")),
     "slim.fastq.gz", None),
    ("regex_3prime", dict(method="regex",
                          pattern=regex.compile(".*(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})$")),
     "slim.fastq.gz", None),
    ("regex_cell", dict(method="regex",
                        pattern=regex.compile("^(?P<cell_1>.{6})(?P<umi_1>.{10})"),
                        extract_cell=True),
     "scrb_seq_fastq.1.gz", None),
    ("regex_retain_umi", dict(method="regex",
                              pattern=regex.compile("^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})"),
                              retain_umi=True),
     "slim.fastq.gz", None),
    ("regex_quality_filter", dict(method="regex",
                                  pattern=regex.compile("^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})"),
                                  quality_encoding="phred33",
                                  quality_filter_threshold=30),
     "slim.fastq.gz", None),
    ("regex_quality_mask", dict(method="regex",
                                pattern=regex.compile("^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})"),
                                quality_encoding="phred33",
                                quality_filter_mask=30),
     "slim.fastq.gz", None),
    ("regex_indrop_fuzzy", dict(
        method="regex",
        pattern=regex.compile(
            "(?P<cell_1>.{8,12})(?P<discard_2>GAGTGATTGCTTGTGACGCCTT{s<=2})"
            "(?P<cell_3>.{8})(?P<umi_1>.{6})T{3}.*"),
        extract_cell=True),
     "indrop.fastq.1.gz", None),
]

PAIRED = [
    ("string_paired", dict(method="string", pattern="NNN", pattern2="NNN"),
     "indrop.fastq.1.gz", "indrop.fastq.2.gz"),
    ("regex_scrb_paired", dict(
        method="regex",
        pattern=regex.compile("^(?P<cell_1>.{6})(?P<umi_1>.{10})"),
        extract_cell=True),
     "scrb_seq_fastq.1.gz", "scrb_seq_fastq.2.gz"),
]


def run_single(label, kwargs, source, limit=120):
    e = em.ExtractFilterAndUpdate(**kwargs)
    rows = []
    for r in reads(source, limit):
        before = rec(r)
        try:
            out = e(Record(r.identifier, r.seq, r.quals))
            rows.append({"input": before,
                         "output": rec(out) if out else None,
                         "error": None})
        except Exception as exc:  # noqa: BLE001 - recording oracle behaviour
            rows.append({"input": before, "output": None,
                         "error": "%s: %s" % (type(exc).__name__, exc)})
    return {"label": label, "cases": rows, "read_counts": dict(e.read_counts)}


def run_paired(label, kwargs, s1, s2, limit=80):
    e = em.ExtractFilterAndUpdate(**kwargs)
    rows = []
    for r1, r2 in zip(reads(s1, limit), reads(s2, limit)):
        b1, b2 = rec(r1), rec(r2)
        try:
            out = e(Record(r1.identifier, r1.seq, r1.quals),
                    Record(r2.identifier, r2.seq, r2.quals))
            rows.append({"input1": b1, "input2": b2,
                         "output1": rec(out[0]) if out else None,
                         "output2": rec(out[1]) if out else None,
                         "error": None})
        except Exception as exc:  # noqa: BLE001
            rows.append({"input1": b1, "input2": b2, "output1": None,
                         "output2": None,
                         "error": "%s: %s" % (type(exc).__name__, exc)})
    return {"label": label, "cases": rows, "read_counts": dict(e.read_counts)}


def main():
    payload = {
        "single": [run_single(l, k, s) for l, k, s, _ in CONFIGS],
        "paired": [run_paired(l, k, a, b) for l, k, a, b in PAIRED],
    }
    total = sum(len(g["cases"]) for g in payload["single"]) + \
        sum(len(g["cases"]) for g in payload["paired"])
    path = os.path.join(OUT, "reference_extract.json")
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=1)
    print("wrote %s (%d read comparisons)" % (path, total))


if __name__ == "__main__":
    main()
