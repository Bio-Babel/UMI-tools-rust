"""Differential test: `umi_tools.extract_methods` vs `umi_tools_rs`.

Runs both implementations over real reads from the package's own FASTQ
fixtures, in one interpreter. See `differential_network.py`'s docstring for
why the extension is imported as `umi_tools_rs`.
"""

import gzip
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _recorder import ResultRecorder  # noqa: E402

import regex  # noqa: E402
import umi_tools.extract_methods as py_extract  # noqa: E402
import umi_tools_rs as rs  # noqa: E402
from umi_tools.umi_methods import Record  # noqa: E402

TESTS = Path(os.environ.get("UMI_TOOLS_SRC", "../umi_tools_src")) / "tests"

INDROP = ("(?P<cell_1>.{8,12})(?P<discard_2>GAGTGATTGCTTGTGACGCCTT{s<=2})"
          "(?P<cell_3>.{8})(?P<umi_1>.{6})T{3}.*")
FIVE_PRIME = "^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})"
THREE_PRIME = ".*(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})$"


def reads(filename, limit):
    path = TESTS / filename
    if not path.exists():
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


def config(label, fixture, method, pattern, **extra):
    """One comparison, with the pattern compiled for Python and raw for Rust."""
    py_kwargs = dict(method=method, **extra)
    rs_kwargs = dict(method=method, **extra)
    if pattern is not None:
        py_kwargs["pattern"] = (regex.compile(pattern) if method == "regex"
                                else pattern)
        rs_kwargs["pattern"] = pattern
    return label, py_kwargs, rs_kwargs, fixture


CONFIGS = [
    config("string_simple", "slim.fastq.gz", "string", "NNNXXXXNN"),
    config("string_3prime", "slim.fastq.gz", "string", "NNNNN", prime3=True),
    config("string_cell", "scrb_seq_fastq.1.gz", "string",
           "CCCCCCNNNNNNNNNN", extract_cell=True),
    config("regex_5prime", "slim.fastq.gz", "regex", FIVE_PRIME),
    config("regex_3prime", "slim.fastq.gz", "regex", THREE_PRIME),
    config("regex_retain_umi", "slim.fastq.gz", "regex", FIVE_PRIME,
           retain_umi=True),
    config("regex_quality_mask", "slim.fastq.gz", "regex", FIVE_PRIME,
           quality_encoding="phred33", quality_filter_mask=30),
    config("regex_indrop_fuzzy", "indrop.fastq.1.gz", "regex", INDROP,
           extract_cell=True),
]


def compare_config(rec, label, py_kwargs, rs_kwargs, fixture, limit=150):
    py_extractor = py_extract.ExtractFilterAndUpdate(**py_kwargs)
    rs_extractor = rs.ExtractFilterAndUpdate(**rs_kwargs)

    mismatches = 0
    compared = 0
    for read in reads(fixture, limit):
        py_out = py_extractor(Record(read.identifier, read.seq, read.quals))
        rs_out = rs_extractor.process(read.identifier, read.seq, read.quals)
        compared += 1

        if py_out is None:
            if rs_out is not None:
                mismatches += 1
            continue
        if rs_out is None:
            mismatches += 1
            continue
        if (py_out.identifier, py_out.seq, py_out.quals) != tuple(rs_out):
            mismatches += 1

    rec.record("ExtractFilterAndUpdate[%s]" % label, "mismatched_reads",
               value=mismatches, threshold=0, tier=1, passed=mismatches == 0)
    rec.log("ExtractFilterAndUpdate[%s]" % label,
            "%d reads compared" % compared)

    # The event tallies are printed at the end of a run, so they are output too.
    py_counts = dict(py_extractor.read_counts)
    rs_counts = rs_extractor.get_read_counts()
    differing = sum(1 for k in set(py_counts) | set(rs_counts)
                    if py_counts.get(k, 0) != rs_counts.get(k, 0))
    rec.record("read_counts[%s]" % label, "differing_events", value=differing,
               threshold=0, tier=1, passed=differing == 0)


def compare_helpers(rec):
    quality_cases = [
        ("ACGT", "I!IJ", "phred33", 30),
        ("ACGTACGT", "IIIIIIII", "phred33", 30),
        ("ACG", "abc", "solexa", 10),
        ("ACG", "abc", "phred64", 10),
    ]
    mask_mismatches = below_mismatches = 0
    for umi, quals, encoding, threshold in quality_cases:
        if py_extract.mask_umi(umi, quals, encoding, threshold) != \
                rs.mask_umi(umi, quals, encoding, threshold):
            mask_mismatches += 1
        if bool(py_extract.umi_below_threshold(quals, encoding, threshold)) != \
                rs.umi_below_threshold(quals, encoding, threshold):
            below_mismatches += 1
    rec.record("mask_umi", "mismatches", value=mask_mismatches, threshold=0,
               tier=1, passed=mask_mismatches == 0)
    rec.record("umi_below_threshold", "mismatches", value=below_mismatches,
               threshold=0, tier=1, passed=below_mismatches == 0)

    splits = 0
    for seq, quals, umi_b, cell_b, discard_b, retain in [
        ("ACGTACGT", "12345678", [0, 1], [2, 3], [4], False),
        ("ACGTACGT", "12345678", [0, 1], [2, 3], [4], True),
        ("ACGT", "1234", [], [], [], False),
    ]:
        py = py_extract.extractSeqAndQuals(
            seq, quals, set(umi_b), set(cell_b), set(discard_b), retain)
        got = rs.extract_seq_and_quals(
            seq, quals, umi_b, cell_b, discard_b, retain)
        if tuple(py) != tuple(got):
            splits += 1
    rec.record("extractSeqAndQuals", "mismatches", value=splits, threshold=0,
               tier=1, passed=splits == 0)


def compare_patterns(rec):
    """The regex group spans, which `extract` slices the read by."""
    subjects = [r.seq for r in reads("slim.fastq.gz", 60)]
    subjects += [r.seq for r in reads("indrop.fastq.1.gz", 60)]
    subjects += ["", "A", "ACGT", "N" * 30]

    mismatches = 0
    compared = 0
    for pattern in (FIVE_PRIME, THREE_PRIME, INDROP):
        compiled = regex.compile(pattern)
        for subject in subjects:
            m = compiled.match(subject)
            py_spans = ({k: list(m.span(k)) for k in compiled.groupindex
                         if m.span(k) != (-1, -1)} if m else None)
            rs_spans = rs.match_pattern(pattern, subject)
            if rs_spans is not None:
                rs_spans = {k: list(v) for k, v in rs_spans.items()}
            compared += 1
            if py_spans != rs_spans:
                mismatches += 1
    rec.record("match_pattern", "mismatched_subjects", value=mismatches,
               threshold=0, tier=1, passed=mismatches == 0)
    rec.log("match_pattern", "%d pattern/subject pairs compared" % compared)


def main():
    rec = ResultRecorder("extract")
    compare_helpers(rec)
    compare_patterns(rec)
    for label, py_kwargs, rs_kwargs, fixture in CONFIGS:
        compare_config(rec, label, py_kwargs, rs_kwargs, fixture)
    rec.save()
    return 0


if __name__ == "__main__":
    sys.exit(main())
