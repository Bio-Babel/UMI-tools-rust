"""Export reference values for core::sam and core::bundles.

Runs the live sam_methods against the checked-in BAM fixtures and records the
per-read coordinate arithmetic plus the full bundle structure, so a Rust
divergence can be attributed to a single read rather than to the aggregate.
"""
import json
import os

import pysam

import umi_tools.sam_methods as sm

OUT = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(os.environ.get("UMI_TOOLS_SRC", "../umi_tools_src"), "tests")


class Options:
    """A stand-in for the optparse Values object get_bundles reads."""

    def __init__(self, **kwargs):
        defaults = dict(
            get_umi_method="read_id", umi_sep="_", per_cell=False, per_gene=False,
            per_contig=False, gene_tag=None, assigned_tag=None, skip_regex=None,
            paired=False, spliced=False, ignore_tlen=False, ignore_umi=False,
            read_length=False, soft_clip_threshold=4, subset=None,
            mapping_quality=None, detection_method=None, whole_contig=False,
            unmapped_reads="discard", unpaired_reads="discard",
            chimeric_pairs="discard", umi_tag="RX", cell_tag=None,
            umi_tag_split=None, umi_tag_delim=None, cell_tag_split="-",
            cell_tag_delim=None,
        )
        defaults.update(kwargs)
        for k, v in defaults.items():
            setattr(self, k, v)


def dump_read_positions(bam, limit=400):
    rows = []
    f = pysam.AlignmentFile(os.path.join(TESTS, bam))
    for n, read in enumerate(f.fetch()):
        if n >= limit:
            break
        if read.is_unmapped:
            continue
        start, pos, spliced = sm.get_read_position(read, 4)
        rows.append({
            "qname": read.qname,
            "reference_name": read.reference_name,
            "flag": read.flag,
            "pos": read.pos,
            "aend": read.aend,
            "cigar": [[op, ln] for op, ln in read.cigar],
            "is_reverse": bool(read.is_reverse),
            "mapq": read.mapping_quality,
            "start": int(start),
            "grouping_pos": int(pos),
            "is_spliced": (False if spliced is False else int(spliced)),
            "barcode_umi": sm.get_barcode_read_id(read)[0].decode(),
        })
    f.close()
    return rows


def dump_bundles(bam, options, flags, limit=4000):
    """Run the real bundler and record every bundle it produced."""
    f = pysam.AlignmentFile(os.path.join(TESTS, bam))
    bundler = sm.get_bundles(options, **flags)
    rows = []
    reads = []
    for n, r in enumerate(f.fetch()):
        if n >= limit:
            break
        reads.append(r)
    for bundle, key, status in bundler(iter(reads)):
        if status == "single_read":
            rows.append({"status": "single_read",
                         "qname": bundle.qname})
            continue
        umis = {}
        for umi in bundle:
            entry = {"count": bundle[umi]["count"]}
            if "read" in bundle[umi]:
                read = bundle[umi]["read"]
                entry["read"] = (read.qname if not isinstance(read, list)
                                 else sorted(x.qname for x in read))
            umis[umi.decode()] = entry
        # `key` is ((coordinate-tuple or gene), cell)
        inner, cell = key
        rows.append({
            "status": "bundle",
            "key_inner": (list(inner) if isinstance(inner, tuple) else inner),
            "key_cell": (cell.decode() if isinstance(cell, bytes) else cell),
            "umis": umis,
        })
    f.close()
    return {"bundles": rows,
            "read_events": dict(bundler.read_events),
            # The Rust side must consume exactly these records, so the
            # truncation point and its endpoints are recorded rather than
            # left to coincide.
            "n_reads": len(reads),
            "first_qname": reads[0].qname if reads else None,
            "last_qname": reads[-1].qname if reads else None}


def main():
    payload = {
        "read_positions": {
            "chr19.bam": dump_read_positions("chr19.bam"),
            "chr19_1mb.bam": dump_read_positions("chr19_1mb.bam"),
        },
        "find_splice": [
            {"cigar": cg, "result": (False if (r := sm.find_splice(cg)) is False
                                     else int(r))}
            for cg in [[(0, 50)], [(4, 3), (0, 47)], [(0, 20), (3, 100), (0, 30)],
                       [(0, 20), (4, 10)], [(5, 2), (0, 48)],
                       [(4, 5), (0, 20), (3, 100)],
                       [(0, 10), (1, 5), (2, 3), (3, 100)]]
        ],
        "bundles_default": dump_bundles(
            "chr19.bam", Options(), dict(only_count_reads=True)),
        "bundles_count_only_1mb": dump_bundles(
            "chr19_1mb.bam", Options(), dict(only_count_reads=True)),
        "bundles_mapq20": dump_bundles(
            "chr19.bam", Options(mapping_quality=20), dict(only_count_reads=True)),
        "bundles_all_reads": dump_bundles(
            "chr19_1mb.bam", Options(), dict(all_reads=True, return_read2=True)),
    }
    path = os.path.join(OUT, "reference_sam.json")
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=1)
    n_pos = sum(len(v) for v in payload["read_positions"].values())
    n_bundles = sum(len(payload[k]["bundles"]) for k in payload
                    if k.startswith("bundles_"))
    print("wrote %s (%d read positions, %d bundles)" % (path, n_pos, n_bundles))


if __name__ == "__main__":
    main()
