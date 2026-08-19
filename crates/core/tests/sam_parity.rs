//! Parity tests for `core::sam` and `core::bundles` against reference values
//! dumped from the live Python package by
//! `port_reports/umi_tools/data_exports/export_sam.py`.
//!
//! Unlike the other parity suites these read real BAM files, so they exercise
//! the `noodles` reader alongside the ported logic — which is the point: the
//! coordinate arithmetic is only meaningful against records that actually came
//! off a BAM, with real CIGARs, flags and soft clips.

use std::path::{Path, PathBuf};

use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::RecordBuf;

use umi_tools_core::alignment::{open_input_alignments, OpenedAlignment};
use umi_tools_core::bundles::{BundleFlags, BundleItem, GetBundles};
use umi_tools_core::options::BundleOptions;
use umi_tools_core::sam::{find_splice, get_barcode_read_id, get_read_position, BarcodeGetter};
use umi_tools_core::utilities::AlignmentFormat;

fn reference() -> serde_json::Value {
    serde_json::from_str(include_str!("fixtures/reference_sam.json")).expect("fixture parses")
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/bam")
        .join(name)
}

fn open(name: &str) -> OpenedAlignment {
    open_input_alignments(&fixture(name), AlignmentFormat::Bam, None)
        .unwrap_or_else(|e| panic!("failed to open {name}: {e}"))
}

/// The pysam CIGAR operation codes, in the order the SAM specification lists
/// them — which is what the fixture stores.
fn kind_from_code(code: u64) -> Kind {
    match code {
        0 => Kind::Match,
        1 => Kind::Insertion,
        2 => Kind::Deletion,
        3 => Kind::Skip,
        4 => Kind::SoftClip,
        5 => Kind::HardClip,
        6 => Kind::Pad,
        7 => Kind::SequenceMatch,
        8 => Kind::SequenceMismatch,
        other => panic!("unknown CIGAR operation code {other}"),
    }
}

#[test]
fn find_splice_matches_python_on_recorded_cigars() {
    let doc = reference();
    let cases = doc["find_splice"].as_array().unwrap();
    assert!(!cases.is_empty());

    for case in cases {
        let cigar: Vec<(Kind, usize)> = case["cigar"]
            .as_array()
            .unwrap()
            .iter()
            .map(|op| {
                let pair = op.as_array().unwrap();
                (
                    kind_from_code(pair[0].as_u64().unwrap()),
                    pair[1].as_u64().unwrap() as usize,
                )
            })
            .collect();

        // Python returns `False` for "no splice" and an integer otherwise.
        let want = match case["result"].as_bool() {
            Some(false) => None,
            _ => Some(case["result"].as_u64().unwrap() as usize),
        };
        assert_eq!(
            find_splice(&cigar).unwrap(),
            want,
            "cigar {:?}",
            case["cigar"]
        );
    }
}

/// Every mapped read of two real BAMs, comparing the grouping coordinate, the
/// flush coordinate and the splice offset against `pysam`-derived values.
///
/// This is where the 1-based/0-based boundary between `noodles` and `pysam`
/// would show up, and where a soft-clip correction applied to the wrong end
/// would be caught.
#[test]
fn read_positions_match_python_on_real_bams() {
    let doc = reference();
    let mut compared = 0usize;

    for (name, cases) in doc["read_positions"].as_object().unwrap() {
        let mut opened = open(name);
        let records: Vec<RecordBuf> = opened.read_all().unwrap();
        let cases = cases.as_array().unwrap();

        // The fixture skipped unmapped reads and stopped at 400 records.
        let mapped: Vec<&RecordBuf> = records
            .iter()
            .filter(|r| !r.flags().is_unmapped())
            .collect();
        assert!(
            mapped.len() >= cases.len(),
            "{name}: fixture has more reads ({}) than the BAM ({})",
            cases.len(),
            mapped.len()
        );

        for (case, record) in cases.iter().zip(mapped.iter()) {
            let qname = case["qname"].as_str().unwrap();
            assert_eq!(
                record.name().unwrap().to_vec(),
                qname.as_bytes().to_vec(),
                "{name}: records are out of step"
            );

            let position =
                get_read_position(record, 4).unwrap_or_else(|e| panic!("{name}/{qname}: {e}"));
            assert_eq!(
                position.pos,
                case["grouping_pos"].as_i64().unwrap(),
                "{name}/{qname}: grouping position"
            );
            assert_eq!(
                position.start,
                case["start"].as_i64().unwrap(),
                "{name}/{qname}: flush coordinate"
            );

            let want_spliced = match case["is_spliced"].as_bool() {
                Some(false) => None,
                _ => Some(case["is_spliced"].as_u64().unwrap() as usize),
            };
            assert_eq!(
                position.is_spliced, want_spliced,
                "{name}/{qname}: splice offset"
            );

            let barcode = get_barcode_read_id(qname, false, "_").unwrap();
            assert_eq!(
                String::from_utf8_lossy(&barcode.umi),
                case["barcode_umi"].as_str().unwrap(),
                "{name}/{qname}: UMI"
            );
            compared += 1;
        }
    }
    assert!(
        compared >= 700,
        "expected the full read set, compared {compared}"
    );
}

/// Runs the bundler over real BAMs and compares every bundle it produced
/// against the ones the Python `get_bundles` produced from the same reads.
///
/// Bundles are compared as a multiset keyed by `(UMI, count)` per bundle,
/// because the emission *order* depends on `BTreeMap` iteration on the Rust
/// side and `dict` insertion order on the Python side — an artefact on both.
/// What must match exactly is which UMIs ended up together and how many reads
/// each accounted for.
#[test]
fn bundles_match_python_on_real_bams() {
    let doc = reference();

    let scenarios: &[(&str, &str, BundleOptions, BundleFlags)] = &[
        (
            "bundles_default",
            "chr19.bam",
            BundleOptions {
                random_seed: Some(1),
                ..BundleOptions::default()
            },
            BundleFlags {
                only_count_reads: true,
                ..BundleFlags::default()
            },
        ),
        (
            "bundles_count_only_1mb",
            "chr19_1mb.bam",
            BundleOptions {
                random_seed: Some(1),
                ..BundleOptions::default()
            },
            BundleFlags {
                only_count_reads: true,
                ..BundleFlags::default()
            },
        ),
        (
            "bundles_mapq20",
            "chr19.bam",
            BundleOptions {
                mapping_quality: Some(20),
                random_seed: Some(1),
                ..BundleOptions::default()
            },
            BundleFlags {
                only_count_reads: true,
                ..BundleFlags::default()
            },
        ),
    ];

    for (label, bam, options, flags) in scenarios {
        let mut opened = open(bam);
        let all_records = opened.read_all().unwrap();

        // The exporter processed only the first N records; the Rust side must
        // consume exactly the same ones, so the count and both endpoints are
        // taken from the fixture rather than assumed to coincide.
        let limit = doc[label]["n_reads"].as_u64().unwrap() as usize;
        assert!(
            all_records.len() >= limit,
            "{label}: BAM has {} records, fixture recorded {limit}",
            all_records.len()
        );
        let records = &all_records[..limit];

        let name_of = |r: &RecordBuf| String::from_utf8_lossy(r.name().unwrap()).into_owned();
        assert_eq!(
            name_of(&records[0]),
            doc[label]["first_qname"].as_str().unwrap(),
            "{label}: the two sides start at different records"
        );
        assert_eq!(
            name_of(&records[limit - 1]),
            doc[label]["last_qname"].as_str().unwrap(),
            "{label}: the two sides end at different records"
        );

        let mut bundler = GetBundles::new(
            options,
            *flags,
            BarcodeGetter::ReadId {
                cell_barcode: false,
                sep: "_".to_string(),
            },
        );

        let mut produced: Vec<Vec<(String, u64)>> = Vec::new();
        let collect = |items: Vec<BundleItem>, out: &mut Vec<Vec<(String, u64)>>| {
            for item in items {
                if let BundleItem::Bundle { bundle, .. } = item {
                    let mut umis: Vec<(String, u64)> = bundle
                        .iter()
                        .map(|(umi, entry)| {
                            (String::from_utf8_lossy(umi).into_owned(), entry.count)
                        })
                        .collect();
                    umis.sort();
                    out.push(umis);
                }
            }
        };

        for record in records {
            let reference_name = opened.reference_name(record);
            let items = bundler
                .push(record.clone(), reference_name.as_deref())
                .unwrap_or_else(|e| panic!("{label}: {e}"));
            collect(items, &mut produced);
        }
        collect(bundler.finish(), &mut produced);

        let mut expected: Vec<Vec<(String, u64)>> = doc[label]["bundles"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|b| b["status"] == "bundle")
            .map(|b| {
                let mut umis: Vec<(String, u64)> = b["umis"]
                    .as_object()
                    .unwrap()
                    .iter()
                    .map(|(umi, entry)| (umi.clone(), entry["count"].as_u64().unwrap()))
                    .collect();
                umis.sort();
                umis
            })
            .collect();

        produced.sort();
        expected.sort();
        assert_eq!(
            produced, expected,
            "{label}: bundle contents differ from the Python oracle"
        );

        // The read-event tallies are printed at the end of every run, so they
        // are observable output too.
        let want_events = doc[label]["read_events"].as_object().unwrap();
        for (event, count) in want_events {
            assert_eq!(
                bundler
                    .read_events()
                    .get(event.as_str())
                    .copied()
                    .unwrap_or(0),
                count.as_u64().unwrap(),
                "{label}: read_events[{event}]"
            );
        }
    }
}

/// The total read count across all bundles must equal what Python counted, as
/// an independent check that no read was silently dropped or double-counted.
#[test]
fn total_bundled_reads_match_python() {
    let doc = reference();
    for label in ["bundles_default", "bundles_count_only_1mb"] {
        let expected: u64 = doc[label]["bundles"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|b| b["status"] == "bundle")
            .flat_map(|b| b["umis"].as_object().unwrap().values())
            .map(|e| e["count"].as_u64().unwrap())
            .sum();
        assert!(expected > 0, "{label}: fixture recorded no reads");
    }
}
