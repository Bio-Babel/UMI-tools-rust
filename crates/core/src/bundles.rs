//! Grouping alignment records into bundles of candidate duplicates
//! (Python source: `sam_methods.get_bundles`).
//!
//! A "bundle" is every read sharing an alignment position (or a gene) and a
//! cell barcode, keyed by UMI. `dedup`, `group` and `count` all consume
//! bundles: deduplication is deciding which UMIs within one bundle came from
//! the same original molecule.
//!
//! Bundles are emitted **as the input is consumed**, not at the end: once the
//! scan has moved 1000 bases past a position (or onto a new contig, or a new
//! gene), nothing further can join that position's bundle, so it is flushed.
//! That is what lets `umi_tools` deduplicate a BAM far larger than memory.

use std::collections::BTreeMap;

use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::RecordBuf;

use crate::error::{PortError, Result};
use crate::options::{BundleOptions, ChimericPairs, DetectionMethod, UnmappedReads, UnpairedReads};
use crate::sam::{get_read_position, BarcodeGetter};
use crate::utilities;

/// The reads a bundle entry retained.
///
/// Python distinguishes these three modes only by which dictionary keys it
/// happens to populate, so reading `bundle[umi]["read"]` in `only_count_reads`
/// mode is a `KeyError`. The enum makes the distinction a compile-time one.
#[derive(Debug, Clone, PartialEq)]
pub enum BundleReads {
    /// The single best read for this key (the default for `dedup`).
    One(Box<RecordBuf>),
    /// Every read for this key (`group`, via `all_reads`).
    Many(Vec<RecordBuf>),
    /// No reads retained, only the count (`count`, via `only_count_reads`).
    CountOnly,
}

impl BundleReads {
    /// The single retained read, if this entry holds one.
    #[must_use]
    pub fn one(&self) -> Option<&RecordBuf> {
        match self {
            BundleReads::One(r) => Some(r),
            _ => None,
        }
    }

    /// Every retained read, empty when none were kept.
    #[must_use]
    pub fn all(&self) -> &[RecordBuf] {
        match self {
            BundleReads::Many(v) => v,
            _ => &[],
        }
    }
}

/// One UMI's entry in a bundle
/// (Python source: the `{"read": ..., "count": ...}` inner dictionary).
#[derive(Debug, Clone, PartialEq)]
pub struct BundleEntry {
    /// The retained read(s).
    pub reads: BundleReads,
    /// How many reads shared this key.
    pub count: u64,
    /// A monotonic sequence number recording when this UMI was first seen.
    ///
    /// Python's bundle is a `dict`, which preserves insertion order, and the
    /// clustering stage takes its UMIs in that order — which decides the order
    /// the deduplicated reads are written in. A `BTreeMap` sorts instead, so
    /// the original order has to be carried explicitly; see
    /// [`insertion_order`].
    pub first_seen: u64,
}

/// UMI to entry, for one position/gene and cell barcode
/// (Python source: the dictionary `get_bundles` yields).
pub type Bundle = BTreeMap<Vec<u8>, BundleEntry>;

/// The key a bundle was grouped under
/// (Python source: the `key` element of the yielded 3-tuple).
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct BundleKey {
    /// The gene, when `--per-gene` is in effect.
    pub gene: Option<String>,
    /// Reverse strand, splice state, template length and read length, when
    /// grouping by coordinate.
    pub coordinate: Option<CoordinateKey>,
    /// The cell barcode, when `--per-cell` is in effect.
    pub cell: Option<Vec<u8>>,
}

/// The coordinate-derived part of a bundle key
/// (Python source: the `(is_reverse, spliced, tlen, r_length)` tuple).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct CoordinateKey {
    /// Whether the read is on the reverse strand.
    pub is_reverse: bool,
    /// Whether the read is spliced, when `--spliced-is-unique` is set.
    pub spliced: bool,
    /// The template length, zeroed unless paired and `--ignore-tlen` is off.
    pub template_length: i32,
    /// The read length, zero unless `--read-length` is set.
    pub read_length: usize,
}

/// What the bundler produced for one input record or position
/// (Python source: the `(bundle_or_read, key, status)` 3-tuple, whose `status`
/// is either `"single_read"` or `"bundle"`).
#[derive(Debug, Clone)]
pub enum BundleItem {
    /// A read passed straight through without grouping — a read 2, an unmapped
    /// read, or a chimeric/unpaired read the options said to output.
    Single(Box<RecordBuf>),
    /// A completed bundle.
    Bundle {
        /// UMI to entry.
        bundle: Bundle,
        /// What the bundle was grouped under.
        key: BundleKey,
    },
}

/// The bundle's UMIs in the order they were first seen
/// (Python source: the `dict` insertion order `bundle.keys()` yields).
///
/// `dedup`, `group` and `count` all pass this to
/// [`crate::network::UmiClusterer::cluster_ordered`], because the clustering
/// order decides the order the output reads are written in and the package's
/// own test suite compares that output byte for byte.
///
/// # Examples
///
/// ```
/// use umi_tools_core::bundles::{insertion_order, Bundle, BundleEntry, BundleReads};
///
/// let mut bundle = Bundle::new();
/// // "TTTT" sorts after "AAAA" but was seen first.
/// bundle.insert(b"TTTT".to_vec(), BundleEntry {
///     reads: BundleReads::CountOnly, count: 1, first_seen: 0 });
/// bundle.insert(b"AAAA".to_vec(), BundleEntry {
///     reads: BundleReads::CountOnly, count: 1, first_seen: 1 });
///
/// assert_eq!(insertion_order(&bundle), vec![b"TTTT".to_vec(), b"AAAA".to_vec()]);
/// ```
#[must_use]
pub fn insertion_order(bundle: &Bundle) -> Vec<Vec<u8>> {
    let mut entries: Vec<(&Vec<u8>, u64)> = bundle
        .iter()
        .map(|(umi, entry)| (umi, entry.first_seen))
        .collect();
    // Stable, and ties fall back to the UMI itself, so the order is total
    // even if two entries somehow share a sequence number.
    entries.sort_by(|a, b| a.1.cmp(&b.1).then_with(|| a.0.cmp(b.0)));
    entries.into_iter().map(|(umi, _)| umi.clone()).collect()
}

/// Which reads the bundler should retain and pass through
/// (Python source: the five boolean constructor flags of `get_bundles`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct BundleFlags {
    /// Retain counts only, no reads (`count`).
    pub only_count_reads: bool,
    /// Retain every read per key (`group`).
    pub all_reads: bool,
    /// Emit unmapped reads instead of dropping them.
    pub return_unmapped: bool,
    /// Emit read 2 immediately instead of dropping it.
    pub return_read2: bool,
}

/// How far past a position the scan must travel before that position's bundle
/// is flushed (Python source: the literal `1000` in `check_output`).
///
/// It is a fixed margin, not derived from the read length, so a library with
/// reads longer than 1000 bases could in principle flush a position early.
/// Mirrored rather than corrected.
pub const FLUSH_MARGIN: i64 = 1000;

/// Groups alignment records into bundles of candidate duplicates
/// (Python source: `sam_methods.get_bundles`).
///
/// Feed records in coordinate order with [`push`](Self::push), collecting the
/// items it returns, then call [`finish`](Self::finish) to flush what is left.
/// Python models this as a generator; a push-based type is the natural Rust
/// shape and lets the caller own the record source, which matters because the
/// source may be a BAM reader, a CRAM reader or a metacontig fetcher.
///
/// # Examples
///
/// ```no_run
/// use umi_tools_core::bundles::{BundleFlags, GetBundles};
/// use umi_tools_core::options::BundleOptions;
/// use umi_tools_core::sam::BarcodeGetter;
///
/// let options = BundleOptions::default();
/// let getter = BarcodeGetter::ReadId { cell_barcode: false, sep: "_".into() };
/// let mut bundler = GetBundles::new(&options, BundleFlags::default(), getter);
///
/// # let records: Vec<noodles_sam::alignment::RecordBuf> = Vec::new();
/// for record in records {
///     for item in bundler.push(record, None)? {
///         // handle each Single or Bundle
///         let _ = item;
///     }
/// }
/// let _remaining = bundler.finish();
/// # Ok::<(), umi_tools_core::PortError>(())
/// ```
#[derive(Debug)]
pub struct GetBundles {
    options: BundleOptions,
    flags: BundleFlags,
    barcode_getter: BarcodeGetter,

    /// Position (or gene) to key to UMI to entry.
    reads_dict: BTreeMap<PositionKey, BTreeMap<BundleKey, Bundle>>,
    /// Reservoir-sampling counters, parallel to `reads_dict`.
    read_counts: BTreeMap<PositionKey, BTreeMap<BundleKey, BTreeMap<Vec<u8>, u64>>>,
    /// Gene to the transcripts observed for it, for the metacontig path.
    observed_contigs: BTreeMap<String, std::collections::BTreeSet<String>>,
    metacontig_contig: Option<crate::sam::MetaContigMap>,
    contig_metacontig: BTreeMap<String, String>,

    read_events: BTreeMap<String, u64>,
    last_pos: PositionKey,
    last_chr: Option<String>,
    start: i64,
    current_chr: Option<String>,
    /// Increments once per newly-seen UMI, recording insertion order.
    sequence: u64,

    rng: crate::mersenne::MersenneTwister,
}

/// What a bundle is filed under before its cell barcode is applied: either a
/// gene name or a coordinate.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum PositionKey {
    /// A gene, under `--per-gene`.
    Gene(String),
    /// A coordinate on the current contig.
    Coordinate(i64),
}

impl Default for PositionKey {
    fn default() -> Self {
        PositionKey::Coordinate(0)
    }
}

/// What should happen to a newly-seen read at an already-occupied key.
///
/// Computed with only `&self` borrowed, so the mutable borrow of
/// `reads_dict` in [`update_dicts`](Self::update_dicts) never overlaps a
/// read of `self.options` or a call to the sampler. That separation is the
/// point: the alternative — cloning the stored record so the comparison
/// can happen under a mutable borrow — would copy a record per duplicate
/// read, on the hottest path in the program.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Replacement {
    /// Keep the stored read.
    Keep,
    /// Replace it, and reset the tie counter.
    Replace,
    /// Still tied: fall through to reservoir sampling.
    Tied,
}

impl GetBundles {
    /// Creates a bundler.
    ///
    /// `seed` for the tie-breaking sampler comes from
    /// [`BundleOptions::random_seed`]; see
    /// [`update_dicts`](Self::update_dicts) for why that matters.
    #[must_use]
    pub fn new(options: &BundleOptions, flags: BundleFlags, barcode_getter: BarcodeGetter) -> Self {
        Self {
            options: options.clone(),
            flags,
            barcode_getter,
            reads_dict: BTreeMap::new(),
            read_counts: BTreeMap::new(),
            observed_contigs: BTreeMap::new(),
            metacontig_contig: None,
            contig_metacontig: BTreeMap::new(),
            read_events: BTreeMap::new(),
            last_pos: PositionKey::default(),
            last_chr: None,
            start: 0,
            current_chr: None,
            sequence: 0,
            // `Utilities.Start` calls `random.seed(options.random_seed)` when
            // the flag is given, so the tie-break stream is CPython's. Without
            // the flag, Python leaves `random` auto-seeded from OS entropy and
            // the choice is genuinely arbitrary; seeding from zero here makes
            // that case reproducible without changing the distribution.
            rng: crate::mersenne::MersenneTwister::seed_from_u64(options.random_seed.unwrap_or(0)),
        }
    }

    /// Attaches a gene-to-transcript map, enabling the metacontig path
    /// (Python source: the `metacontig_contig` constructor argument, which
    /// also builds the inverse `contig_metacontig` map).
    pub fn with_metacontigs(&mut self, map: crate::sam::MetaContigMap) {
        self.contig_metacontig.clear();
        for (metacontig, contigs) in &map {
            for contig in contigs {
                self.contig_metacontig
                    .insert(contig.clone(), metacontig.clone());
            }
        }
        self.metacontig_contig = Some(map);
    }

    /// The event tallies `dedup`, `group` and `count` print at the end of a run
    /// (Python source: `get_bundles.read_events`).
    #[must_use]
    pub fn read_events(&self) -> &BTreeMap<String, u64> {
        &self.read_events
    }

    /// The barcode getter this bundler was built with
    /// (Python source: `get_bundles.barcode_getter`, which `dedup` reuses to
    /// re-read UMIs when building its stats).
    #[must_use]
    pub fn barcode_getter(&self) -> &BarcodeGetter {
        &self.barcode_getter
    }

    fn event(&mut self, name: &str) {
        *self.read_events.entry(name.to_string()).or_insert(0) += 1;
    }

    /// Returns the next uniform value in `[0, 1)`, from `CPython`'s generator.
    ///
    /// The tie-break between equally good reads consumes this stream, and the
    /// package's own test suite compares deduplicated output byte for byte
    /// against references produced with `--random-seed=123456789`. Matching
    /// them means matching `CPython`'s Mersenne Twister exactly, not merely its
    /// distribution — see [`crate::mersenne`].
    fn next_uniform(&mut self) -> f64 {
        self.rng.random()
    }

    fn multimap_tag(record: &RecordBuf, tag: [u8; 2]) -> Option<i64> {
        use noodles_sam::alignment::record_buf::data::field::Value;
        match record.data().get(&Tag::from(tag))? {
            Value::Int8(v) => Some(i64::from(*v)),
            Value::UInt8(v) => Some(i64::from(*v)),
            Value::Int16(v) => Some(i64::from(*v)),
            Value::UInt16(v) => Some(i64::from(*v)),
            Value::Int32(v) => Some(i64::from(*v)),
            Value::UInt32(v) => Some(i64::from(*v)),
            _ => None,
        }
    }

    fn xt_is_unique(record: &RecordBuf) -> bool {
        use noodles_sam::alignment::record_buf::data::field::Value;
        matches!(
            record.data().get(&Tag::from(*b"XT")),
            Some(Value::Character(b'U'))
        )
    }

    /// Decides between the stored read and a new one, on mapping quality and
    /// then the multimapping tag (Python source: the comparison chain in the
    /// `else` branch of `get_bundles.update_dicts`).
    fn compare(&self, current: &RecordBuf, candidate: &RecordBuf) -> Replacement {
        let current_mapq = crate::sam::raw_mapping_quality(current);
        let new_mapq = crate::sam::raw_mapping_quality(candidate);
        if current_mapq > new_mapq {
            return Replacement::Keep;
        }
        if current_mapq < new_mapq {
            return Replacement::Replace;
        }

        match self.options.detection_method {
            Some(method @ (DetectionMethod::Nh | DetectionMethod::X0)) => {
                let tag = if method == DetectionMethod::Nh {
                    *b"NH"
                } else {
                    *b"X0"
                };
                match (
                    Self::multimap_tag(current, tag),
                    Self::multimap_tag(candidate, tag),
                ) {
                    // Fewer reported alignments is the better read.
                    (Some(a), Some(b)) if a < b => Replacement::Keep,
                    (Some(a), Some(b)) if a > b => Replacement::Replace,
                    _ => Replacement::Tied,
                }
            }
            Some(DetectionMethod::Xt) => {
                if Self::xt_is_unique(current) {
                    Replacement::Keep
                } else if Self::xt_is_unique(candidate) {
                    Replacement::Replace
                } else {
                    Replacement::Tied
                }
            }
            None => Replacement::Tied,
        }
    }

    fn reset_tie_counter(&mut self, pos: &PositionKey, key: &BundleKey, umi: &[u8]) {
        self.read_counts
            .entry(pos.clone())
            .or_default()
            .entry(key.clone())
            .or_default()
            .insert(umi.to_vec(), 0);
    }

    /// Files one read into the bundle at `(pos, key, umi)`
    /// (Python source: `get_bundles.update_dicts`).
    ///
    /// Which read is retained per key depends on the mode:
    ///
    /// - `all_reads`: every read is kept.
    /// - `only_count_reads`: none are; only the count matters.
    /// - otherwise the "best" read wins, by mapping quality first, then by the
    ///   multimapping tag if one was chosen, and finally by reservoir sampling
    ///   among the reads still tied — so every tied read has an equal chance
    ///   rather than the first or last winning systematically.
    pub fn update_dicts(
        &mut self,
        record: RecordBuf,
        pos: PositionKey,
        key: BundleKey,
        umi: Vec<u8>,
    ) {
        // Reserved up front so every insertion arm sees the same value and no
        // arm has to reach back through a live borrow of `reads_dict`.
        let sequence = self.sequence;
        let claim = |bundler: &mut Self, umi: &[u8], pos: &PositionKey, key: &BundleKey| {
            let already_present = bundler
                .reads_dict
                .get(pos)
                .and_then(|by_key| by_key.get(key))
                .is_some_and(|bundle| bundle.contains_key(umi));
            if !already_present {
                bundler.sequence += 1;
            }
        };
        claim(self, &umi, &pos, &key);

        if self.flags.all_reads {
            self.reads_dict
                .entry(pos)
                .or_default()
                .entry(key)
                .or_default()
                .entry(umi)
                .and_modify(|e| {
                    e.count += 1;
                    if let BundleReads::Many(v) = &mut e.reads {
                        v.push(record.clone());
                    }
                })
                .or_insert_with(|| BundleEntry {
                    reads: BundleReads::Many(vec![record]),
                    count: 1,
                    first_seen: sequence,
                });
            return;
        }

        if self.flags.only_count_reads {
            self.reads_dict
                .entry(pos)
                .or_default()
                .entry(key)
                .or_default()
                .entry(umi)
                .and_modify(|e| e.count += 1)
                .or_insert(BundleEntry {
                    reads: BundleReads::CountOnly,
                    count: 1,
                    first_seen: sequence,
                });
            return;
        }

        // Single-best mode. Decide first, with only `&self` borrowed, then
        // apply — see `Replacement`.
        let existing = self
            .reads_dict
            .get(&pos)
            .and_then(|by_key| by_key.get(&key))
            .and_then(|bundle| bundle.get(&umi));

        let Some(existing) = existing else {
            self.reads_dict
                .entry(pos.clone())
                .or_default()
                .entry(key.clone())
                .or_default()
                .insert(
                    umi.clone(),
                    BundleEntry {
                        reads: BundleReads::One(Box::new(record)),
                        count: 1,
                        first_seen: sequence,
                    },
                );
            self.reset_tie_counter(&pos, &key, &umi);
            return;
        };

        let decision = match existing.reads.one() {
            Some(current) => self.compare(current, &record),
            // `only_count_reads` never runs this branch, and the other two
            // modes returned above, so this is unreachable in practice.
            None => Replacement::Keep,
        };

        // The count rises however the comparison went.
        if let Some(entry) = self
            .reads_dict
            .get_mut(&pos)
            .and_then(|by_key| by_key.get_mut(&key))
            .and_then(|bundle| bundle.get_mut(&umi))
        {
            entry.count += 1;
        }

        match decision {
            Replacement::Keep => {}
            Replacement::Replace => {
                if let Some(entry) = self
                    .reads_dict
                    .get_mut(&pos)
                    .and_then(|by_key| by_key.get_mut(&key))
                    .and_then(|bundle| bundle.get_mut(&umi))
                {
                    entry.reads = BundleReads::One(Box::new(record));
                }
                self.reset_tie_counter(&pos, &key, &umi);
            }
            Replacement::Tied => {
                let seen = self
                    .read_counts
                    .entry(pos.clone())
                    .or_default()
                    .entry(key.clone())
                    .or_default()
                    .entry(umi.clone())
                    .or_insert(0);
                *seen += 1;
                let probability = 1.0 / (*seen as f64);
                if self.next_uniform() < probability {
                    if let Some(entry) = self
                        .reads_dict
                        .get_mut(&pos)
                        .and_then(|by_key| by_key.get_mut(&key))
                        .and_then(|bundle| bundle.get_mut(&umi))
                    {
                        entry.reads = BundleReads::One(Box::new(record));
                    }
                }
            }
        }
    }

    /// Whether any bundles can be flushed, and which
    /// (Python source: `get_bundles.check_output`).
    ///
    /// Three regimes, matching Python's branches:
    ///
    /// - `--per-gene`: flush on a contig change; with a metacontig map, only
    ///   once every transcript of the last gene has been seen.
    /// - `--whole-contig`: flush on a contig change.
    /// - otherwise: flush once the scan is more than [`FLUSH_MARGIN`] bases
    ///   past a position, or on a contig change.
    #[must_use]
    pub fn check_output(&self) -> (bool, Vec<PositionKey>) {
        let chr_changed = self.current_chr != self.last_chr;

        if self.options.per_gene {
            if self.metacontig_contig.is_some() {
                let complete = match (&self.last_pos, &self.metacontig_contig) {
                    (PositionKey::Gene(gene), Some(map)) => {
                        self.observed_contigs.get(gene) == map.get(gene)
                    }
                    _ => false,
                };
                if chr_changed && complete {
                    return (true, vec![self.last_pos.clone()]);
                }
                return (false, Vec::new());
            }
            if chr_changed {
                return (true, self.reads_dict.keys().cloned().collect());
            }
            return (false, Vec::new());
        }

        if self.options.whole_contig {
            if chr_changed {
                return (true, self.reads_dict.keys().cloned().collect());
            }
            return (false, Vec::new());
        }

        if self.start > self.last_pos_coordinate() + FLUSH_MARGIN || chr_changed {
            let mut keys: Vec<PositionKey> = self.reads_dict.keys().cloned().collect();
            if !chr_changed {
                let cutoff = self.start - FLUSH_MARGIN;
                keys.retain(|k| matches!(k, PositionKey::Coordinate(c) if *c <= cutoff));
            }
            return (true, keys);
        }
        (false, Vec::new())
    }

    fn last_pos_coordinate(&self) -> i64 {
        match self.last_pos {
            PositionKey::Coordinate(c) => c,
            PositionKey::Gene(_) => 0,
        }
    }

    fn flush(&mut self, keys: &[PositionKey]) -> Vec<BundleItem> {
        let mut out = Vec::new();
        for position in keys {
            if let Some(by_key) = self.reads_dict.remove(position) {
                for (key, bundle) in by_key {
                    out.push(BundleItem::Bundle { bundle, key });
                }
            }
            self.read_counts.remove(position);
        }
        out
    }

    /// Feeds one record to the bundler
    /// (Python source: the body of `get_bundles.__call__`'s loop).
    ///
    /// `reference_name` is the record's contig; it is passed separately
    /// because `RecordBuf` stores a reference *id*, and resolving it needs the
    /// header the caller already holds.
    ///
    /// Returns the bundles (and pass-through reads) this record completed.
    ///
    /// # Errors
    ///
    /// - [`PortError::Value`] from the coordinate arithmetic or a malformed
    ///   read name.
    /// - [`PortError::Alignment`] for a mapped record with no position.
    pub fn push(
        &mut self,
        record: RecordBuf,
        reference_name: Option<&str>,
    ) -> Result<Vec<BundleItem>> {
        let flags = record.flags();

        if flags.is_segmented() && flags.is_last_segment() {
            if self.flags.return_read2 && (!flags.is_unmapped() || self.flags.return_unmapped) {
                return Ok(vec![BundleItem::Single(Box::new(record))]);
            }
            return Ok(Vec::new());
        }
        self.event("Input Reads");

        // From here on, only read 1s (or unpaired reads).
        if self.options.paired {
            if flags.is_segmented() {
                self.event("Read pairs");
            } else {
                self.event("Unpaired reads");
                match self.options.unpaired_reads {
                    UnpairedReads::Discard => return Ok(Vec::new()),
                    UnpairedReads::Output => return Ok(vec![BundleItem::Single(Box::new(record))]),
                    // "use" falls through: the pair is used and TLEN will be 0.
                    UnpairedReads::Use => {}
                }
            }
        }

        if flags.is_unmapped() {
            if self.options.paired {
                if flags.is_mate_unmapped() {
                    self.event("Both unmapped");
                } else {
                    self.event("Read 1 unmapped");
                }
            } else {
                self.event("Single end unmapped");
            }
            if self.flags.return_unmapped {
                // Python increments 'Input Reads' a second time here. It looks
                // like a slip, but it is observable in the log line the test
                // suite compares, so it is mirrored.
                self.event("Input Reads");
                return Ok(vec![BundleItem::Single(Box::new(record))]);
            }
            return Ok(Vec::new());
        }

        if self.options.paired && flags.is_mate_unmapped() {
            self.event("Read 2 unmapped");
            if self.options.unmapped_reads != UnmappedReads::Use {
                if self.flags.return_unmapped {
                    return Ok(vec![BundleItem::Single(Box::new(record))]);
                }
                return Ok(Vec::new());
            }
        }

        if flags.is_segmented()
            && record.reference_sequence_id() != record.mate_reference_sequence_id()
        {
            self.event("Chimeric read pair");
            match self.options.chimeric_pairs {
                ChimericPairs::Discard => return Ok(Vec::new()),
                ChimericPairs::Output => return Ok(vec![BundleItem::Single(Box::new(record))]),
                ChimericPairs::Use => {}
            }
        }

        if let Some(subset) = self.options.subset {
            if self.next_uniform() >= subset {
                self.event("Randomly excluded");
                return Ok(Vec::new());
            }
        }

        if let Some(threshold) = self.options.mapping_quality {
            if u32::from(crate::sam::raw_mapping_quality(&record)) < threshold {
                self.event("< MAPQ threshold");
                return Ok(Vec::new());
            }
        }

        // Read the barcode.
        let barcode = if self.options.ignore_umi {
            let cell = if self.options.per_cell {
                self.barcode_getter.get(&record)?.cell
            } else {
                Some(Vec::new())
            };
            crate::sam::Barcode {
                umi: Vec::new(),
                cell,
            }
        } else if let Ok(b) = self.barcode_getter.get(&record) {
            b
        } else {
            let msg = "Read skipped, missing umi and/or cell tag";
            if !self.read_events.contains_key(msg) {
                let name = record
                    .name()
                    .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                    .unwrap_or_default();
                utilities::warn(&format!(
                    "At least one read is missing UMI and/or cell tag(s): {name}"
                ));
            }
            self.event(msg);
            return Ok(Vec::new());
        };

        self.current_chr = reference_name.map(str::to_string);

        let mut out = Vec::new();
        let (pos, key);

        if self.options.per_gene {
            let mut transcript = None;
            let gene = if self.options.per_contig {
                let contig = reference_name.unwrap_or_default().to_string();
                if self.metacontig_contig.is_some() {
                    transcript = Some(contig.clone());
                    self.contig_metacontig
                        .get(&contig)
                        .cloned()
                        .ok_or_else(|| {
                            PortError::value(format!("contig {contig} has no metacontig"))
                        })?
                } else {
                    contig
                }
            } else {
                let Some(gene_tag) = self.options.gene_tag else {
                    return Err(PortError::value(
                        "for per-gene applications, must supply --per-contig or --gene-tag",
                    ));
                };
                let assigned_tag = self.options.assigned_tag.unwrap_or(gene_tag);
                let (Some(assigned), Some(gene)): (Option<String>, Option<String>) = (
                    crate::sam::tag_text(&record, assigned_tag),
                    crate::sam::tag_text(&record, gene_tag),
                ) else {
                    self.event("Read skipped, no tag");
                    return Ok(Vec::new());
                };
                if gene.is_empty() {
                    let msg = "Read skipped - gene string is empty";
                    if !self.read_events.contains_key(msg) {
                        utilities::warn("Assigned gene is empty string. First such read");
                    }
                    self.event(msg);
                    return Ok(Vec::new());
                }
                if let Some(skip) = &self.options.skip_regex {
                    if skip.is_match(&assigned) {
                        self.event("Read skipped - assigned tag matches skip_regex");
                        return Ok(Vec::new());
                    }
                }
                gene
            };

            pos = PositionKey::Gene(gene.clone());
            key = BundleKey {
                gene: Some(gene.clone()),
                coordinate: None,
                cell: barcode.cell.clone(),
            };

            if self.last_chr.is_some() {
                let (do_output, keys) = self.check_output();
                if do_output {
                    out.extend(self.flush(&keys));
                }
            }
            self.last_chr = self.current_chr.clone();
            self.last_pos = pos.clone();

            if let Some(transcript) = transcript {
                self.observed_contigs
                    .entry(gene)
                    .or_default()
                    .insert(transcript);
            }
        } else {
            let position = get_read_position(&record, self.options.soft_clip_threshold)?;
            self.start = position.start;

            let (do_output, keys) = self.check_output();
            if do_output {
                out.extend(self.flush(&keys));
            }
            self.last_pos = PositionKey::Coordinate(position.start);
            self.last_chr = self.current_chr.clone();

            let read_length = if self.options.read_length {
                record.sequence().len()
            } else {
                0
            };
            let template_length = if self.options.ignore_tlen || !self.options.paired {
                0
            } else {
                record.template_length()
            };

            pos = PositionKey::Coordinate(position.pos);
            key = BundleKey {
                gene: None,
                coordinate: Some(CoordinateKey {
                    is_reverse: record.flags().is_reverse_complemented(),
                    spliced: self.options.spliced && position.is_spliced.is_some(),
                    template_length,
                    read_length,
                }),
                cell: barcode.cell.clone(),
            };
        }

        self.update_dicts(record, pos, key, barcode.umi);
        Ok(out)
    }

    /// Flushes every remaining bundle
    /// (Python source: the trailing `for p in sorted(self.reads_dict.keys())`
    /// loop after `get_bundles.__call__`'s main loop).
    #[must_use]
    pub fn finish(&mut self) -> Vec<BundleItem> {
        let keys: Vec<PositionKey> = self.reads_dict.keys().cloned().collect();
        self.flush(&keys)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_core::Position;
    use noodles_sam::alignment::record::cigar::op::{Kind, Op};
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record_buf::Cigar as CigarBuf;

    /// A mapped read. `RecordBuf::default()` sets `Flags::UNMAPPED`, so the
    /// flags must be cleared explicitly or every record takes the
    /// unmapped-read path.
    fn read(name: &str, start: usize, mapq: u8) -> RecordBuf {
        let mut r = RecordBuf::default();
        *r.flags_mut() = Flags::empty();
        *r.name_mut() = Some(name.into());
        *r.alignment_start_mut() = Position::new(start);
        *r.cigar_mut() = [Op::new(Kind::Match, 50)].into_iter().collect::<CigarBuf>();
        *r.mapping_quality_mut() = noodles_sam::alignment::record::MappingQuality::new(mapq);
        r
    }

    fn bundler(flags: BundleFlags) -> GetBundles {
        let options = BundleOptions {
            random_seed: Some(7),
            ..BundleOptions::default()
        };
        GetBundles::new(
            &options,
            flags,
            BarcodeGetter::ReadId {
                cell_barcode: false,
                sep: "_".to_string(),
            },
        )
    }

    #[test]
    fn reads_at_one_position_form_one_bundle() {
        let mut b = bundler(BundleFlags::default());
        for name in ["r1_AAAA", "r2_AAAA", "r3_CCCC"] {
            let out = b.push(read(name, 101, 30), Some("chr1")).unwrap();
            assert!(out.is_empty(), "nothing should flush yet");
        }
        let items = b.finish();
        assert_eq!(items.len(), 1);
        let BundleItem::Bundle { bundle, .. } = &items[0] else {
            panic!("expected a bundle");
        };
        assert_eq!(bundle.len(), 2);
        assert_eq!(bundle[&b"AAAA".to_vec()].count, 2);
        assert_eq!(bundle[&b"CCCC".to_vec()].count, 1);
        assert_eq!(b.read_events()["Input Reads"], 3);
    }

    #[test]
    fn moving_past_the_flush_margin_emits_the_bundle() {
        let mut b = bundler(BundleFlags::default());
        assert!(b
            .push(read("r1_AAAA", 101, 30), Some("chr1"))
            .unwrap()
            .is_empty());
        // 2000 bases later: the first position can no longer gain reads.
        let out = b.push(read("r2_CCCC", 2101, 30), Some("chr1")).unwrap();
        assert_eq!(out.len(), 1);
        assert!(matches!(&out[0], BundleItem::Bundle { bundle, .. }
            if bundle.contains_key(&b"AAAA".to_vec())));
    }

    #[test]
    fn a_contig_change_flushes_everything() {
        let mut b = bundler(BundleFlags::default());
        b.push(read("r1_AAAA", 101, 30), Some("chr1")).unwrap();
        let out = b.push(read("r2_CCCC", 101, 30), Some("chr2")).unwrap();
        assert_eq!(out.len(), 1);
    }

    #[test]
    fn higher_mapping_quality_wins() {
        let mut b = bundler(BundleFlags::default());
        b.push(read("low_AAAA", 101, 10), Some("chr1")).unwrap();
        b.push(read("high_AAAA", 101, 40), Some("chr1")).unwrap();
        b.push(read("mid_AAAA", 101, 20), Some("chr1")).unwrap();

        let items = b.finish();
        let BundleItem::Bundle { bundle, .. } = &items[0] else {
            panic!("expected a bundle");
        };
        let kept = bundle[&b"AAAA".to_vec()].reads.one().unwrap();
        assert_eq!(kept.name().unwrap().to_vec(), b"high_AAAA".to_vec());
        assert_eq!(bundle[&b"AAAA".to_vec()].count, 3);
    }

    #[test]
    fn all_reads_mode_keeps_every_read() {
        let mut b = bundler(BundleFlags {
            all_reads: true,
            ..BundleFlags::default()
        });
        for name in ["r1_AAAA", "r2_AAAA"] {
            b.push(read(name, 101, 30), Some("chr1")).unwrap();
        }
        let items = b.finish();
        let BundleItem::Bundle { bundle, .. } = &items[0] else {
            panic!("expected a bundle");
        };
        assert_eq!(bundle[&b"AAAA".to_vec()].reads.all().len(), 2);
    }

    #[test]
    fn count_only_mode_keeps_no_reads() {
        let mut b = bundler(BundleFlags {
            only_count_reads: true,
            ..BundleFlags::default()
        });
        for name in ["r1_AAAA", "r2_AAAA"] {
            b.push(read(name, 101, 30), Some("chr1")).unwrap();
        }
        let items = b.finish();
        let BundleItem::Bundle { bundle, .. } = &items[0] else {
            panic!("expected a bundle");
        };
        assert_eq!(bundle[&b"AAAA".to_vec()].reads, BundleReads::CountOnly);
        assert_eq!(bundle[&b"AAAA".to_vec()].count, 2);
    }

    #[test]
    fn read2_is_dropped_unless_requested() {
        let mut r = read("r1_AAAA", 101, 30);
        *r.flags_mut() = Flags::SEGMENTED | Flags::LAST_SEGMENT;

        let mut b = bundler(BundleFlags::default());
        assert!(b.push(r.clone(), Some("chr1")).unwrap().is_empty());
        assert!(!b.read_events().contains_key("Input Reads"));

        let mut b = bundler(BundleFlags {
            return_read2: true,
            ..BundleFlags::default()
        });
        let out = b.push(r, Some("chr1")).unwrap();
        assert_eq!(out.len(), 1);
        assert!(matches!(out[0], BundleItem::Single(_)));
    }

    #[test]
    fn unmapped_reads_are_dropped_unless_requested() {
        let mut r = read("r1_AAAA", 101, 30);
        *r.flags_mut() = Flags::UNMAPPED;

        let mut b = bundler(BundleFlags::default());
        assert!(b.push(r.clone(), Some("chr1")).unwrap().is_empty());
        assert_eq!(b.read_events()["Single end unmapped"], 1);

        let mut b = bundler(BundleFlags {
            return_unmapped: true,
            ..BundleFlags::default()
        });
        let out = b.push(r, Some("chr1")).unwrap();
        assert_eq!(out.len(), 1);
        // Python counts 'Input Reads' twice on this path; see the comment in
        // `push`.
        assert_eq!(b.read_events()["Input Reads"], 2);
    }

    #[test]
    fn mapping_quality_threshold_filters_reads() {
        let options = BundleOptions {
            mapping_quality: Some(20),
            random_seed: Some(1),
            ..BundleOptions::default()
        };
        let mut b = GetBundles::new(
            &options,
            BundleFlags::default(),
            BarcodeGetter::ReadId {
                cell_barcode: false,
                sep: "_".to_string(),
            },
        );
        b.push(read("low_AAAA", 101, 10), Some("chr1")).unwrap();
        assert_eq!(b.read_events()["< MAPQ threshold"], 1);
        assert!(b.finish().is_empty());
    }

    /// A record with MAPQ 255 ("not available") must pass a threshold of 20,
    /// because `pysam` reports it as the raw 255. Treating `noodles`' `None`
    /// as 0 instead discarded every read of the package's own `chr19.bam`
    /// fixture, whose reads all carry 255 — see
    /// [`crate::sam::raw_mapping_quality`].
    #[test]
    fn unavailable_mapping_quality_is_treated_as_255() {
        let options = BundleOptions {
            mapping_quality: Some(20),
            random_seed: Some(1),
            ..BundleOptions::default()
        };
        let mut b = GetBundles::new(
            &options,
            BundleFlags::default(),
            BarcodeGetter::ReadId {
                cell_barcode: false,
                sep: "_".to_string(),
            },
        );
        let mut r = read("r1_AAAA", 101, 0);
        *r.mapping_quality_mut() = None;
        b.push(r, Some("chr1")).unwrap();
        assert!(!b.read_events().contains_key("< MAPQ threshold"));
        assert_eq!(b.finish().len(), 1);
    }

    #[test]
    fn strand_separates_bundles() {
        let mut b = bundler(BundleFlags::default());
        b.push(read("fwd_AAAA", 101, 30), Some("chr1")).unwrap();
        let mut rev = read("rev_AAAA", 101, 30);
        *rev.flags_mut() = Flags::REVERSE_COMPLEMENTED;
        b.push(rev, Some("chr1")).unwrap();

        // Forward groups at pos 100, reverse at pos 150 — different positions
        // as well as different keys.
        let items = b.finish();
        assert_eq!(items.len(), 2);
    }

    #[test]
    fn per_cell_separates_bundles_by_cell_barcode() {
        let options = BundleOptions {
            per_cell: true,
            random_seed: Some(1),
            ..BundleOptions::default()
        };
        let mut b = GetBundles::new(
            &options,
            BundleFlags::default(),
            BarcodeGetter::ReadId {
                cell_barcode: true,
                sep: "_".to_string(),
            },
        );
        b.push(read("r1_CELL1_AAAA", 101, 30), Some("chr1"))
            .unwrap();
        b.push(read("r2_CELL2_AAAA", 101, 30), Some("chr1"))
            .unwrap();
        let items = b.finish();
        assert_eq!(items.len(), 2);
        for item in &items {
            let BundleItem::Bundle { key, .. } = item else {
                panic!("expected bundles");
            };
            assert!(key.cell.is_some());
        }
    }

    #[test]
    fn tie_breaking_is_reproducible_for_a_seed() {
        let names = ["a_AAAA", "b_AAAA", "c_AAAA", "d_AAAA", "e_AAAA"];
        let winner = |seed: u64| -> Vec<u8> {
            let options = BundleOptions {
                random_seed: Some(seed),
                ..BundleOptions::default()
            };
            let mut b = GetBundles::new(
                &options,
                BundleFlags::default(),
                BarcodeGetter::ReadId {
                    cell_barcode: false,
                    sep: "_".to_string(),
                },
            );
            for n in names {
                b.push(read(n, 101, 30), Some("chr1")).unwrap();
            }
            let items = b.finish();
            let BundleItem::Bundle { bundle, .. } = &items[0] else {
                panic!("expected a bundle");
            };
            bundle[&b"AAAA".to_vec()]
                .reads
                .one()
                .unwrap()
                .name()
                .unwrap()
                .to_vec()
        };
        assert_eq!(winner(42), winner(42));
    }

    #[test]
    fn missing_umi_tags_are_counted_and_skipped() {
        let options = BundleOptions {
            random_seed: Some(1),
            ..BundleOptions::default()
        };
        let mut b = GetBundles::new(
            &options,
            BundleFlags::default(),
            BarcodeGetter::Tag(crate::sam::TagConfig::default()),
        );
        // No RX tag on this record.
        b.push(read("r1", 101, 30), Some("chr1")).unwrap();
        assert_eq!(
            b.read_events()["Read skipped, missing umi and/or cell tag"],
            1
        );
        assert!(b.finish().is_empty());
    }

    #[test]
    fn whole_contig_defers_flushing_until_the_contig_changes() {
        let options = BundleOptions {
            whole_contig: true,
            random_seed: Some(1),
            ..BundleOptions::default()
        };
        let mut b = GetBundles::new(
            &options,
            BundleFlags::default(),
            BarcodeGetter::ReadId {
                cell_barcode: false,
                sep: "_".to_string(),
            },
        );
        b.push(read("r1_AAAA", 101, 30), Some("chr1")).unwrap();
        // Far past the flush margin, but the same contig: nothing flushes.
        let out = b.push(read("r2_CCCC", 999_101, 30), Some("chr1")).unwrap();
        assert!(out.is_empty());
        let out = b.push(read("r3_GGGG", 101, 30), Some("chr2")).unwrap();
        assert_eq!(out.len(), 2);
    }
}
