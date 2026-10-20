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

use std::collections::{BTreeMap, VecDeque};

use indexmap::IndexMap;

use crate::alignment::BamRecord;
use crate::error::{PortError, Result};
use crate::options::{BundleOptions, ChimericPairs, DetectionMethod, UnmappedReads, UnpairedReads};
use crate::sam::{get_read_position, BarcodeGetter};
use crate::utilities;

/// The reads a bundle entry retained.
///
/// Python distinguishes these three modes only by which dictionary keys it
/// happens to populate, so reading `bundle[umi]["read"]` in `only_count_reads`
/// mode is a `KeyError`. The enum makes the distinction a compile-time one.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum BundleReads {
    /// The single best read for this key (the default for `dedup`).
    One(BamRecord),
    /// Every read for this key (`group`, via `all_reads`).
    Many(Vec<BamRecord>),
    /// The `--group-out` fields of every read (`group` without `--output-bam`).
    Summaries(Vec<ReadSummary>),
    /// No reads retained, only the count (`count`, via `only_count_reads`).
    #[default]
    CountOnly,
}

/// What `group --group-out` reports of a read, kept in place of the read
/// when no reads are written out.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadSummary {
    /// The read name, empty for `*`.
    pub name: Box<[u8]>,
    /// The reference-sequence index.
    pub reference_id: Option<usize>,
    /// The grouping coordinate.
    pub position: i64,
}

impl ReadSummary {
    fn of(record: &BamRecord, position: i64) -> Self {
        Self {
            name: record.name().unwrap_or_default().into(),
            reference_id: record.reference_id(),
            position,
        }
    }
}

impl BundleReads {
    /// The single retained read, if this entry holds one.
    #[must_use]
    pub fn one(&self) -> Option<&BamRecord> {
        match self {
            BundleReads::One(r) => Some(r),
            _ => None,
        }
    }

    /// Every retained read, empty when none were kept.
    #[must_use]
    pub fn all(&self) -> &[BamRecord] {
        match self {
            BundleReads::Many(v) => v,
            _ => &[],
        }
    }

    /// Every retained summary, empty when none were kept.
    #[must_use]
    pub fn summaries(&self) -> &[ReadSummary] {
        match self {
            BundleReads::Summaries(v) => v,
            _ => &[],
        }
    }
}

/// One UMI's entry in a bundle
/// (Python source: the `{"read": ..., "count": ...}` inner dictionary).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BundleEntry {
    /// The retained read(s).
    pub reads: BundleReads,
    /// How many reads shared this key.
    pub count: u64,
}

/// UMI to entry, for one position/gene and cell barcode, in the order the
/// UMIs were first seen (Python source: the dictionary `get_bundles` yields).
///
/// The order is Python's `dict` insertion order, and it matters: the
/// clustering stage takes the UMIs in that order, which decides the order the
/// deduplicated reads are written in.
pub type Bundle = IndexMap<Vec<u8>, BundleEntry>;

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
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CoordinateKey {
    /// Whether the read is on the reverse strand.
    pub is_reverse: bool,
    /// The first splice offset under `--spliced-is-unique`, else 0 (Python's `False`).
    pub spliced: usize,
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
    Single(BamRecord),
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
/// `dedup`, `group` and `count` pass this to
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
/// bundle.insert(b"TTTT".to_vec(), BundleEntry { reads: BundleReads::CountOnly, count: 1 });
/// bundle.insert(b"AAAA".to_vec(), BundleEntry { reads: BundleReads::CountOnly, count: 1 });
///
/// assert_eq!(insertion_order(&bundle), vec![b"TTTT".to_vec(), b"AAAA".to_vec()]);
/// ```
#[must_use]
pub fn insertion_order(bundle: &Bundle) -> Vec<Vec<u8>> {
    bundle.keys().cloned().collect()
}

/// Which reads the bundler should retain and pass through
/// (Python source: the five boolean constructor flags of `get_bundles`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct BundleFlags {
    /// Retain counts only, no reads (`count`).
    pub only_count_reads: bool,
    /// Retain every read per key (`group`).
    pub all_reads: bool,
    /// With `all_reads`, retain a [`ReadSummary`] instead of each read.
    pub summarise_reads: bool,
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

/// Python's bundle key without the gene, plus the UMI.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct SlotKey {
    coordinate: Option<CoordinateKey>,
    barcode: Box<[u8]>,
    cell_len: Option<u32>,
}

impl SlotKey {
    fn new(coordinate: Option<CoordinateKey>, cell: Option<Vec<u8>>, umi: &[u8]) -> Self {
        let cell_len = cell
            .as_ref()
            .map(|c| u32::try_from(c.len()).unwrap_or(u32::MAX));
        let cell = cell.unwrap_or_default();
        let mut barcode = Vec::with_capacity(cell.len() + umi.len());
        barcode.extend_from_slice(&cell);
        barcode.extend_from_slice(umi);
        Self {
            coordinate,
            barcode: barcode.into_boxed_slice(),
            cell_len,
        }
    }

    fn cell(&self) -> Option<&[u8]> {
        self.cell_len.map(|n| &self.barcode[..n as usize])
    }

    fn umi(&self) -> &[u8] {
        &self.barcode[self.cell_len.map_or(0, |n| n as usize)..]
    }
}

/// Python's `reads_dict[pos][key][umi]` and `read_counts[pos][key][umi]`.
#[derive(Debug)]
struct Slot {
    reads: BundleReads,
    count: u64,
    ties: u64,
}

type PositionSlots = IndexMap<SlotKey, Slot>;

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
/// # let records: Vec<umi_tools_core::BamRecord> = Vec::new();
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

    reads_dict: BTreeMap<PositionKey, PositionSlots>,
    ready: VecDeque<Queued>,
    emitting: Option<Emitting>,
    /// Gene to the transcripts observed for it, for the metacontig path.
    observed_contigs: BTreeMap<String, indexmap::IndexSet<String>>,
    metacontig_contig: Option<crate::sam::MetaContigMap>,
    contig_metacontig: BTreeMap<String, String>,

    read_events: BTreeMap<String, u64>,
    last_pos: PositionKey,
    last_chr: Option<usize>,
    start: i64,
    current_chr: Option<usize>,

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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Replacement {
    /// Keep the stored read.
    Keep,
    /// Replace it, and reset the tie counter (higher mapping quality).
    Replace,
    /// Replace it and still draw (Python's multimapping branch has no `return`).
    ReplaceAndDraw,
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
            ready: VecDeque::new(),
            emitting: None,
            observed_contigs: BTreeMap::new(),
            metacontig_contig: None,
            contig_metacontig: BTreeMap::new(),
            read_events: BTreeMap::new(),
            last_pos: PositionKey::default(),
            last_chr: None,
            start: 0,
            current_chr: None,
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

    /// How many `(key, umi)` entries are waiting to be flushed.
    #[must_use]
    pub fn pending(&self) -> usize {
        self.reads_dict.values().map(IndexMap::len).sum()
    }

    fn event(&mut self, name: &str) {
        if let Some(n) = self.read_events.get_mut(name) {
            *n += 1;
        } else {
            self.read_events.insert(name.to_string(), 1);
        }
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

    fn multimap_tag(record: &BamRecord, tag: [u8; 2]) -> Option<i64> {
        use noodles_sam::alignment::record::data::field::Value;
        let view = record.view();
        match view.data().get(&tag)?.ok()? {
            Value::Int8(v) => Some(i64::from(v)),
            Value::UInt8(v) => Some(i64::from(v)),
            Value::Int16(v) => Some(i64::from(v)),
            Value::UInt16(v) => Some(i64::from(v)),
            Value::Int32(v) => Some(i64::from(v)),
            Value::UInt32(v) => Some(i64::from(v)),
            _ => None,
        }
    }

    fn xt_is_unique(record: &BamRecord) -> bool {
        use noodles_sam::alignment::record::data::field::Value;
        let view = record.view();
        match view.data().get(b"XT") {
            Some(Ok(Value::Character(b'U'))) => true,
            Some(Ok(Value::String(s))) => <[u8]>::eq(s.as_ref(), b"U"),
            _ => false,
        }
    }

    /// Decides between the stored read and a new one, on mapping quality and
    /// then the multimapping tag (Python source: the comparison chain in the
    /// `else` branch of `get_bundles.update_dicts`).
    fn compare(
        detection_method: Option<DetectionMethod>,
        current: &BamRecord,
        candidate: &BamRecord,
    ) -> Replacement {
        let current_mapq = current.mapping_quality();
        let new_mapq = candidate.mapping_quality();
        if current_mapq > new_mapq {
            return Replacement::Keep;
        }
        if current_mapq < new_mapq {
            return Replacement::Replace;
        }

        match detection_method {
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
                    (Some(a), Some(b)) if a > b => Replacement::ReplaceAndDraw,
                    _ => Replacement::Tied,
                }
            }
            Some(DetectionMethod::Xt) => {
                if Self::xt_is_unique(current) {
                    Replacement::Keep
                } else if Self::xt_is_unique(candidate) {
                    Replacement::ReplaceAndDraw
                } else {
                    Replacement::Tied
                }
            }
            None => Replacement::Tied,
        }
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
        record: BamRecord,
        pos: PositionKey,
        key: BundleKey,
        umi: Vec<u8>,
    ) -> Result<()> {
        self.file(record, pos, SlotKey::new(key.coordinate, key.cell, &umi))
    }

    fn file(&mut self, record: BamRecord, pos: PositionKey, key: SlotKey) -> Result<()> {
        use indexmap::map::Entry;

        let Self {
            reads_dict,
            flags,
            options,
            rng,
            ..
        } = self;
        // Under `--per-gene` the position is the gene, so use the read's own.
        let summary = if flags.all_reads && flags.summarise_reads {
            let position = match &pos {
                PositionKey::Coordinate(p) => *p,
                PositionKey::Gene(_) => {
                    get_read_position(&record, options.soft_clip_threshold)?.pos
                }
            };
            Some(ReadSummary::of(&record, position))
        } else {
            None
        };

        let slots = reads_dict
            .entry(pos)
            .or_insert_with(|| PositionSlots::with_capacity(1));

        let slot = match slots.entry(key) {
            Entry::Vacant(vacant) => {
                let reads = if let Some(summary) = summary {
                    BundleReads::Summaries(vec![summary])
                } else if flags.all_reads {
                    BundleReads::Many(vec![record])
                } else if flags.only_count_reads {
                    BundleReads::CountOnly
                } else {
                    BundleReads::One(record)
                };
                vacant.insert(Slot {
                    reads,
                    count: 1,
                    ties: 0,
                });
                return Ok(());
            }
            Entry::Occupied(occupied) => occupied.into_mut(),
        };

        slot.count += 1;
        if flags.all_reads {
            match (&mut slot.reads, summary) {
                (BundleReads::Summaries(v), Some(summary)) => v.push(summary),
                (BundleReads::Many(v), None) => v.push(record),
                _ => {}
            }
            return Ok(());
        }
        if flags.only_count_reads {
            return Ok(());
        }

        let decision = match &slot.reads {
            BundleReads::One(current) => Self::compare(options.detection_method, current, &record),
            _ => Replacement::Keep,
        };
        match decision {
            Replacement::Keep => {}
            Replacement::Replace => {
                slot.reads = BundleReads::One(record);
                slot.ties = 0;
            }
            Replacement::ReplaceAndDraw => {
                slot.reads = BundleReads::One(record);
                slot.ties = 1;
                let _ = rng.random();
            }
            Replacement::Tied => {
                slot.ties += 1;
                let probability = 1.0 / (slot.ties as f64);
                if rng.random() < probability {
                    slot.reads = BundleReads::One(record);
                }
            }
        }
        Ok(())
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
    ///   past the previous read, or on a contig change.
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

    fn flush(&mut self, keys: Vec<PositionKey>) {
        for position in keys {
            if let Some(slots) = self.reads_dict.remove(&position) {
                self.ready.push_back(Queued::Position(position, slots));
            }
        }
    }

    fn next_ready(&mut self) -> Option<BundleItem> {
        loop {
            if let Some(bundle) = self.emitting.as_mut().and_then(Emitting::next_bundle) {
                return Some(bundle);
            }
            self.emitting = None;
            match self.ready.pop_front()? {
                Queued::Single(record) => return Some(BundleItem::Single(record)),
                Queued::Position(position, slots) => {
                    self.emitting = Some(Emitting::new(position, slots));
                }
            }
        }
    }

    /// Feeds one record to the bundler
    /// (Python source: the body of `get_bundles.__call__`'s loop).
    ///
    /// `reference_name` is the contig name, the gene under `--per-contig`.
    ///
    /// Returns the bundles (and pass-through reads) this record completed.
    /// They are produced lazily; items not taken are handed out by the next call.
    ///
    /// # Errors
    ///
    /// - [`PortError::Value`] from the coordinate arithmetic or a malformed
    ///   read name.
    /// - [`PortError::Alignment`] for a mapped record with no position, or
    ///   unparseable tags.
    pub fn push(&mut self, record: BamRecord, reference_name: Option<&str>) -> Result<Ready<'_>> {
        self.file_read(record, reference_name)?;
        Ok(Ready { bundler: self })
    }

    fn file_read(&mut self, record: BamRecord, reference_name: Option<&str>) -> Result<()> {
        let flags = record.flags();

        if flags.is_segmented() && flags.is_last_segment() {
            if self.flags.return_read2 && (!flags.is_unmapped() || self.flags.return_unmapped) {
                {
                    self.ready.push_back(Queued::Single(record));
                    return Ok(());
                }
            }
            return Ok(());
        }
        self.event("Input Reads");

        // From here on, only read 1s (or unpaired reads).
        if self.options.paired {
            if flags.is_segmented() {
                self.event("Read pairs");
            } else {
                self.event("Unpaired reads");
                match self.options.unpaired_reads {
                    UnpairedReads::Discard => return Ok(()),
                    // Python has no `continue` here: the read is bundled too.
                    UnpairedReads::Output => {
                        self.ready.push_back(Queued::Single(record.clone()));
                    }
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
                {
                    self.ready.push_back(Queued::Single(record));
                    return Ok(());
                }
            }
            return Ok(());
        }

        if self.options.paired && flags.is_mate_unmapped() {
            self.event("Read 2 unmapped");
            if self.options.unmapped_reads != UnmappedReads::Use {
                if self.flags.return_unmapped {
                    {
                        self.ready.push_back(Queued::Single(record));
                        return Ok(());
                    }
                }
                return Ok(());
            }
        }

        if flags.is_segmented() && record.reference_id() != record.mate_reference_id() {
            self.event("Chimeric read pair");
            match self.options.chimeric_pairs {
                ChimericPairs::Discard => return Ok(()),
                ChimericPairs::Output => {
                    self.ready.push_back(Queued::Single(record));
                    return Ok(());
                }
                ChimericPairs::Use => {}
            }
        }

        // `if self.options.subset:`, so 0 turns subsetting off.
        if let Some(subset) = self.options.subset.filter(|&s| s != 0.0) {
            if self.next_uniform() >= subset {
                self.event("Randomly excluded");
                return Ok(());
            }
        }

        if let Some(threshold) = self.options.mapping_quality {
            if u32::from(record.mapping_quality()) < threshold {
                self.event("< MAPQ threshold");
                return Ok(());
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
                let name = String::from_utf8_lossy(record.name().unwrap_or_default()).into_owned();
                utilities::warn(&format!(
                    "At least one read is missing UMI and/or cell tag(s): {name}"
                ));
            }
            self.event(msg);
            return Ok(());
        };

        self.current_chr = record.reference_id();

        let (pos, coordinate);

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
                let (Some(assigned), Some(gene)) = (
                    crate::sam::tag_text(&record, assigned_tag)?,
                    crate::sam::tag_text(&record, gene_tag)?,
                ) else {
                    self.event("Read skipped, no tag");
                    return Ok(());
                };
                if gene.is_empty() {
                    let msg = "Read skipped - gene string is empty";
                    if !self.read_events.contains_key(msg) {
                        utilities::warn("Assigned gene is empty string. First such read");
                    }
                    self.event(msg);
                    return Ok(());
                }
                if let Some(skip) = &self.options.skip_regex {
                    if skip.is_match(&assigned) {
                        self.event("Read skipped - assigned tag matches skip_regex");
                        return Ok(());
                    }
                }
                gene
            };

            if self.last_chr.is_some() {
                let (do_output, keys) = self.check_output();
                if do_output {
                    self.flush(keys);
                }
            }
            self.last_chr = self.current_chr;
            if let Some(transcript) = transcript {
                self.observed_contigs
                    .entry(gene.clone())
                    .or_default()
                    .insert(transcript);
            }
            pos = PositionKey::Gene(gene);
            self.last_pos = pos.clone();
            coordinate = None;
        } else {
            let position = get_read_position(&record, self.options.soft_clip_threshold)?;
            self.start = position.start;

            let (do_output, keys) = self.check_output();
            if do_output {
                self.flush(keys);
            }
            self.last_pos = PositionKey::Coordinate(position.start);
            self.last_chr = self.current_chr;

            let read_length = if self.options.read_length {
                record.sequence_length()
            } else {
                0
            };
            let template_length = if self.options.ignore_tlen || !self.options.paired {
                0
            } else {
                record.template_length()
            };

            pos = PositionKey::Coordinate(position.pos);
            coordinate = Some(CoordinateKey {
                is_reverse: flags.is_reverse_complemented(),
                spliced: if self.options.spliced {
                    position.is_spliced.unwrap_or(0)
                } else {
                    0
                },
                template_length,
                read_length,
            });
        }

        self.file(
            record,
            pos,
            SlotKey::new(coordinate, barcode.cell, &barcode.umi),
        )
    }

    /// Flushes every remaining bundle
    /// (Python source: the trailing `for p in sorted(self.reads_dict.keys())`
    /// loop after `get_bundles.__call__`'s main loop).
    pub fn finish(&mut self) -> Ready<'_> {
        while let Some((position, slots)) = self.reads_dict.pop_first() {
            self.ready.push_back(Queued::Position(position, slots));
        }
        Ready { bundler: self }
    }
}

/// The items a [`GetBundles::push`] or [`GetBundles::finish`] made ready.
pub struct Ready<'a> {
    bundler: &'a mut GetBundles,
}

impl Iterator for Ready<'_> {
    type Item = BundleItem;

    fn next(&mut self) -> Option<BundleItem> {
        self.bundler.next_ready()
    }
}

#[derive(Debug)]
enum Queued {
    Single(BamRecord),
    Position(PositionKey, PositionSlots),
}

#[derive(Debug)]
struct Emitting {
    gene: Option<String>,
    entries: std::iter::Peekable<std::vec::IntoIter<(SlotKey, Slot)>>,
}

impl Emitting {
    fn new(position: PositionKey, slots: PositionSlots) -> Self {
        let gene = match position {
            PositionKey::Gene(gene) => Some(gene),
            PositionKey::Coordinate(_) => None,
        };
        let mut entries: Vec<(SlotKey, Slot)> = slots.into_iter().collect();
        // Stable, so each key's UMIs stay in first-seen order.
        entries.sort_by(|(a, _), (b, _)| (a.coordinate, a.cell()).cmp(&(b.coordinate, b.cell())));
        Self {
            gene,
            entries: entries.into_iter().peekable(),
        }
    }

    fn next_bundle(&mut self) -> Option<BundleItem> {
        let (first, slot) = self.entries.next()?;
        let mut bundle = Bundle::with_capacity(1);
        bundle.insert(
            first.umi().to_vec(),
            BundleEntry {
                reads: slot.reads,
                count: slot.count,
            },
        );
        while let Some((next, slot)) = self
            .entries
            .next_if(|(next, _)| next.coordinate == first.coordinate && next.cell() == first.cell())
        {
            bundle.insert(
                next.umi().to_vec(),
                BundleEntry {
                    reads: slot.reads,
                    count: slot.count,
                },
            );
        }
        Some(BundleItem::Bundle {
            bundle,
            key: BundleKey {
                gene: self.gene.clone(),
                coordinate: first.coordinate,
                cell: first.cell().map(<[u8]>::to_vec),
            },
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_core::Position;
    use noodles_sam::alignment::record::cigar::op::{Kind, Op};
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record_buf::Cigar as CigarBuf;
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    use std::num::NonZeroUsize;

    fn on(r: RecordBuf, contig: &str) -> BamRecord {
        let mut header = noodles_sam::Header::builder();
        for name in ["chr1", "chr2"] {
            header = header.add_reference_sequence(
                name,
                Map::<ReferenceSequence>::new(NonZeroUsize::new(10_000_000).unwrap()),
            );
        }
        let mut r = r;
        *r.reference_sequence_id_mut() = Some(usize::from(contig == "chr2"));
        BamRecord::from_alignment_record(&header.build(), &r).unwrap()
    }

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
            let out = b
                .push(on(read(name, 101, 30), "chr1"), Some("chr1"))
                .unwrap()
                .collect::<Vec<_>>();
            assert!(out.is_empty(), "nothing should flush yet");
        }
        let items = b.finish().collect::<Vec<_>>();
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
            .push(on(read("r1_AAAA", 101, 30), "chr1"), Some("chr1"))
            .unwrap()
            .collect::<Vec<_>>()
            .is_empty());
        // 2000 bases later: the first position can no longer gain reads.
        let out = b
            .push(on(read("r2_CCCC", 2101, 30), "chr1"), Some("chr1"))
            .unwrap()
            .collect::<Vec<_>>();
        assert_eq!(out.len(), 1);
        assert!(matches!(&out[0], BundleItem::Bundle { bundle, .. }
            if bundle.contains_key(&b"AAAA".to_vec())));
    }

    #[test]
    fn a_contig_change_flushes_everything() {
        let mut b = bundler(BundleFlags::default());
        b.push(on(read("r1_AAAA", 101, 30), "chr1"), Some("chr1"))
            .unwrap()
            .for_each(drop);
        let out = b
            .push(on(read("r2_CCCC", 101, 30), "chr2"), Some("chr2"))
            .unwrap()
            .collect::<Vec<_>>();
        assert_eq!(out.len(), 1);
    }

    #[test]
    fn higher_mapping_quality_wins() {
        let mut b = bundler(BundleFlags::default());
        b.push(on(read("low_AAAA", 101, 10), "chr1"), Some("chr1"))
            .unwrap()
            .for_each(drop);
        b.push(on(read("high_AAAA", 101, 40), "chr1"), Some("chr1"))
            .unwrap()
            .for_each(drop);
        b.push(on(read("mid_AAAA", 101, 20), "chr1"), Some("chr1"))
            .unwrap()
            .for_each(drop);

        let items = b.finish().collect::<Vec<_>>();
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
            b.push(on(read(name, 101, 30), "chr1"), Some("chr1"))
                .unwrap()
                .for_each(drop);
        }
        let items = b.finish().collect::<Vec<_>>();
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
            b.push(on(read(name, 101, 30), "chr1"), Some("chr1"))
                .unwrap()
                .for_each(drop);
        }
        let items = b.finish().collect::<Vec<_>>();
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
        assert!(b
            .push(on(r.clone(), "chr1"), Some("chr1"))
            .unwrap()
            .collect::<Vec<_>>()
            .is_empty());
        assert!(!b.read_events().contains_key("Input Reads"));

        let mut b = bundler(BundleFlags {
            return_read2: true,
            ..BundleFlags::default()
        });
        let out = b
            .push(on(r, "chr1"), Some("chr1"))
            .unwrap()
            .collect::<Vec<_>>();
        assert_eq!(out.len(), 1);
        assert!(matches!(out[0], BundleItem::Single(_)));
    }

    #[test]
    fn unmapped_reads_are_dropped_unless_requested() {
        let mut r = read("r1_AAAA", 101, 30);
        *r.flags_mut() = Flags::UNMAPPED;

        let mut b = bundler(BundleFlags::default());
        assert!(b
            .push(on(r.clone(), "chr1"), Some("chr1"))
            .unwrap()
            .collect::<Vec<_>>()
            .is_empty());
        assert_eq!(b.read_events()["Single end unmapped"], 1);

        let mut b = bundler(BundleFlags {
            return_unmapped: true,
            ..BundleFlags::default()
        });
        let out = b
            .push(on(r, "chr1"), Some("chr1"))
            .unwrap()
            .collect::<Vec<_>>();
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
        b.push(on(read("low_AAAA", 101, 10), "chr1"), Some("chr1"))
            .unwrap()
            .for_each(drop);
        assert_eq!(b.read_events()["< MAPQ threshold"], 1);
        assert!(b.finish().collect::<Vec<_>>().is_empty());
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
        b.push(on(r, "chr1"), Some("chr1")).unwrap().for_each(drop);
        assert!(!b.read_events().contains_key("< MAPQ threshold"));
        assert_eq!(b.finish().collect::<Vec<_>>().len(), 1);
    }

    #[test]
    fn strand_separates_bundles() {
        let mut b = bundler(BundleFlags::default());
        b.push(on(read("fwd_AAAA", 101, 30), "chr1"), Some("chr1"))
            .unwrap()
            .for_each(drop);
        let mut rev = read("rev_AAAA", 101, 30);
        *rev.flags_mut() = Flags::REVERSE_COMPLEMENTED;
        b.push(on(rev, "chr1"), Some("chr1"))
            .unwrap()
            .for_each(drop);

        // Forward groups at pos 100, reverse at pos 150 — different positions
        // as well as different keys.
        let items = b.finish().collect::<Vec<_>>();
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
        b.push(on(read("r1_CELL1_AAAA", 101, 30), "chr1"), Some("chr1"))
            .unwrap()
            .for_each(drop);
        b.push(on(read("r2_CELL2_AAAA", 101, 30), "chr1"), Some("chr1"))
            .unwrap()
            .for_each(drop);
        let items = b.finish().collect::<Vec<_>>();
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
                b.push(on(read(n, 101, 30), "chr1"), Some("chr1"))
                    .unwrap()
                    .for_each(drop);
            }
            let items = b.finish().collect::<Vec<_>>();
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
        b.push(on(read("r1", 101, 30), "chr1"), Some("chr1"))
            .unwrap()
            .for_each(drop);
        assert_eq!(
            b.read_events()["Read skipped, missing umi and/or cell tag"],
            1
        );
        assert!(b.finish().collect::<Vec<_>>().is_empty());
    }

    fn with_nh(mut r: RecordBuf, nh: u8) -> RecordBuf {
        use noodles_sam::alignment::record::data::field::Tag;
        use noodles_sam::alignment::record_buf::data::field::Value;
        r.data_mut().insert(Tag::from(*b"NH"), Value::UInt8(nh));
        r
    }

    fn read_id_bundler(options: &BundleOptions) -> GetBundles {
        GetBundles::new(
            options,
            BundleFlags::default(),
            BarcodeGetter::ReadId {
                cell_barcode: false,
                sep: "_".to_string(),
            },
        )
    }

    #[test]
    fn a_better_multimapping_tag_replaces_and_still_draws() {
        let options = BundleOptions {
            detection_method: Some(DetectionMethod::Nh),
            random_seed: Some(3),
            ..BundleOptions::default()
        };
        let mut rng = crate::mersenne::MersenneTwister::seed_from_u64(3);
        let _consumed_by_b = rng.random();
        let c_wins = rng.random() < 0.5;

        let mut b = read_id_bundler(&options);
        for (name, nh) in [("a_AAAA", 2), ("b_AAAA", 1), ("c_AAAA", 1)] {
            b.push(on(with_nh(read(name, 101, 30), nh), "chr1"), Some("chr1"))
                .unwrap()
                .for_each(drop);
        }
        let items = b.finish().collect::<Vec<_>>();
        let BundleItem::Bundle { bundle, .. } = &items[0] else {
            panic!("expected a bundle");
        };
        let kept = bundle[&b"AAAA".to_vec()]
            .reads
            .one()
            .unwrap()
            .name()
            .unwrap();
        assert_eq!(
            kept,
            if c_wins {
                &b"c_AAAA"[..]
            } else {
                &b"b_AAAA"[..]
            }
        );
    }

    #[test]
    fn spliced_reads_are_keyed_by_where_the_splice_starts() {
        let options = BundleOptions {
            spliced: true,
            random_seed: Some(1),
            ..BundleOptions::default()
        };
        let mut b = read_id_bundler(&options);
        for (name, ops) in [
            (
                "a_AAAA",
                vec![
                    Op::new(Kind::Match, 10),
                    Op::new(Kind::Skip, 100),
                    Op::new(Kind::Match, 40),
                ],
            ),
            (
                "b_AAAA",
                vec![
                    Op::new(Kind::Match, 20),
                    Op::new(Kind::Skip, 100),
                    Op::new(Kind::Match, 30),
                ],
            ),
            (
                "c_AAAA",
                vec![
                    Op::new(Kind::Match, 20),
                    Op::new(Kind::Skip, 50),
                    Op::new(Kind::Match, 30),
                ],
            ),
        ] {
            let mut r = read(name, 101, 30);
            *r.cigar_mut() = ops.into_iter().collect::<CigarBuf>();
            b.push(on(r, "chr1"), Some("chr1")).unwrap().for_each(drop);
        }
        let splices: Vec<usize> = b
            .finish()
            .map(|item| match item {
                BundleItem::Bundle { key, .. } => key.coordinate.unwrap().spliced,
                BundleItem::Single(_) => panic!("expected bundles"),
            })
            .collect();
        assert_eq!(splices, vec![10, 20]);
    }

    #[test]
    fn a_subset_of_zero_keeps_every_read() {
        let options = BundleOptions {
            subset: Some(0.0),
            random_seed: Some(1),
            ..BundleOptions::default()
        };
        let mut b = read_id_bundler(&options);
        b.push(on(read("a_AAAA", 101, 30), "chr1"), Some("chr1"))
            .unwrap()
            .for_each(drop);
        assert!(!b.read_events().contains_key("Randomly excluded"));
        assert_eq!(b.finish().count(), 1);
    }

    #[test]
    fn output_unpaired_reads_are_also_bundled() {
        let options = BundleOptions {
            paired: true,
            unpaired_reads: UnpairedReads::Output,
            random_seed: Some(1),
            ..BundleOptions::default()
        };
        let mut b = read_id_bundler(&options);
        let out = b
            .push(on(read("a_AAAA", 101, 30), "chr1"), Some("chr1"))
            .unwrap()
            .collect::<Vec<_>>();
        assert!(matches!(out[..], [BundleItem::Single(_)]));
        assert_eq!(b.finish().count(), 1);
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
        b.push(on(read("r1_AAAA", 101, 30), "chr1"), Some("chr1"))
            .unwrap()
            .for_each(drop);
        // Far past the flush margin, but the same contig: nothing flushes.
        let out = b
            .push(on(read("r2_CCCC", 999_101, 30), "chr1"), Some("chr1"))
            .unwrap()
            .collect::<Vec<_>>();
        assert!(out.is_empty());
        let out = b
            .push(on(read("r3_GGGG", 101, 30), "chr2"), Some("chr2"))
            .unwrap()
            .collect::<Vec<_>>();
        assert_eq!(out.len(), 2);
    }
}
