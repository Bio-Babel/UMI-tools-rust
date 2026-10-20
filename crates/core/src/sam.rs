//! Reading UMIs and cell barcodes out of alignment records, and the
//! coordinate arithmetic that groups them
//! (Python source: `umi_tools/sam_methods.py`).
//!
//! `umi_tools` accepts the barcode in three places, selected by
//! `--extract-umi-method`: appended to the read name by `umi_tools extract`,
//! in a SAM tag written by an aligner or a 10x pipeline, or in the
//! `umis`-style `:UMI_<seq>:` read-name convention.

use std::collections::{BTreeMap, BTreeSet, HashSet};

use indexmap::{IndexMap, IndexSet};
use std::io::BufRead;
use std::path::Path;

use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::Value;

use crate::alignment::{AlignmentSink, AlignmentSource, BamRecord, OpenedAlignment, ReferenceScan};
use crate::error::{PortError, Result};
use crate::record::TagValue;
use crate::utilities;

/// The raw mapping quality byte, as `pysam`'s `read.mapping_quality` reports it.
///
/// `noodles` models the SAM specification's "mapping quality is not available"
/// sentinel (255) as `None`, which is the more honest type — but `pysam`
/// returns the raw byte, and `umi_tools` compares it numerically against
/// `--mapping-quality` and against other reads' values. Mapping `None` back to
/// 255 is what keeps those comparisons equivalent.
///
/// This is not a hypothetical: every read in the package's own `chr19.bam`
/// fixture carries MAPQ 255, so treating `None` as 0 makes
/// `--mapping-quality=20` discard the entire file. Caught by the differential
/// test against the Python oracle, not by reading the source.
///
/// # Examples
///
/// ```
/// use noodles_sam::{alignment::RecordBuf, Header};
/// use umi_tools_core::{sam::raw_mapping_quality, BamRecord};
/// // A record with no mapping quality set reads back as 255, not 0.
/// let record = BamRecord::from_alignment_record(&Header::default(), &RecordBuf::default())?;
/// assert_eq!(raw_mapping_quality(&record), 255);
/// # Ok::<(), std::io::Error>(())
/// ```
#[must_use]
pub fn raw_mapping_quality(record: &BamRecord) -> u8 {
    record.mapping_quality()
}

/// A UMI and, optionally, the cell barcode it was found with
/// (Python source: the `(umi, cell)` tuple every barcode getter returns).
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default)]
pub struct Barcode {
    /// The UMI.
    pub umi: Vec<u8>,
    /// The cell barcode, absent unless `--per-cell` was requested.
    pub cell: Option<Vec<u8>>,
}

/// Extracts the barcode from a read name written by `umi_tools extract`
/// (Python source: `sam_methods.get_barcode_read_id`).
///
/// The UMI is the last `sep`-delimited field of the read name and, with
/// `cell_barcode`, the cell barcode is the second-to-last — the layout
/// [`crate::extract::add_barcodes_to_identifier`] writes.
///
/// # Errors
///
/// Returns [`PortError::Value`] when the name has too few fields, with
/// Python's message. Python wraps the whole body in a bare `except:` and
/// re-raises `ValueError`, so *any* failure surfaces the same way.
///
/// # Examples
///
/// ```
/// use umi_tools_core::sam::get_barcode_read_id;
/// let bc = get_barcode_read_id("read1_ACGT", false, "_")?;
/// assert_eq!(bc.umi, b"ACGT");
/// assert!(bc.cell.is_none());
///
/// let bc = get_barcode_read_id("read1_CELL_ACGT", true, "_")?;
/// assert_eq!(bc.cell.as_deref(), Some(&b"CELL"[..]));
/// # Ok::<(), umi_tools_core::PortError>(())
/// ```
pub fn get_barcode_read_id(read_name: &str, cell_barcode: bool, sep: &str) -> Result<Barcode> {
    let fields: Vec<&str> = read_name.split(sep).collect();
    let bad = || {
        PortError::value(
            "Could not extract UMI +/- cell barcode from the read\
             ID, please check UMI is encoded in the read name",
        )
    };
    let umi = fields.last().ok_or_else(bad)?.as_bytes().to_vec();
    let cell = if cell_barcode {
        let n = fields.len();
        if n < 2 {
            return Err(bad());
        }
        Some(fields[n - 2].as_bytes().to_vec())
    } else {
        None
    };
    Ok(Barcode { umi, cell })
}

/// How the barcode tags should be read
/// (Python source: the eight keyword parameters of `get_barcode_tag`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TagConfig {
    /// Also read a cell barcode.
    pub cell_barcode: bool,
    /// The two-character tag holding the UMI (`--umi-tag`, default `RX`).
    pub umi_tag: [u8; 2],
    /// The tag holding the cell barcode (`--cell-tag`).
    pub cell_tag: Option<[u8; 2]>,
    /// Keep only the part of the UMI before this separator.
    pub umi_tag_split: Option<String>,
    /// Remove every occurrence of this delimiter from the UMI.
    pub umi_tag_delim: Option<String>,
    /// Keep only the part of the cell barcode before this separator.
    ///
    /// Defaults to `-`, because 10x pipelines append a GEM well suffix
    /// (`GATAGATACCTAGATA-1`) that is not part of the barcode.
    pub cell_tag_split: Option<String>,
    /// Remove every occurrence of this delimiter from the cell barcode.
    pub cell_tag_delim: Option<String>,
}

impl Default for TagConfig {
    /// Python's defaults: `umi_tag='RX'`, `cell_tag_split="-"`, the rest unset.
    fn default() -> Self {
        Self {
            cell_barcode: false,
            umi_tag: *b"RX",
            cell_tag: None,
            umi_tag_split: None,
            umi_tag_delim: None,
            cell_tag_split: Some("-".to_string()),
            cell_tag_delim: None,
        }
    }
}

/// Reads a SAM tag as text, whatever scalar type it was stored as
/// (Python source: `pysam`'s `read.get_tag`, which likewise coerces).
///
/// Returns `None` when the tag is absent or holds an array, which has no
/// single textual value.
pub fn tag_text(record: &BamRecord, tag: [u8; 2]) -> Result<Option<String>> {
    let view = record.view();
    let data = view.data();
    let Some(value) = data.get(&tag) else {
        return Ok(None);
    };
    let value = value.map_err(|e| PortError::alignment(e.to_string()))?;
    Ok(match value {
        Value::String(s) | Value::Hex(s) => Some(String::from_utf8_lossy(s.as_ref()).into_owned()),
        Value::Character(c) => Some((c as char).to_string()),
        Value::Int8(v) => Some(v.to_string()),
        Value::UInt8(v) => Some(v.to_string()),
        Value::Int16(v) => Some(v.to_string()),
        Value::UInt16(v) => Some(v.to_string()),
        Value::Int32(v) => Some(v.to_string()),
        Value::UInt32(v) => Some(v.to_string()),
        Value::Float(v) => Some(v.to_string()),
        Value::Array(_) => None,
    })
}

fn apply_split_and_delim(value: String, split: Option<&String>, delim: Option<&String>) -> String {
    let mut value = value;
    if let Some(sep) = split {
        if !sep.is_empty() {
            value = value.split(sep.as_str()).next().unwrap_or("").to_string();
        }
    }
    if let Some(d) = delim {
        if !d.is_empty() {
            value = value.split(d.as_str()).collect::<Vec<_>>().concat();
        }
    }
    value
}

/// Extracts the barcode from SAM tags
/// (Python source: `sam_methods.get_barcode_tag`).
///
/// The split and delimiter post-processing runs in Python's order: split
/// first, then delimiter removal. The cell barcode's split defaults to `-`
/// because 10x appends a GEM well number to it.
///
/// # Errors
///
/// Returns [`PortError::Value`] when a required tag is missing, matching
/// Python, whose `read.get_tag` raises `KeyError` — which the callers in
/// [`crate::bundles`] catch to count "Read skipped, missing umi and/or cell
/// tag", exactly as Python does.
pub fn get_barcode_tag(record: &BamRecord, config: &TagConfig) -> Result<Barcode> {
    let missing = |what: &str| {
        PortError::value(format!(
            "Could not extract UMI +/- cell barcode from the read tag: {what} not present"
        ))
    };

    let raw_umi = tag_text(record, config.umi_tag)?
        .ok_or_else(|| missing(&String::from_utf8_lossy(&config.umi_tag)))?;
    let umi = apply_split_and_delim(
        raw_umi,
        config.umi_tag_split.as_ref(),
        config.umi_tag_delim.as_ref(),
    );

    let cell = if config.cell_barcode {
        let cell_tag = config
            .cell_tag
            .ok_or_else(|| PortError::value("Need to supply the --cell-tag option"))?;
        let raw = tag_text(record, cell_tag)?
            .ok_or_else(|| missing(&String::from_utf8_lossy(&cell_tag)))?;
        // Python guards the post-processing with `if cell and ...`, so an
        // empty cell barcode is passed through untouched.
        let processed = if raw.is_empty() {
            raw
        } else {
            apply_split_and_delim(
                raw,
                config.cell_tag_split.as_ref(),
                config.cell_tag_delim.as_ref(),
            )
        };
        // Python's `if cell:` before `.encode()` means an empty barcode stays
        // `""`, which is falsy and therefore never encoded — it reaches the
        // caller as the empty string, not as None.
        Some(processed.into_bytes())
    } else {
        None
    };

    Ok(Barcode {
        umi: umi.into_bytes(),
        cell,
    })
}

/// Extracts the UMI from a read name, ignoring any cell barcode
/// (Python source: `sam_methods.get_umi_read_string`).
///
/// Used by `count_tab`, whose input is a flat file of read identifiers rather
/// than an alignment.
///
/// # Errors
///
/// Returns [`PortError::Value`] when the name has no fields at all.
pub fn get_umi_read_string(read_id: &str, sep: &str) -> Result<Barcode> {
    let umi = read_id
        .split(sep)
        .last()
        .ok_or_else(|| {
            PortError::value(
                "Could not extract UMI from the read ID, please\
                 check UMI is encoded in the read name",
            )
        })?
        .as_bytes()
        .to_vec();
    Ok(Barcode { umi, cell: None })
}

/// Extracts the cell barcode and UMI from a read name
/// (Python source: `sam_methods.get_cell_umi_read_string`).
///
/// # Errors
///
/// Returns [`PortError::Value`] when the name has fewer than two fields, with
/// Python's message including the offending read ID.
pub fn get_cell_umi_read_string(read_id: &str, sep: &str) -> Result<Barcode> {
    let fields: Vec<&str> = read_id.split(sep).collect();
    if fields.len() < 2 {
        return Err(PortError::value(format!(
            "Could not extract UMI or CB from the read ID, please\
             check UMI and CB are encoded in the read name:{read_id}"
        )));
    }
    Ok(Barcode {
        umi: fields[fields.len() - 1].as_bytes().to_vec(),
        cell: Some(fields[fields.len() - 2].as_bytes().to_vec()),
    })
}

/// Extracts the barcode from a `umis`-style read name
/// (Python source: `sam_methods.get_barcode_umis`).
///
/// The `umis` toolkit writes `:UMI_<seq>:` and `:CELL_<seq>:` fields into the
/// read name; this scans the colon-separated fields for those prefixes.
///
/// # Errors
///
/// Returns [`PortError::Value`] when no `UMI_` field is present, matching
/// Python, which raises a bare `ValueError()` and catches it into the same
/// message.
pub fn get_barcode_umis(read_name: &str, cell_barcode: bool) -> Result<Barcode> {
    let mut umi = None;
    let mut cell = None;
    for element in read_name.split(':') {
        if let Some(rest) = element.strip_prefix("UMI_") {
            umi = Some(rest.as_bytes().to_vec());
        } else if cell_barcode {
            if let Some(rest) = element.strip_prefix("CELL_") {
                cell = Some(rest.as_bytes().to_vec());
            }
        }
    }
    let umi = umi.ok_or_else(|| {
        PortError::value("Could not extract UMI +/- cell barcode from the read tag")
    })?;
    Ok(Barcode { umi, cell })
}

/// How a barcode should be read out of a record
/// (Python source: the `functools.partial` bound to
/// `get_bundles.barcode_getter`, selected by `--extract-umi-method`).
///
/// Python stores a partially-applied function; an enum makes the three
/// possibilities explicit and keeps the type concrete, which a `Box<dyn Fn>`
/// would not.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BarcodeGetter {
    /// The barcode is appended to the read name.
    ReadId {
        /// Also read a cell barcode.
        cell_barcode: bool,
        /// The field separator (`--umi-separator`).
        sep: String,
    },
    /// The barcode is in SAM tags.
    Tag(TagConfig),
    /// The read name follows the `umis` toolkit's convention.
    Umis {
        /// Also read a cell barcode.
        cell_barcode: bool,
    },
}

impl BarcodeGetter {
    /// Reads the barcode out of `record`.
    ///
    /// # Errors
    ///
    /// Propagates the underlying getter's errors.
    pub fn get(&self, record: &BamRecord) -> Result<Barcode> {
        let name = || String::from_utf8_lossy(record.name().unwrap_or_default()).into_owned();
        match self {
            BarcodeGetter::ReadId { cell_barcode, sep } => {
                get_barcode_read_id(&name(), *cell_barcode, sep)
            }
            BarcodeGetter::Tag(config) => get_barcode_tag(record, config),
            BarcodeGetter::Umis { cell_barcode } => get_barcode_umis(&name(), *cell_barcode),
        }
    }
}

/// Which barcode fields a `count_tab` input line carries
/// (Python source: the `bc_getter` partial `count_tab.main` builds).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StringBarcodeGetter {
    /// `read_id[sep]UMI`.
    UmiOnly,
    /// `read_id[sep]CB[sep]UMI`.
    CellAndUmi,
}

impl StringBarcodeGetter {
    /// Reads the barcode out of a read identifier.
    ///
    /// # Errors
    ///
    /// Propagates [`get_umi_read_string`] or [`get_cell_umi_read_string`].
    pub fn get(self, read_id: &str, sep: &str) -> Result<Barcode> {
        match self {
            StringBarcodeGetter::UmiOnly => get_umi_read_string(read_id, sep),
            StringBarcodeGetter::CellAndUmi => get_cell_umi_read_string(read_id, sep),
        }
    }
}

/// One cell's UMI counts within a gene.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CellEntry {
    /// The cell barcode, `None` when the input carries none.
    pub cell: Option<Vec<u8>>,
    /// The cell's UMIs, in the order they were first seen.
    pub umis: Vec<Vec<u8>>,
    /// UMI to read count.
    pub counts: BTreeMap<Vec<u8>, u64>,
}

/// A gene's per-cell UMI counts, in the order the cells were first seen
/// (Python source: the `defaultdict(Counter)` `get_gene_count_tab` yields).
///
/// Order is not cosmetic here. `count_tab` writes one row per cell as it goes,
/// with no final sort, so the row order is the dictionary's insertion order —
/// and the UMI order within a cell decides the clustering order, which decides
/// the counts themselves when UMIs tie. A `BTreeMap` alone would silently sort
/// both, so the insertion order is carried explicitly.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CellUmiCounts {
    entries: Vec<CellEntry>,
    index: BTreeMap<Option<Vec<u8>>, usize>,
}

impl CellUmiCounts {
    /// Records one observation, creating the cell and UMI entries if new.
    pub fn observe(&mut self, cell: Option<Vec<u8>>, umi: Vec<u8>) {
        let position = if let Some(&i) = self.index.get(&cell) {
            i
        } else {
            let i = self.entries.len();
            self.index.insert(cell.clone(), i);
            self.entries.push(CellEntry {
                cell,
                ..CellEntry::default()
            });
            i
        };
        let entry = &mut self.entries[position];
        let slot = entry.counts.entry(umi.clone()).or_insert(0);
        if *slot == 0 {
            entry.umis.push(umi);
        }
        *slot += 1;
    }

    /// The cells, in first-seen order.
    #[must_use]
    pub fn cells(&self) -> &[CellEntry] {
        &self.entries
    }

    /// How many distinct cells were seen.
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether no observations were recorded.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// One cell's entry, by barcode.
    #[must_use]
    pub fn get(&self, cell: &Option<Vec<u8>>) -> Option<&CellEntry> {
        self.index.get(cell).map(|&i| &self.entries[i])
    }
}

/// Yields per-gene UMI counts from a `count_tab` input stream
/// (Python source: `sam_methods.get_gene_count_tab`).
///
/// The input is a two-column TSV of `read_id<TAB>gene`, **sorted by gene** —
/// this is a streaming group-by, so an unsorted input silently splits a gene
/// across several yields, exactly as it does in Python.
///
/// # Errors
///
/// Each item is a [`Result`]. A line without exactly two tab-separated columns
/// yields [`PortError::Assertion`] with Python's message; a read identifier
/// the getter cannot parse yields [`PortError::Value`].
pub struct GeneCountTab<R: BufRead> {
    reader: R,
    getter: StringBarcodeGetter,
    sep: String,
    gene: Option<String>,
    counts: CellUmiCounts,
    finished: bool,
}

impl<R: BufRead> GeneCountTab<R> {
    /// Reads `reader`, parsing barcodes with `getter` and `sep`.
    pub fn new(reader: R, getter: StringBarcodeGetter, sep: impl Into<String>) -> Self {
        Self {
            reader,
            getter,
            sep: sep.into(),
            gene: None,
            counts: CellUmiCounts::default(),
            finished: false,
        }
    }

    fn next_group(&mut self) -> Result<Option<(String, CellUmiCounts)>> {
        loop {
            let mut line = String::new();
            let n = self
                .reader
                .read_line(&mut line)
                .map_err(|e| PortError::io(Path::new("<count_tab input>"), e))?;
            if n == 0 {
                // Python yields the final gene unconditionally after the loop,
                // even when the file was empty — in which case it yields
                // `(None, Counter())`. That degenerate yield has no gene to
                // report, so it is simply the end of the stream here.
                self.finished = true;
                return Ok(self
                    .gene
                    .take()
                    .map(|g| (g, std::mem::take(&mut self.counts))));
            }

            let trimmed = line.trim_end_matches(['\n', '\r']).trim();
            let values: Vec<&str> = trimmed.split('\t').collect();
            if values.len() != 2 {
                return Err(PortError::assertion(format!(
                    "line: {line} does not contain 2 columns"
                )));
            }
            let (read_id, assigned_gene) = (values[0], values[1]);

            let mut finished_group = None;
            if self.gene.as_deref() != Some(assigned_gene) {
                if let Some(previous) = self.gene.take() {
                    finished_group = Some((previous, std::mem::take(&mut self.counts)));
                }
                self.gene = Some(assigned_gene.to_string());
                self.counts = CellUmiCounts::default();
            }

            let barcode = self.getter.get(read_id, &self.sep)?;
            self.counts.observe(barcode.cell, barcode.umi);

            if finished_group.is_some() {
                return Ok(finished_group);
            }
        }
    }
}

impl<R: BufRead> Iterator for GeneCountTab<R> {
    type Item = Result<(String, CellUmiCounts)>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished && self.gene.is_none() {
            return None;
        }
        match self.next_group() {
            Ok(Some(group)) => Some(Ok(group)),
            Ok(None) => None,
            Err(e) => {
                self.finished = true;
                self.gene = None;
                Some(Err(e))
            }
        }
    }
}

/// Gene to the transcripts belonging to it
/// (Python source: the `defaultdict(set)` `getMetaContig2contig` returns).
///
/// Transcripts keep the file's order; Python's `set` order varies per run,
/// and only membership is compared.
pub type MetaContigMap = IndexMap<String, IndexSet<String>>;

/// Reads a gene-to-transcript map, keeping only transcripts present in the
/// alignment (Python source: `sam_methods.getMetaContig2contig`).
///
/// Python takes the whole `AlignmentFile` and reads `bamfile.references` off
/// it; this takes the reference names directly, so the function has no
/// dependency on the alignment reader and can be tested without one.
///
/// Parsing stops at the first blank line, and `#`-prefixed lines are skipped —
/// both behaviours are Python's.
///
/// # Errors
///
/// - [`PortError::Io`] if the file cannot be read.
/// - [`PortError::Format`] for a line without exactly two tab-separated
///   columns, where Python raises `ValueError` from tuple unpacking.
pub fn get_metacontig2contig(references: &BTreeSet<String>, path: &Path) -> Result<MetaContigMap> {
    let reader = utilities::open_file(path)?;
    let mut map = MetaContigMap::new();
    for line in reader.lines() {
        let line = line.map_err(|e| PortError::io(path, e))?;
        if line.starts_with('#') {
            continue;
        }
        if line.trim().is_empty() {
            break;
        }
        let fields: Vec<&str> = line.trim_end_matches(['\n', '\r']).split('\t').collect();
        if fields.len() != 2 {
            return Err(PortError::format(format!(
                "{}: expected 2 tab-separated columns, got {}",
                path.display(),
                fields.len()
            )));
        }
        let (gene, transcript) = (fields[0], fields[1]);
        if references.contains(transcript) {
            map.entry(gene.to_string())
                .or_default()
                .insert(transcript.to_string());
        }
    }
    Ok(map)
}

/// Where a splice begins, as an offset from the alignment start
/// (Python source: `sam_methods.find_splice`).
///
/// A soft clip at the **end** of the read counts as a splice; one at the
/// start does not, and instead shifts the offset. `None` means no splice was
/// found, which Python signals by returning `False` — an `int | False` union
/// that [`Option<usize>`] models exactly.
///
/// # Errors
///
/// Returns [`PortError::Value`] for an unrecognised CIGAR operation, with
/// Python's message `"Bad Cigar operation: <op>"`.
///
/// # Examples
///
/// ```
/// use noodles_sam::alignment::record::cigar::op::Kind;
/// use umi_tools_core::sam::find_splice;
///
/// // 20M100N30M: the splice starts 20 bases in.
/// let cigar = [(Kind::Match, 20), (Kind::Skip, 100), (Kind::Match, 30)];
/// assert_eq!(find_splice(&cigar)?, Some(20));
///
/// // 50M: no splice.
/// assert_eq!(find_splice(&[(Kind::Match, 50)])?, None);
/// # Ok::<(), umi_tools_core::PortError>(())
/// ```
pub fn find_splice(cigar: &[(Kind, usize)]) -> Result<Option<usize>> {
    let mut offset = 0usize;
    let mut ops = cigar;

    // A leading soft clip is not a splice; it shifts the coordinate instead.
    if let Some(&(Kind::SoftClip, len)) = ops.first() {
        offset = len;
        ops = &ops[1..];
    }

    for &(kind, bases) in ops {
        match kind {
            // N or S: the splice starts here.
            Kind::Skip | Kind::SoftClip => return Ok(Some(offset)),
            // M, D, = or X: reference-consuming.
            Kind::Match | Kind::Deletion | Kind::SequenceMatch | Kind::SequenceMismatch => {
                offset += bases;
            }
            // I, H, P: not reference-consuming.
            Kind::Insertion | Kind::HardClip | Kind::Pad => {}
        }
    }
    Ok(None)
}

/// Where a read starts and where it is grouped
/// (Python source: the `(start, pos, is_spliced)` tuple `get_read_position`
/// returns).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadPosition {
    /// The read's leftmost coordinate, used to decide when a bundle can be
    /// flushed.
    pub start: i64,
    /// The coordinate reads are grouped on: the 5-prime end in read
    /// orientation, adjusted for soft clipping.
    pub pos: i64,
    /// Where the first splice begins, if the read is spliced.
    pub is_spliced: Option<usize>,
}

/// Computes a read's grouping coordinate, correcting for soft clipping
/// (Python source: `sam_methods.get_read_position`).
///
/// Deduplication groups reads by their 5-prime end, which is the alignment
/// start for a forward read and the alignment end for a reverse one.
/// Soft-clipped bases are added back, so a read whose first bases failed to
/// align still groups with its untrimmed copies — without that correction,
/// PCR duplicates differing only in adapter trimming would be counted twice.
///
/// # Errors
///
/// - [`PortError::Alignment`] for a record with no alignment start or end,
///   which cannot be positioned; Python would raise `TypeError` on the
///   arithmetic.
/// - [`PortError::Value`] from [`find_splice`] for an unrecognised CIGAR
///   operation.
pub fn get_read_position(record: &BamRecord, soft_clip_threshold: usize) -> Result<ReadPosition> {
    let view = record.view();
    let cigar: Vec<(Kind, usize)> = view
        .cigar()
        .iter()
        .map(|op| op.map(|op| (op.kind(), op.len())))
        .collect::<std::io::Result<_>>()
        .map_err(|e| PortError::alignment(e.to_string()))?;
    let has_skip = cigar.iter().any(|(k, _)| *k == Kind::Skip);

    let alignment_start = record
        .position()
        .ok_or_else(|| PortError::alignment("record has no alignment start"))?;

    if record.flags().is_reverse_complemented() {
        let span: usize = cigar
            .iter()
            .filter(|(k, _)| k.consumes_reference())
            .map(|(_, n)| n)
            .sum();
        let mut pos = alignment_start + i64::try_from(span).unwrap_or(i64::MAX);
        if pos == 0 {
            return Err(PortError::alignment("record has no alignment end"));
        }
        if let Some(&(Kind::SoftClip, len)) = cigar.last() {
            pos += i64::try_from(len).unwrap_or(0);
        }
        // Python tests the *first* operation here, `read.cigar[0]`.
        let leading_clip = matches!(cigar.first(), Some(&(Kind::SoftClip, len))
            if len > soft_clip_threshold);
        let is_spliced = if has_skip || leading_clip {
            // The 5-prime end of a reverse read is the right-hand end, so the
            // CIGAR is walked backwards.
            let reversed: Vec<(Kind, usize)> = cigar.iter().rev().copied().collect();
            find_splice(&reversed)?
        } else {
            None
        };
        Ok(ReadPosition {
            start: alignment_start,
            pos,
            is_spliced,
        })
    } else {
        let mut pos = alignment_start;
        if let Some(&(Kind::SoftClip, len)) = cigar.first() {
            pos -= i64::try_from(len).unwrap_or(0);
        }
        let trailing_clip = matches!(cigar.last(), Some(&(Kind::SoftClip, len))
            if len > soft_clip_threshold);
        let is_spliced = if has_skip || trailing_clip {
            find_splice(&cigar)?
        } else {
            None
        };
        Ok(ReadPosition {
            start: pos,
            pos,
            is_spliced,
        })
    }
}

/// Every read of each metacontig's contigs in turn, tagged `MC` (Python
/// source: `sam_methods.metafetcher`); the input must be indexed.
pub fn metafetch<'a>(
    reader: &'a mut OpenedAlignment,
    metacontigs: &MetaContigMap,
    tag: [u8; 2],
) -> Result<MetaFetcher<'a>> {
    if !reader.has_index()? {
        return Err(PortError::value("fetch called on bamfile without index"));
    }
    let ids: BTreeMap<String, usize> = reader
        .reference_names()
        .into_iter()
        .enumerate()
        .map(|(id, name)| (name, id))
        .collect();
    let mut plan = Vec::new();
    for (metacontig, contigs) in metacontigs {
        for contig in contigs {
            if let Some(&id) = ids.get(contig) {
                plan.push((id, metacontig.clone()));
            }
        }
    }
    Ok(MetaFetcher {
        reader,
        plan: plan.into_iter(),
        current: None,
        tag,
    })
}

/// See [`metafetch`].
pub struct MetaFetcher<'a> {
    reader: &'a mut OpenedAlignment,
    plan: std::vec::IntoIter<(usize, String)>,
    current: Option<(ReferenceScan, String)>,
    tag: [u8; 2],
}

impl Iterator for MetaFetcher<'_> {
    type Item = Result<BamRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some((scan, metacontig)) = self.current.as_mut() {
                match scan.next(self.reader) {
                    Some(Ok(mut record)) => {
                        return Some(
                            record
                                .set_tag(self.tag, TagValue::String(metacontig.as_bytes()))
                                .map(|()| record)
                                .map_err(|e| PortError::alignment(e.to_string())),
                        );
                    }
                    Some(Err(e)) => return Some(Err(e)),
                    None => self.current = None,
                }
            }
            let (id, metacontig) = self.plan.next()?;
            match self.reader.start_reference(id) {
                Ok(scan) => self.current = Some((scan, metacontig)),
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

/// Writes a deduplicated alignment together with each retained read's mate
/// (Python source: `sam_methods.TwoPassPairWriter`).
///
/// Deduplication selects read 1s, but a paired-end output is only usable if
/// each retained read 1 is accompanied by its read 2. Finding that mate by
/// seeking would be prohibitive, so the original takes two passes: it notes
/// the mates it still owes, and flushes them whenever the scan moves to a new
/// contig — by which point every mate on the previous contig has gone by.
/// Whatever is still owed at the end is found in a final full pass.
///
/// The passes reread the input; Python's per-contig `fetch(reference=chrom)`
/// becomes one forward scan of a second reader, the same on sorted input.
///
/// # A note on `close`
///
/// The Python original warns that "if `close()` is not called, at least a
/// contig's worth of mates will be missing". Here [`close`](Self::close)
/// consumes the writer, so forgetting it is a compile error rather than a
/// silently truncated file.
pub struct TwoPassPairWriter {
    source: AlignmentSource,
    sink: AlignmentSink,
    mates: Option<MateScan>,
    /// Mates still owed, keyed as Python keys them.
    owed: HashSet<MateKey>,
    chrom: Option<usize>,
}

/// Python's `(query_name, next_reference_name, next_reference_start)` of a read 1.
type MateKey = (Box<[u8]>, Option<usize>, Option<i64>);

struct MateScan {
    reader: crate::alignment::OpenedAlignment,
    next: Option<BamRecord>,
    last: usize,
}

impl TwoPassPairWriter {
    /// Creates a writer over `sink` that finds mates in `source`.
    #[must_use]
    pub fn new(source: AlignmentSource, sink: AlignmentSink) -> Self {
        Self {
            source,
            sink,
            mates: None,
            owed: HashSet::new(),
            chrom: None,
        }
    }

    /// Writes one read, flushing owed mates when the contig changes
    /// (Python source: `TwoPassPairWriter.write`).
    ///
    /// An unmapped read, or one whose mate is unmapped, is written straight
    /// through and no mate is owed for it.
    ///
    /// # Errors
    ///
    /// Returns [`PortError::Alignment`] if the record cannot be written, or a
    /// mate pass cannot read the input.
    pub fn write(&mut self, record: BamRecord, unmapped: bool) -> Result<()> {
        if unmapped || record.flags().is_mate_unmapped() {
            return self.sink.write(record);
        }

        let contig = record.reference_id();
        if self.chrom != contig {
            self.write_mates()?;
            self.chrom = contig;
        }

        self.owed.insert((
            record.name().unwrap_or_default().into(),
            record.mate_reference_id(),
            record.mate_position(),
        ));
        self.sink.write(record)
    }

    /// Writes the owed mates found on the current contig
    /// (Python source: `TwoPassPairWriter.write_mates`).
    pub fn write_mates(&mut self) -> Result<()> {
        // Python's first call scans for an empty set of mates: a no-op.
        let Some(chrom) = self.chrom else {
            return Ok(());
        };
        utilities::debug(&format!(
            "Dumping {} mates for contig {chrom}",
            self.owed.len()
        ));

        // Restart the scan only if a contig is revisited (unsorted input).
        if self.mates.as_ref().is_none_or(|scan| scan.last >= chrom) {
            self.mates = Some(MateScan {
                reader: self.source.open()?,
                next: None,
                last: chrom,
            });
        }
        let scan = self.mates.as_mut().expect("opened above");
        scan.last = chrom;
        loop {
            let record = match scan.next.take() {
                Some(record) => record,
                None => match scan.reader.read_record()? {
                    Some(record) => record,
                    None => break,
                },
            };
            match record.reference_id() {
                Some(id) if id < chrom => {}
                Some(id) if id == chrom => {
                    if Self::is_owed_mate(&mut self.owed, &record) {
                        self.sink.write(record)?;
                    }
                }
                _ => {
                    scan.next = Some(record);
                    break;
                }
            }
        }
        utilities::debug(&format!("{} mates remaining", self.owed.len()));
        Ok(())
    }

    /// Whether `record` is a mate still owed, crossing it off if so.
    fn is_owed_mate(owed: &mut HashSet<MateKey>, record: &BamRecord) -> bool {
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_mate_unmapped() || flags.is_first_segment() {
            return false;
        }
        let key: MateKey = (
            record.name().unwrap_or_default().into(),
            record.reference_id(),
            record.position(),
        );
        owed.remove(&key)
    }

    /// Flushes the last contig's mates, then searches the whole file for
    /// whatever is still owed (Python source: `TwoPassPairWriter.close`).
    ///
    /// # Errors
    ///
    /// Returns [`PortError::Alignment`] if a record cannot be read or written,
    /// or the output cannot be finalised.
    pub fn close(mut self) -> Result<()> {
        self.write_mates()?;
        self.mates = None;
        utilities::info(&format!(
            "Searching for mates for {} unmatched alignments",
            self.owed.len()
        ));
        if !self.owed.is_empty() {
            let mut reader = self.source.open()?;
            for record in reader.records(true) {
                let record = record?;
                if Self::is_owed_mate(&mut self.owed, &record) {
                    self.sink.write(record)?;
                }
            }
        }
        utilities::info(&format!("{} mates never found", self.owed.len()));
        self.sink.finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::alignment::record::data::field::Tag;
    use noodles_sam::alignment::record_buf::data::field::Value;
    use noodles_sam::alignment::RecordBuf;

    fn bam(r: &RecordBuf) -> BamRecord {
        BamRecord::from_alignment_record(&noodles_sam::Header::default(), r).unwrap()
    }
    use noodles_core::Position;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record_buf::Cigar as CigarBuf;

    fn record(name: &str) -> RecordBuf {
        let mut r = RecordBuf::default();
        *r.name_mut() = Some(name.into());
        r
    }

    fn positioned(name: &str, start: usize, ops: &[(Kind, usize)], reverse: bool) -> RecordBuf {
        let mut r = record(name);
        *r.alignment_start_mut() = Position::new(start);
        *r.cigar_mut() = ops
            .iter()
            .map(|&(k, n)| Op::new(k, n))
            .collect::<CigarBuf>();
        if reverse {
            *r.flags_mut() = Flags::REVERSE_COMPLEMENTED;
        }
        r
    }

    #[test]
    fn read_id_barcodes_match_python() {
        let bc = get_barcode_read_id("read1_ACGT", false, "_").unwrap();
        assert_eq!(bc.umi, b"ACGT");
        assert_eq!(bc.cell, None);

        let bc = get_barcode_read_id("read1_CELL_ACGT", true, "_").unwrap();
        assert_eq!(bc.umi, b"ACGT");
        assert_eq!(bc.cell.as_deref(), Some(&b"CELL"[..]));
    }

    #[test]
    fn read_id_barcodes_honour_the_separator() {
        let bc = get_barcode_read_id("read1-ACGT", false, "-").unwrap();
        assert_eq!(bc.umi, b"ACGT");
        // With the wrong separator the whole name becomes the UMI, exactly as
        // Python's split does.
        let bc = get_barcode_read_id("read1-ACGT", false, "_").unwrap();
        assert_eq!(bc.umi, b"read1-ACGT");
    }

    #[test]
    fn read_id_cell_barcode_needs_two_fields() {
        assert!(get_barcode_read_id("readonly", true, "_").is_err());
        assert!(get_barcode_read_id("readonly", false, "_").is_ok());
    }

    #[test]
    fn umi_read_string_helpers_match_python() {
        assert_eq!(get_umi_read_string("abc_ACGT", "_").unwrap().umi, b"ACGT");
        let bc = get_cell_umi_read_string("abc_CELL_ACGT", "_").unwrap();
        assert_eq!(bc.umi, b"ACGT");
        assert_eq!(bc.cell.as_deref(), Some(&b"CELL"[..]));
        assert!(get_cell_umi_read_string("solo", "_").is_err());
    }

    #[test]
    fn umis_convention_barcodes() {
        let bc = get_barcode_umis("READ:UMI_ACGT:CELL_TTTT:1", true).unwrap();
        assert_eq!(bc.umi, b"ACGT");
        assert_eq!(bc.cell.as_deref(), Some(&b"TTTT"[..]));

        let bc = get_barcode_umis("READ:UMI_ACGT:CELL_TTTT", false).unwrap();
        assert_eq!(bc.cell, None);

        assert!(get_barcode_umis("READ:CELL_TTTT", true).is_err());
    }

    #[test]
    fn tag_barcodes_apply_split_then_delim() {
        let mut r = record("read");
        r.data_mut()
            .insert(Tag::from(*b"RX"), Value::String("ACGT-1".into()));
        r.data_mut()
            .insert(Tag::from(*b"CB"), Value::String("TTTTTT-1".into()));

        let config = TagConfig {
            cell_barcode: true,
            cell_tag: Some(*b"CB"),
            ..TagConfig::default()
        };
        let bc = get_barcode_tag(&bam(&r), &config).unwrap();
        // The UMI has no split configured, so the suffix stays.
        assert_eq!(bc.umi, b"ACGT-1");
        // The cell barcode splits on the default "-", dropping the GEM well.
        assert_eq!(bc.cell.as_deref(), Some(&b"TTTTTT"[..]));

        let config = TagConfig {
            umi_tag_split: Some("-".to_string()),
            ..config
        };
        assert_eq!(get_barcode_tag(&bam(&r), &config).unwrap().umi, b"ACGT");
    }

    #[test]
    fn tag_delimiter_removes_every_occurrence() {
        let mut r = record("read");
        r.data_mut()
            .insert(Tag::from(*b"RX"), Value::String("AC-GT-TT".into()));
        let config = TagConfig {
            umi_tag_delim: Some("-".to_string()),
            ..TagConfig::default()
        };
        assert_eq!(get_barcode_tag(&bam(&r), &config).unwrap().umi, b"ACGTTT");
    }

    #[test]
    fn missing_tag_is_an_error() {
        let r = record("read");
        assert!(get_barcode_tag(&bam(&r), &TagConfig::default()).is_err());
    }

    #[test]
    fn cell_tag_without_a_tag_name_is_an_error() {
        let mut r = record("read");
        r.data_mut()
            .insert(Tag::from(*b"RX"), Value::String("ACGT".into()));
        let config = TagConfig {
            cell_barcode: true,
            cell_tag: None,
            ..TagConfig::default()
        };
        assert!(get_barcode_tag(&bam(&r), &config).is_err());
    }

    #[test]
    fn find_splice_matches_python() {
        assert_eq!(find_splice(&[(Kind::Match, 50)]).unwrap(), None);
        assert_eq!(
            find_splice(&[(Kind::SoftClip, 3), (Kind::Match, 47)]).unwrap(),
            None
        );
        assert_eq!(
            find_splice(&[(Kind::Match, 20), (Kind::Skip, 100), (Kind::Match, 30)]).unwrap(),
            Some(20)
        );
        assert_eq!(
            find_splice(&[(Kind::Match, 20), (Kind::SoftClip, 10)]).unwrap(),
            Some(20)
        );
        assert_eq!(
            find_splice(&[(Kind::HardClip, 2), (Kind::Match, 48)]).unwrap(),
            None
        );
    }

    #[test]
    fn find_splice_leading_clip_shifts_the_offset() {
        assert_eq!(
            find_splice(&[(Kind::SoftClip, 5), (Kind::Match, 20), (Kind::Skip, 100)]).unwrap(),
            Some(25)
        );
    }

    #[test]
    fn find_splice_counts_only_reference_consuming_ops() {
        assert_eq!(
            find_splice(&[
                (Kind::Match, 10),
                (Kind::Insertion, 5),
                (Kind::Deletion, 3),
                (Kind::Skip, 100)
            ])
            .unwrap(),
            Some(13)
        );
    }

    #[test]
    fn forward_read_position_uses_the_alignment_start() {
        // noodles is 1-based, pysam 0-based: position 101 -> pos 100.
        let r = positioned("r", 101, &[(Kind::Match, 50)], false);
        let p = get_read_position(&bam(&r), 4).unwrap();
        assert_eq!(p.pos, 100);
        assert_eq!(p.start, 100);
        assert_eq!(p.is_spliced, None);
    }

    #[test]
    fn forward_read_position_backs_out_a_leading_soft_clip() {
        let r = positioned("r", 101, &[(Kind::SoftClip, 6), (Kind::Match, 44)], false);
        let p = get_read_position(&bam(&r), 4).unwrap();
        assert_eq!(p.pos, 94);
        assert_eq!(p.start, 94);
    }

    #[test]
    fn reverse_read_position_uses_the_alignment_end() {
        let r = positioned("r", 101, &[(Kind::Match, 50)], true);
        let p = get_read_position(&bam(&r), 4).unwrap();
        // 1-based end 150; pysam's `aend` is the 0-based half-open end, 150.
        assert_eq!(p.pos, 150);
        assert_eq!(p.start, 100);
    }

    #[test]
    fn reverse_read_position_adds_a_trailing_soft_clip() {
        let r = positioned("r", 101, &[(Kind::Match, 44), (Kind::SoftClip, 6)], true);
        let p = get_read_position(&bam(&r), 4).unwrap();
        assert_eq!(p.pos, 150);
    }

    #[test]
    fn spliced_reads_report_the_splice_offset() {
        let r = positioned(
            "r",
            101,
            &[(Kind::Match, 20), (Kind::Skip, 100), (Kind::Match, 30)],
            false,
        );
        assert_eq!(get_read_position(&bam(&r), 4).unwrap().is_spliced, Some(20));
    }

    #[test]
    fn a_soft_clip_under_the_threshold_is_not_a_splice() {
        let r = positioned("r", 101, &[(Kind::Match, 46), (Kind::SoftClip, 4)], false);
        assert_eq!(get_read_position(&bam(&r), 4).unwrap().is_spliced, None);
        let r = positioned("r", 101, &[(Kind::Match, 45), (Kind::SoftClip, 5)], false);
        assert_eq!(get_read_position(&bam(&r), 4).unwrap().is_spliced, Some(45));
    }

    #[test]
    fn a_reverse_read_is_spliced_by_its_leading_soft_clip() {
        let r = positioned("r", 101, &[(Kind::SoftClip, 5), (Kind::Match, 45)], true);
        assert_eq!(get_read_position(&bam(&r), 4).unwrap().is_spliced, Some(45));
        let r = positioned("r", 101, &[(Kind::Match, 45), (Kind::SoftClip, 5)], true);
        assert_eq!(get_read_position(&bam(&r), 4).unwrap().is_spliced, None);
    }

    #[test]
    fn unpositioned_records_are_alignment_errors() {
        let r = record("r");
        assert!(matches!(
            get_read_position(&bam(&r), 4),
            Err(PortError::Alignment(_))
        ));
    }

    #[test]
    fn gene_count_tab_groups_by_gene() {
        let data = "r1_AAA\tgeneA\nr2_AAA\tgeneA\nr3_CCC\tgeneA\nr4_GGG\tgeneB\n";
        let groups: Vec<(String, CellUmiCounts)> =
            GeneCountTab::new(data.as_bytes(), StringBarcodeGetter::UmiOnly, "_")
                .collect::<Result<Vec<_>>>()
                .unwrap();
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].0, "geneA");
        let cell = groups[0].1.get(&None).expect("the no-cell entry");
        assert_eq!(cell.counts[&b"AAA".to_vec()], 2);
        assert_eq!(cell.counts[&b"CCC".to_vec()], 1);
        // First-seen order, not sorted order.
        assert_eq!(cell.umis, vec![b"AAA".to_vec(), b"CCC".to_vec()]);
        assert_eq!(groups[1].0, "geneB");
    }

    #[test]
    fn gene_count_tab_per_cell() {
        let data = "r1_CELL1_AAA\tgeneA\nr2_CELL2_AAA\tgeneA\n";
        let groups: Vec<(String, CellUmiCounts)> =
            GeneCountTab::new(data.as_bytes(), StringBarcodeGetter::CellAndUmi, "_")
                .collect::<Result<Vec<_>>>()
                .unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].1.len(), 2);
        let cells = groups[0].1.cells();
        assert_eq!(cells[0].cell.as_deref(), Some(&b"CELL1"[..]));
        assert_eq!(cells[0].counts[&b"AAA".to_vec()], 1);
        assert_eq!(cells[1].cell.as_deref(), Some(&b"CELL2"[..]));
    }

    #[test]
    fn gene_count_tab_rejects_malformed_lines() {
        let data = "r1_AAA\tgeneA\textra\n";
        let err = GeneCountTab::new(data.as_bytes(), StringBarcodeGetter::UmiOnly, "_")
            .collect::<Result<Vec<_>>>()
            .unwrap_err();
        assert!(matches!(err, PortError::Assertion(ref m)
            if m.contains("does not contain 2 columns")));
    }

    #[test]
    fn gene_count_tab_on_an_empty_stream_yields_nothing() {
        let groups: Vec<(String, CellUmiCounts)> =
            GeneCountTab::new(&b""[..], StringBarcodeGetter::UmiOnly, "_")
                .collect::<Result<Vec<_>>>()
                .unwrap();
        assert!(groups.is_empty());
    }

    #[test]
    fn metacontig_map_keeps_only_known_transcripts() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("map.tsv");
        std::fs::write(
            &path,
            "# header\ngeneA\ttx1\ngeneA\ttx2\ngeneB\ttx3\ngeneC\tmissing\n",
        )
        .unwrap();
        let refs: BTreeSet<String> = ["tx1", "tx2", "tx3"]
            .iter()
            .map(|s| (*s).to_string())
            .collect();
        let map = get_metacontig2contig(&refs, &path).unwrap();
        assert_eq!(map.len(), 2);
        assert_eq!(map["geneA"].len(), 2);
        assert!(!map.contains_key("geneC"));
    }

    #[test]
    fn metacontig_map_stops_at_a_blank_line() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("map.tsv");
        std::fs::write(&path, "geneA\ttx1\n\ngeneB\ttx2\n").unwrap();
        let refs: BTreeSet<String> = ["tx1", "tx2"].iter().map(|s| (*s).to_string()).collect();
        let map = get_metacontig2contig(&refs, &path).unwrap();
        assert_eq!(map.len(), 1);
    }

    #[test]
    fn metafetch_reads_genes_in_map_order_and_tags_them() {
        use crate::alignment::{open_input_alignments, AlignmentSink};
        use crate::utilities::AlignmentFormat;
        use noodles_core::Position;
        use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
        use std::num::NonZeroUsize;

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("tx.bam");
        let mut header = noodles_sam::Header::builder();
        for name in ["t1", "t2", "t3"] {
            header = header.add_reference_sequence(
                name,
                Map::<ReferenceSequence>::new(NonZeroUsize::new(100).unwrap()),
            );
        }
        let header = header.build();
        // Sorted, so the header says so, which indexing requires.
        let mut sink = AlignmentSink::open(
            Some(&path),
            AlignmentFormat::Bam,
            &header,
            None,
            None,
            true,
            None,
        )
        .unwrap();
        for (name, tid) in [("a", 0), ("b", 1), ("c", 2)] {
            let mut r = record(name);
            *r.flags_mut() = noodles_sam::alignment::record::Flags::empty();
            *r.reference_sequence_id_mut() = Some(tid);
            *r.alignment_start_mut() = Position::new(5);
            sink.write(BamRecord::from_alignment_record(&header, &r).unwrap())
                .unwrap();
        }
        sink.finish().unwrap();

        let map: MetaContigMap = [
            (
                "g2".to_string(),
                ["t3", "t1"].iter().map(|s| (*s).to_string()).collect(),
            ),
            (
                "g1".to_string(),
                ["t2"].iter().map(|s| (*s).to_string()).collect(),
            ),
        ]
        .into_iter()
        .collect();

        let mut opened = open_input_alignments(&path, AlignmentFormat::Bam, None).unwrap();
        let err = metafetch(&mut opened, &map, *b"MC").err().unwrap();
        assert!(
            matches!(err, PortError::Value(ref m) if m == "fetch called on bamfile without index")
        );

        let index = noodles_bam::fs::index(&path).unwrap();
        noodles_bam::bai::fs::write(dir.path().join("tx.bam.bai"), &index).unwrap();
        let mut opened = open_input_alignments(&path, AlignmentFormat::Bam, None).unwrap();
        let got: Vec<(String, Option<String>)> = metafetch(&mut opened, &map, *b"MC")
            .unwrap()
            .map(|r| {
                let r = r.unwrap();
                (
                    String::from_utf8_lossy(r.name().unwrap()).into_owned(),
                    tag_text(&r, *b"MC").unwrap(),
                )
            })
            .collect();
        let expect = |n: &str, g: &str| (n.to_string(), Some(g.to_string()));
        assert_eq!(
            got,
            vec![expect("c", "g2"), expect("a", "g2"), expect("b", "g1")]
        );
    }

    #[test]
    fn barcode_getter_dispatch() {
        let mut r = record("read1_CELL_ACGT");
        let g = BarcodeGetter::ReadId {
            cell_barcode: true,
            sep: "_".to_string(),
        };
        assert_eq!(g.get(&bam(&r)).unwrap().umi, b"ACGT");

        r.data_mut()
            .insert(Tag::from(*b"RX"), Value::String("TTTT".into()));
        let g = BarcodeGetter::Tag(TagConfig::default());
        assert_eq!(g.get(&bam(&r)).unwrap().umi, b"TTTT");

        let r = record("READ:UMI_GGGG");
        let g = BarcodeGetter::Umis {
            cell_barcode: false,
        };
        assert_eq!(g.get(&bam(&r)).unwrap().umi, b"GGGG");
    }
}
