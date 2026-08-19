//! Barcode extraction, quality filtering and read rewriting
//! (Python source: `umi_tools/extract_methods.py`).
//!
//! `umi_tools extract` moves the UMI (and optionally the cell barcode) out of
//! a read's sequence and into its identifier, so downstream aligners carry the
//! barcode through to the BAM. Two ways of saying where the barcode is:
//!
//! - **string** ([`ExtractMethod::String`]) — a positional mask over the read's
//!   first (or last) *n* bases, where `N` marks UMI, `C` marks cell barcode and
//!   `X` marks bases to keep on the read.
//! - **regex** ([`ExtractMethod::Regex`]) — a regular expression whose named
//!   groups mark the same three roles, allowing variable-length barcodes and
//!   linker sequences. See [`crate::pattern`].

use std::collections::{BTreeMap, BTreeSet};

use crate::error::{PortError, Result};
use crate::pattern::Pattern;
use crate::umi_methods::{QualityEncoding, Record};

/// Appends the UMI, and optionally the cell barcode, to a read identifier
/// (Python source: `extract_methods.addBarcodesToIdentifier`).
///
/// The barcodes go on the **first whitespace-separated field** only, so a
/// read's description survives intact. With a cell barcode the order is
/// `<id><sep><cell><sep><umi>`, matching what the downstream `--per-cell`
/// parsers in [`crate::sam`] expect.
///
/// An empty `cell` means "no cell barcode", exactly as Python's `if cell == "":`
/// test does — an empty string is not written as a separator-delimited field.
///
/// # Examples
///
/// ```
/// use umi_tools_core::extract::add_barcodes_to_identifier;
/// assert_eq!(add_barcodes_to_identifier("r1 desc", "AAA", "", "_"), "r1_AAA desc");
/// assert_eq!(add_barcodes_to_identifier("r1 desc", "AAA", "CCC", "_"), "r1_CCC_AAA desc");
/// ```
#[must_use]
pub fn add_barcodes_to_identifier(identifier: &str, umi: &str, cell: &str, sep: &str) -> String {
    let mut fields: Vec<&str> = identifier.split(' ').collect();
    let head = if cell.is_empty() {
        format!("{}{sep}{umi}", fields.first().copied().unwrap_or(""))
    } else {
        format!(
            "{}{sep}{cell}{sep}{umi}",
            fields.first().copied().unwrap_or("")
        )
    };
    if fields.is_empty() {
        return head;
    }
    fields[0] = &head;
    fields.join(" ")
}

/// What [`extract_seq_and_quals`] separated a read into
/// (Python source: the 4-tuple `extractSeqAndQuals` returns).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SeqQuals {
    /// The sequence left on the read.
    pub new_seq: String,
    /// The quality string left on the read.
    pub new_quals: String,
    /// Quality characters at the UMI positions.
    pub umi_quals: String,
    /// Quality characters at the cell-barcode positions.
    pub cell_quals: String,
}

/// Splits a read's sequence and qualities by base role
/// (Python source: `extract_methods.extractSeqAndQuals`).
///
/// Each position falls into exactly one branch, in Python's order: a position
/// that is neither discarded nor cell-barcode is sequence or UMI; a
/// cell-barcode position contributes only its quality; a discarded position
/// contributes nothing.
///
/// # `retain_umi` widens `umi_quals` beyond the UMI
///
/// With `retain_umi` set, Python takes the branch
///
/// ```text
/// if retain_umi:
///     new_quals += qual; new_seq += base; umi_quals += qual
/// ```
///
/// for **every** non-discarded, non-cell position — not only the UMI ones. So
/// `umi_quals` then covers the retained sequence as well. That is what the
/// original does and what any downstream quality filter therefore sees;
/// verified by running it, and mirrored here.
///
/// # Examples
///
/// ```
/// use std::collections::BTreeSet;
/// use umi_tools_core::extract::extract_seq_and_quals;
///
/// let umi: BTreeSet<usize> = [0, 1].into_iter().collect();
/// let cell: BTreeSet<usize> = [2, 3].into_iter().collect();
/// let discard: BTreeSet<usize> = [4].into_iter().collect();
///
/// let out = extract_seq_and_quals("ACGTACGT", "12345678", &umi, &cell, &discard, false);
/// assert_eq!(out.new_seq, "CGT");
/// assert_eq!(out.umi_quals, "12");
/// assert_eq!(out.cell_quals, "34");
/// ```
#[must_use]
pub fn extract_seq_and_quals(
    seq: &str,
    quals: &str,
    umi_bases: &BTreeSet<usize>,
    cell_bases: &BTreeSet<usize>,
    discard_bases: &BTreeSet<usize>,
    retain_umi: bool,
) -> SeqQuals {
    let mut out = SeqQuals::default();
    for (ix, (base, qual)) in seq.chars().zip(quals.chars()).enumerate() {
        if !discard_bases.contains(&ix) && !cell_bases.contains(&ix) {
            if retain_umi {
                out.new_quals.push(qual);
                out.new_seq.push(base);
                out.umi_quals.push(qual);
            } else if umi_bases.contains(&ix) {
                out.umi_quals.push(qual);
            } else {
                out.new_quals.push(qual);
                out.new_seq.push(base);
            }
        } else if cell_bases.contains(&ix) {
            out.cell_quals.push(qual);
        }
    }
    out
}

/// Which UMI quality scores fall below `threshold`
/// (Python source: `extract_methods.get_below_threshold`).
///
/// The score is the quality character's ASCII value minus the encoding's
/// offset, exactly as Python's `x - RANGES[quality_encoding][0]`.
///
/// # Examples
///
/// ```
/// use umi_tools_core::extract::get_below_threshold;
/// use umi_tools_core::umi_methods::QualityEncoding;
/// assert_eq!(
///     get_below_threshold("I!I", QualityEncoding::Phred33, 30),
///     vec![false, true, false]
/// );
/// ```
#[must_use]
pub fn get_below_threshold(quals: &str, encoding: QualityEncoding, threshold: i32) -> Vec<bool> {
    let offset = i32::from(encoding.offset());
    quals
        .bytes()
        .map(|q| i32::from(q) - offset < threshold)
        .collect()
}

/// Whether **any** UMI base falls below `threshold`
/// (Python source: `extract_methods.umi_below_threshold`).
///
/// An empty quality string yields `false`, matching Python's `any([])`.
#[must_use]
pub fn umi_below_threshold(quals: &str, encoding: QualityEncoding, threshold: i32) -> bool {
    get_below_threshold(quals, encoding, threshold)
        .into_iter()
        .any(|below| below)
}

/// Replaces low-quality UMI bases with `N`
/// (Python source: `extract_methods.mask_umi`).
///
/// Positions are paired with quality characters by `zip`, so a UMI longer than
/// its quality string is truncated to the shorter of the two — Python's `zip`
/// behaves the same way.
///
/// # Examples
///
/// ```
/// use umi_tools_core::extract::mask_umi;
/// use umi_tools_core::umi_methods::QualityEncoding;
/// assert_eq!(mask_umi("ACG", "I!I", QualityEncoding::Phred33, 30), "ANG");
/// ```
#[must_use]
pub fn mask_umi(umi: &str, quals: &str, encoding: QualityEncoding, threshold: i32) -> String {
    get_below_threshold(quals, encoding, threshold)
        .into_iter()
        .zip(umi.chars())
        .map(|(below, base)| if below { 'N' } else { base })
        .collect()
}

/// Which barcode roles [`extract_barcodes`] should pull out of a match
/// (Python source: the four boolean parameters of `ExtractBarcodes`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ExtractBarcodeOpts {
    /// Collect groups named `umi_*`.
    pub extract_umi: bool,
    /// Collect groups named `cell_*`.
    pub extract_cell: bool,
    /// Drop the bases covered by groups named `discard_*`.
    pub discard: bool,
    /// Keep UMI bases on the read as well as extracting them.
    pub retain_umi: bool,
}

/// What [`extract_barcodes`] pulled out of one read
/// (Python source: the 6-tuple `ExtractBarcodes` returns).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Barcodes {
    /// The concatenated cell barcode.
    pub cell_barcode: String,
    /// Quality characters at the cell-barcode positions.
    ///
    /// Python declares this in its return tuple but always returns the empty
    /// string it was initialised with — `cell_quals` from `extractSeqAndQuals`
    /// is computed and then dropped on the floor. Mirrored, so a caller
    /// reading this field sees what Python's caller sees.
    pub cell_barcode_quals: String,
    /// The concatenated UMI.
    pub umi: String,
    /// Quality characters at the UMI positions.
    pub umi_quals: String,
    /// The sequence left on the read.
    pub new_seq: String,
    /// The quality string left on the read.
    pub new_quals: String,
}

/// Pulls barcodes out of a regex match against a read
/// (Python source: `extract_methods.ExtractBarcodes`).
///
/// Groups are visited in **sorted name order** (Python's
/// `for k in sorted(list(groupdict))`), so `cell_1` precedes `cell_3` and the
/// concatenated barcode is assembled left to right by name, not by position in
/// the read. That distinction matters for the inDrop patterns, whose cell
/// barcode is split across `cell_1` and `cell_3` with a linker between.
///
/// # Errors
///
/// Returns [`PortError::Value`] when neither `extract_umi` nor `extract_cell`
/// is set, matching Python's `U.error("must set either extract_cell and/or
/// extract_umi to true")` — which exits the process there; see `05_design.md`
/// deviation D6.
pub fn extract_barcodes(
    read: &Record,
    matched: &crate::pattern::Match,
    opts: ExtractBarcodeOpts,
) -> Result<Barcodes> {
    if !opts.extract_cell && !opts.extract_umi {
        return Err(PortError::value(
            "must set either extract_cell and/or extract_umi to true",
        ));
    }

    let mut cell_barcode = String::new();
    let mut umi = String::new();
    let mut cell_bases: BTreeSet<usize> = BTreeSet::new();
    let mut umi_bases: BTreeSet<usize> = BTreeSet::new();
    let mut discard_bases: BTreeSet<usize> = BTreeSet::new();

    // `groups` is a BTreeMap, so iteration is already in sorted name order.
    for (name, &(start, end)) in &matched.groups {
        let text = &read.seq[start.min(read.seq.len())..end.min(read.seq.len())];
        if opts.extract_cell && name.starts_with("cell_") {
            cell_barcode.push_str(text);
            cell_bases.extend(start..end);
        } else if opts.extract_umi && name.starts_with("umi_") {
            umi.push_str(text);
            umi_bases.extend(start..end);
        } else if opts.discard && name.starts_with("discard_") {
            discard_bases.extend(start..end);
        }
    }

    let split = extract_seq_and_quals(
        &read.seq,
        &read.quals,
        &umi_bases,
        &cell_bases,
        &discard_bases,
        opts.retain_umi,
    );

    Ok(Barcodes {
        cell_barcode,
        // Python initialises this to "" and never reassigns it; see the field's
        // own documentation.
        cell_barcode_quals: String::new(),
        umi,
        umi_quals: split.umi_quals,
        new_seq: split.new_seq,
        new_quals: split.new_quals,
    })
}

/// How the barcode's location is specified (Python source: the
/// `--extract-method` choice, and the `self.method` string
/// `ExtractFilterAndUpdate` branches on).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ExtractMethod {
    /// A positional `N`/`C`/`X` mask.
    #[default]
    String,
    /// A regular expression with named groups.
    Regex,
}

impl std::str::FromStr for ExtractMethod {
    type Err = PortError;

    /// # Errors
    ///
    /// Returns [`PortError::Value`] for anything but `string` or `regex`.
    fn from_str(s: &str) -> Result<Self> {
        match s {
            "string" => Ok(ExtractMethod::String),
            "regex" => Ok(ExtractMethod::Regex),
            other => Err(PortError::value(format!("unknown extract method: {other}"))),
        }
    }
}

/// What to do when `--either-read` finds a UMI on both reads
/// (Python source: the `--either-read-resolve` choice).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EitherReadResolve {
    /// Throw the pair away.
    #[default]
    Discard,
    /// Keep the read whose UMI has the higher minimum quality.
    Quality,
}

impl std::str::FromStr for EitherReadResolve {
    type Err = PortError;

    /// # Errors
    ///
    /// Returns [`PortError::Value`] for anything but `discard` or `quality`,
    /// matching Python's `raise ValueError("unexpected value for
    /// either_read_resolve")`.
    fn from_str(s: &str) -> Result<Self> {
        match s {
            "discard" => Ok(EitherReadResolve::Discard),
            "quality" => Ok(EitherReadResolve::Quality),
            _ => Err(PortError::value("unexpected value for either_read_resolve")),
        }
    }
}

/// A parsed string-method barcode mask.
///
/// The mask says, per position, whether the base is UMI (`N`), cell barcode
/// (`C`) or sample sequence to keep (`X`).
#[derive(Debug, Clone, PartialEq, Eq)]
struct StringMask {
    length: usize,
    umi_bases: Vec<usize>,
    bc_bases: Vec<usize>,
    cell_bases: Vec<usize>,
}

impl StringMask {
    /// # Errors
    ///
    /// Returns [`PortError::Value`] for any character other than `N`, `X` or
    /// `C`, with Python's message.
    fn parse(pattern: &str, label: &str) -> Result<Self> {
        let stray: String = pattern
            .chars()
            .filter(|c| !matches!(c, 'N' | 'X' | 'C'))
            .collect();
        if !stray.is_empty() {
            return Err(PortError::value(format!(
                "barcode {label} ({pattern}) should only contain N/X/C characters"
            )));
        }
        Ok(Self {
            length: pattern.chars().count(),
            umi_bases: pattern
                .char_indices()
                .filter(|(_, c)| *c == 'N')
                .map(|(i, _)| i)
                .collect(),
            bc_bases: pattern
                .char_indices()
                .filter(|(_, c)| *c == 'X')
                .map(|(i, _)| i)
                .collect(),
            cell_bases: pattern
                .char_indices()
                .filter(|(_, c)| *c == 'C')
                .map(|(i, _)| i)
                .collect(),
        })
    }
}

/// Everything [`ExtractFilterAndUpdate`] needs, replacing Python's fourteen
/// positional constructor parameters.
///
/// Python's `ExtractFilterAndUpdate.__init__` takes fourteen arguments, most
/// of them booleans, and `extract.py` passes all fourteen positionally. A
/// config struct with [`Default`] makes each one nameable at the call site and
/// makes "leave it alone" the default, which is what all but two or three call
/// sites want.
#[derive(Debug, Clone, Default)]
pub struct ExtractConfig {
    /// How the barcode's location is specified.
    pub method: ExtractMethod,
    /// The read-1 pattern (an `N`/`C`/`X` mask or a regex, per `method`).
    pub pattern: Option<String>,
    /// The read-2 pattern.
    pub pattern2: Option<String>,
    /// Take the string mask from the 3-prime end instead of the 5-prime end.
    pub prime3: bool,
    /// Extract cell barcodes as well as UMIs.
    pub extract_cell: bool,
    /// The quality encoding, required by the quality filter and mask.
    pub quality_encoding: Option<QualityEncoding>,
    /// Discard reads with any UMI base below this score.
    pub quality_filter_threshold: Option<i32>,
    /// Replace UMI bases below this score with `N`.
    pub quality_filter_mask: Option<i32>,
    /// Filter UMIs against a whitelist.
    pub filter_umi_barcode: bool,
    /// Filter cell barcodes against a whitelist.
    pub filter_cell_barcode: bool,
    /// Keep the UMI on the read as well as in the identifier.
    pub retain_umi: bool,
    /// Accept a UMI found on either read of a pair.
    pub either_read: bool,
    /// How to resolve a UMI found on both reads.
    pub either_read_resolve: EitherReadResolve,
    /// The separator placed between the identifier and the barcodes.
    pub umi_separator: String,
}

impl ExtractConfig {
    /// A configuration with Python's defaults, matching
    /// `ExtractFilterAndUpdate.__init__`'s default arguments — notably
    /// `umi_separator="_"`, which [`Default`] alone would leave empty.
    #[must_use]
    pub fn new() -> Self {
        Self {
            umi_separator: "_".to_string(),
            ..Self::default()
        }
    }
}

/// The barcodes and rewritten sequences one read (or pair) yielded
/// (Python source: the 7-tuple `getBarcodes` returns).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct BarcodeValues {
    /// The cell barcode, empty when none was extracted.
    pub cell: String,
    /// The UMI.
    pub umi: String,
    /// Quality characters at the UMI positions.
    pub umi_quals: String,
    /// Read 1's sequence after extraction.
    pub new_seq: String,
    /// Read 1's qualities after extraction.
    pub new_quals: String,
    /// Read 2's sequence after extraction.
    pub new_seq2: String,
    /// Read 2's qualities after extraction.
    pub new_quals2: String,
}

/// Extracts barcodes from reads, filters them, and rewrites the reads
/// (Python source: `extract_methods.ExtractFilterAndUpdate`).
///
/// A stateful functor: it counts what it saw, which `extract` and `whitelist`
/// report at the end of a run. That is why [`process`](Self::process) takes
/// `&mut self` where Python's `__call__` looks pure.
///
/// # Examples
///
/// ```
/// use umi_tools_core::extract::{ExtractConfig, ExtractFilterAndUpdate};
/// use umi_tools_core::umi_methods::Record;
///
/// let mut extractor = ExtractFilterAndUpdate::new(ExtractConfig {
///     pattern: Some("NNNXXXXNN".to_string()),
///     ..ExtractConfig::new()
/// })?;
///
/// let mut read = Record::new("read1 d", "ACGTTTTGG", "IIIIIIIII");
/// assert!(extractor.process(&mut read, None)?);
/// assert_eq!(read.to_string(), "@read1_ACGGG d\nTTTT\n+\nIIII");
/// # Ok::<(), umi_tools_core::PortError>(())
/// ```
#[derive(Debug)]
pub struct ExtractFilterAndUpdate {
    config: ExtractConfig,
    mask: Option<StringMask>,
    mask2: Option<StringMask>,
    regex: Option<Pattern>,
    regex2: Option<Pattern>,
    read_counts: BTreeMap<String, u64>,

    /// Cell barcodes accepted without correction.
    pub cell_whitelist: Option<BTreeSet<String>>,
    /// Cell barcodes correctable to a whitelisted one.
    pub false_to_true_map: Option<BTreeMap<String, String>>,
    /// Cell barcodes to reject outright.
    pub cell_blacklist: Option<BTreeSet<String>>,
    /// UMIs accepted without correction.
    pub umi_whitelist: Option<BTreeSet<String>>,
    /// UMIs correctable to a whitelisted one; a `None` value marks a UMI
    /// within threshold of more than one whitelist entry.
    pub umi_false_to_true_map: Option<BTreeMap<String, Option<String>>>,
    /// Per-whitelisted-UMI tallies of error-free and error-corrected reads.
    pub umi_whitelist_counts: BTreeMap<String, UmiWhitelistCount>,
}

/// Error-free and error-corrected read tallies for one whitelisted UMI
/// (Python source: the inner `collections.Counter` of
/// `ReadExtractor.umi_whitelist_counts`, whose keys are `"no_error"` and
/// `"error"`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct UmiWhitelistCount {
    /// Reads whose UMI was whitelisted as observed.
    pub no_error: u64,
    /// Reads whose UMI was corrected to this one.
    pub error: u64,
}

impl ExtractFilterAndUpdate {
    /// Builds an extractor from `config`.
    ///
    /// String masks and regexes are both compiled here, matching Python, which
    /// resolves the string method's base positions in `__init__` and compiles
    /// the regexes in `Utilities.validateExtractOptions` before construction.
    ///
    /// # Errors
    ///
    /// - [`PortError::Value`] for a string mask containing anything but
    ///   `N`/`X`/`C`, with Python's message.
    /// - [`PortError::Format`] for an unparseable regex.
    pub fn new(config: ExtractConfig) -> Result<Self> {
        let (mask, mask2, regex, regex2) = match config.method {
            ExtractMethod::String => (
                config
                    .pattern
                    .as_deref()
                    .map(|p| StringMask::parse(p, "pattern"))
                    .transpose()?,
                config
                    .pattern2
                    .as_deref()
                    .map(|p| StringMask::parse(p, "pattern2"))
                    .transpose()?,
                None,
                None,
            ),
            ExtractMethod::Regex => (
                None,
                None,
                config.pattern.as_deref().map(Pattern::parse).transpose()?,
                config.pattern2.as_deref().map(Pattern::parse).transpose()?,
            ),
        };
        Ok(Self {
            config,
            mask,
            mask2,
            regex,
            regex2,
            read_counts: BTreeMap::new(),
            cell_whitelist: None,
            false_to_true_map: None,
            cell_blacklist: None,
            umi_whitelist: None,
            umi_false_to_true_map: None,
            umi_whitelist_counts: BTreeMap::new(),
        })
    }

    /// The event tallies accumulated so far
    /// (Python source: `ExtractFilterAndUpdate.getReadCounts`).
    #[must_use]
    pub fn read_counts(&self) -> &BTreeMap<String, u64> {
        &self.read_counts
    }

    fn count(&mut self, event: &str) {
        *self.read_counts.entry(event.to_string()).or_insert(0) += 1;
    }

    /// Splits a sequence at the mask boundary, 5-prime or 3-prime
    /// (Python source: `_extract_5prime` / `_extract_3prime`).
    ///
    /// Returns `(barcode_region, remaining_sequence)`. The 3-prime form takes
    /// the mask from the **end** and leaves the front, which is why the two
    /// halves come back in the same order either way.
    fn split(&self, sequence: &str, mask: &StringMask) -> (String, String) {
        let chars: Vec<char> = sequence.chars().collect();
        let n = mask.length.min(chars.len());
        if self.config.prime3 {
            let cut = chars.len() - n;
            (chars[cut..].iter().collect(), chars[..cut].iter().collect())
        } else {
            (chars[..n].iter().collect(), chars[n..].iter().collect())
        }
    }

    /// Rejoins the extracted sample bases with the remaining sequence
    /// (Python source: `_joiner_5prime` / `_joiner_3prime`).
    fn join(&self, sequence: &str, sample: &str) -> String {
        if self.config.prime3 {
            format!("{sequence}{sample}")
        } else {
            format!("{sample}{sequence}")
        }
    }

    fn pick(text: &str, positions: &[usize]) -> String {
        let chars: Vec<char> = text.chars().collect();
        positions.iter().filter_map(|&i| chars.get(i)).collect()
    }

    /// Extracts barcodes using the string mask
    /// (Python source: `_getBarcodesString`).
    fn get_barcodes_string(&self, read1: &Record, read2: Option<&Record>) -> BarcodeValues {
        let mut out = BarcodeValues::default();

        if let Some(mask) = &self.mask {
            let (bc, sequence) = self.split(&read1.seq, mask);
            let (bc_qual, seq_qual) = self.split(&read1.quals, mask);
            out.umi_quals = Self::pick(&bc_qual, &mask.umi_bases);
            out.umi = Self::pick(&bc, &mask.umi_bases);
            out.cell = Self::pick(&bc, &mask.cell_bases);
            let sample = Self::pick(&bc, &mask.bc_bases);
            let sample_qual = Self::pick(&bc_qual, &mask.bc_bases);
            out.new_seq = self.join(&sequence, &sample);
            out.new_quals = self.join(&seq_qual, &sample_qual);
        }

        if let (Some(mask2), Some(read2)) = (&self.mask2, read2) {
            let (bc, sequence) = self.split(&read2.seq, mask2);
            let (bc_qual, seq_qual) = self.split(&read2.quals, mask2);
            out.umi_quals
                .push_str(&Self::pick(&bc_qual, &mask2.umi_bases));
            out.umi.push_str(&Self::pick(&bc, &mask2.umi_bases));
            out.cell.push_str(&Self::pick(&bc, &mask2.cell_bases));
            let sample = Self::pick(&bc, &mask2.bc_bases);
            let sample_qual = Self::pick(&bc_qual, &mask2.bc_bases);
            out.new_seq2 = self.join(&sequence, &sample);
            out.new_quals2 = self.join(&seq_qual, &sample_qual);
        }

        out
    }

    /// Extracts barcodes using the regexes
    /// (Python source: `_getBarcodesRegex`).
    ///
    /// Both patterns are tested before either is used, so a read pair is not
    /// half-extracted when only one side matches — the comment in the Python
    /// source says as much.
    ///
    /// # Errors
    ///
    /// Propagates [`extract_barcodes`]' errors.
    fn get_barcodes_regex(
        &mut self,
        read1: &Record,
        read2: Option<&Record>,
    ) -> Result<Option<BarcodeValues>> {
        let opts = ExtractBarcodeOpts {
            extract_umi: true,
            extract_cell: self.config.extract_cell,
            discard: true,
            retain_umi: self.config.retain_umi,
        };

        let match1 = self.regex.as_ref().map(|p| p.match_at_start(&read1.seq));
        if let Some(m) = &match1 {
            if m.is_none() {
                self.count("regex does not match read1");
                if !self.config.either_read {
                    return Ok(None);
                }
            } else {
                self.count("regex matches read1");
            }
        }
        let match1 = match1.flatten();

        let match2 = match (read2, self.regex2.as_ref()) {
            (Some(r2), Some(p2)) => {
                let m = p2.match_at_start(&r2.seq);
                if m.is_none() {
                    self.count("regex does not match read2");
                    if self.config.either_read {
                        if match1.is_none() {
                            return Ok(None);
                        }
                    } else {
                        return Ok(None);
                    }
                } else {
                    self.count("regex matches read2");
                }
                m
            }
            _ => None,
        };

        if self.config.either_read {
            if match1.is_some()
                && match2.is_some()
                && self.config.either_read_resolve == EitherReadResolve::Discard
            {
                self.count("regex matches both. discarded");
                return Ok(None);
            }

            let b1 = match &match1 {
                Some(m) => Some(extract_barcodes(read1, m, opts)?),
                None => None,
            };
            let b2 = match (&match2, read2) {
                (Some(m), Some(r2)) => Some(extract_barcodes(r2, m, opts)?),
                _ => None,
            };

            return Ok(match (b1, b2) {
                (Some(b1), Some(b2)) => {
                    let Some(encoding) = self.config.quality_encoding else {
                        return Err(PortError::value(
                            "--either-read-resolve=quality requires --quality-encoding",
                        ));
                    };
                    let offset = i32::from(encoding.offset());
                    let min_of =
                        |q: &str| -> Option<i32> { q.bytes().map(|c| i32::from(c) - offset).min() };
                    let (Some(q1), Some(q2)) = (min_of(&b1.umi_quals), min_of(&b2.umi_quals))
                    else {
                        return Err(PortError::value(
                            "cannot compare read qualities: a UMI has no quality scores",
                        ));
                    };
                    // Note the asymmetry Python documents: the UMI is chosen by
                    // whichever read had the better quality, but the sequences
                    // of *both* reads are the post-extraction ones.
                    if q1 >= q2 {
                        self.count("regex matches both. read1 used");
                        Some(BarcodeValues {
                            cell: b1.cell_barcode,
                            umi: b1.umi,
                            umi_quals: b1.umi_quals,
                            new_seq: b1.new_seq,
                            new_quals: b1.new_quals,
                            new_seq2: b2.new_seq,
                            new_quals2: b2.new_quals,
                        })
                    } else {
                        self.count("regex matches both. read2 used");
                        Some(BarcodeValues {
                            cell: b2.cell_barcode,
                            umi: b2.umi,
                            umi_quals: b2.umi_quals,
                            new_seq: b1.new_seq,
                            new_quals: b1.new_quals,
                            new_seq2: b2.new_seq,
                            new_quals2: b2.new_quals,
                        })
                    }
                }
                (Some(b1), None) => {
                    self.count("regex only matches read1");
                    Some(BarcodeValues {
                        cell: b1.cell_barcode,
                        umi: b1.umi,
                        umi_quals: b1.umi_quals,
                        new_seq: b1.new_seq,
                        new_quals: b1.new_quals,
                        new_seq2: String::new(),
                        new_quals2: String::new(),
                    })
                }
                (None, Some(b2)) => {
                    self.count("regex only matches read2");
                    Some(BarcodeValues {
                        cell: b2.cell_barcode,
                        umi: b2.umi,
                        umi_quals: b2.umi_quals,
                        new_seq: String::new(),
                        new_quals: String::new(),
                        new_seq2: b2.new_seq,
                        new_quals2: b2.new_quals,
                    })
                }
                (None, None) => None,
            });
        }

        // The ordinary path: extract from whichever patterns were supplied.
        let mut out = BarcodeValues::default();
        if self.regex.is_some() {
            let Some(m) = &match1 else { return Ok(None) };
            let b = extract_barcodes(read1, m, opts)?;
            out.cell = b.cell_barcode;
            out.umi = b.umi;
            out.umi_quals = b.umi_quals;
            out.new_seq = b.new_seq;
            out.new_quals = b.new_quals;
        }
        if let (Some(r2), Some(_), Some(m2)) = (read2, self.regex2.as_ref(), &match2) {
            let b2 = extract_barcodes(r2, m2, opts)?;
            out.cell.push_str(&b2.cell_barcode);
            out.umi.push_str(&b2.umi);
            out.umi_quals.push_str(&b2.umi_quals);
            out.new_seq2 = b2.new_seq;
            out.new_quals2 = b2.new_quals;
        }
        Ok(Some(out))
    }

    /// Extracts barcodes from a read or read pair
    /// (Python source: the `getBarcodes` dispatch attribute).
    ///
    /// Returns `None` when the read should be discarded — for the regex
    /// method, when a required pattern did not match.
    ///
    /// # Errors
    ///
    /// Propagates [`extract_barcodes`]' errors.
    pub fn get_barcodes(
        &mut self,
        read1: &Record,
        read2: Option<&Record>,
    ) -> Result<Option<BarcodeValues>> {
        match self.config.method {
            ExtractMethod::String => Ok(Some(self.get_barcodes_string(read1, read2))),
            ExtractMethod::Regex => self.get_barcodes_regex(read1, read2),
        }
    }

    /// Extracts only the cell barcode, without rewriting the read
    /// (Python source: the `getCellBarcode` dispatch attribute).
    ///
    /// Used by `whitelist`, which needs to tally cell barcodes and has no
    /// interest in the rewritten sequence.
    ///
    /// # Errors
    ///
    /// Propagates [`extract_barcodes`]' errors.
    pub fn get_cell_barcode(
        &mut self,
        read1: &Record,
        read2: Option<&Record>,
    ) -> Result<Option<String>> {
        match self.config.method {
            ExtractMethod::String => {
                let mut cell = String::new();
                if let Some(mask) = &self.mask {
                    let (bc, _) = self.split(&read1.seq, mask);
                    cell.push_str(&Self::pick(&bc, &mask.cell_bases));
                }
                if let (Some(mask2), Some(read2)) = (&self.mask2, read2) {
                    let (bc, _) = self.split(&read2.seq, mask2);
                    cell.push_str(&Self::pick(&bc, &mask2.cell_bases));
                }
                Ok(Some(cell))
            }
            ExtractMethod::Regex => {
                let m1 = self.regex.as_ref().map(|p| p.match_at_start(&read1.seq));
                let m2 = match (self.regex2.as_ref(), read2) {
                    (Some(p2), Some(r2)) => Some(p2.match_at_start(&r2.seq)),
                    _ => None,
                };
                // Python: `if not ((self.pattern and not match1) or
                //                  (self.pattern2 and not match2))`
                if matches!(m1, Some(None)) || matches!(m2, Some(None)) {
                    return Ok(None);
                }
                let opts = ExtractBarcodeOpts {
                    extract_cell: true,
                    extract_umi: false,
                    ..ExtractBarcodeOpts::default()
                };
                let mut cell = String::new();
                if let Some(Some(m)) = &m1 {
                    cell.push_str(&extract_barcodes(read1, m, opts)?.cell_barcode);
                }
                if let (Some(Some(m)), Some(r2)) = (&m2, read2) {
                    cell.push_str(&extract_barcodes(r2, m, opts)?.cell_barcode);
                }
                Ok(Some(cell))
            }
        }
    }

    /// Whether the UMI should be discarded for low quality
    /// (Python source: `ExtractFilterAndUpdate.filterQuality`).
    ///
    /// # Errors
    ///
    /// Returns [`PortError::Value`] when no quality encoding was configured,
    /// where Python raises `KeyError: None` from `RANGES[None]`.
    pub fn filter_quality(&mut self, umi_quals: &str) -> Result<bool> {
        let Some(threshold) = self.config.quality_filter_threshold else {
            return Ok(false);
        };
        let Some(encoding) = self.config.quality_encoding else {
            return Err(PortError::value(
                "must provide a quality encoding (--quality-encoding) to filter UMIs by quality",
            ));
        };
        if umi_below_threshold(umi_quals, encoding, threshold) {
            self.count("filtered: umi quality");
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Masks low-quality UMI bases with `N`
    /// (Python source: `ExtractFilterAndUpdate.maskQuality`).
    ///
    /// # Errors
    ///
    /// Returns [`PortError::Value`] when no quality encoding was configured.
    pub fn mask_quality(&mut self, umi: &str, umi_quals: &str) -> Result<String> {
        let Some(threshold) = self.config.quality_filter_mask else {
            return Ok(umi.to_string());
        };
        let Some(encoding) = self.config.quality_encoding else {
            return Err(PortError::value(
                "must provide a quality encoding (--quality-encoding) to mask low quality bases",
            ));
        };
        let masked = mask_umi(umi, umi_quals, encoding, threshold);
        if masked != umi {
            self.count("UMI masked");
        }
        Ok(masked)
    }

    /// Filters a cell barcode against the whitelist, correcting it if possible
    /// (Python source: `ExtractFilterAndUpdate.filterCellBarcode`).
    ///
    /// Returns `None` when the read should be discarded. The blacklist is
    /// checked **twice**, before and after correction, so a barcode that
    /// corrects *into* the blacklist is still rejected — Python does the same,
    /// and the second check's counter name (`"Cell barcode corrected to
    /// barcode blacklist"`) says why.
    #[must_use]
    pub fn filter_cell_barcode(&mut self, cell: &str) -> Option<String> {
        if self
            .cell_blacklist
            .as_ref()
            .is_some_and(|b| b.contains(cell))
        {
            self.count("Cell barcode in blacklist");
            return None;
        }

        let mut cell = cell.to_string();
        let whitelisted = self
            .cell_whitelist
            .as_ref()
            .is_some_and(|w| w.contains(&cell));

        if !whitelisted {
            if let Some(map) = self.false_to_true_map.as_ref() {
                if let Some(corrected) = map.get(&cell) {
                    cell = corrected.clone();
                    self.count("False cell barcode. Error-corrected");
                } else {
                    self.count("Filtered cell barcode. Not correctable");
                    return None;
                }
            } else {
                self.count("Filtered cell barcode");
                return None;
            }
        }

        if self
            .cell_blacklist
            .as_ref()
            .is_some_and(|b| b.contains(&cell))
        {
            self.count("Cell barcode corrected to barcode blacklist");
            return None;
        }
        Some(cell)
    }

    /// Filters a UMI against the whitelist, correcting it if possible
    /// (Python source: `ExtractFilterAndUpdate.filterUMIBarcode`).
    ///
    /// Returns `None` when the read should be discarded. The three
    /// not-correctable cases are counted separately, as Python counts them:
    /// within threshold of more than one whitelist entry, within threshold of
    /// none, and no correction map at all.
    #[must_use]
    pub fn filter_umi_barcode(&mut self, umi: &str) -> Option<String> {
        if self.umi_whitelist.as_ref().is_some_and(|w| w.contains(umi)) {
            self.umi_whitelist_counts
                .entry(umi.to_string())
                .or_default()
                .no_error += 1;
            return Some(umi.to_string());
        }

        if let Some(map) = self.umi_false_to_true_map.as_ref() {
            match map.get(umi) {
                Some(Some(corrected)) => {
                    let corrected = corrected.clone();
                    self.count("False UMI barcode. Error-corrected");
                    self.umi_whitelist_counts
                        .entry(corrected.clone())
                        .or_default()
                        .error += 1;
                    Some(corrected)
                }
                Some(None) => {
                    self.count(
                        "False UMI barcode. Not correctable - within threshold \
                     to more than one whitelisted UMI",
                    );
                    None
                }
                None => {
                    self.count(
                        "False UMI barcode. Not correctable - not within \
                     threshold to whitelisted UMI",
                    );
                    None
                }
            }
        } else {
            self.count("Filtered umi barcode");
            None
        }
    }

    /// Extracts, filters and rewrites a read or read pair
    /// (Python source: `ExtractFilterAndUpdate.__call__`).
    ///
    /// Returns `true` when the read was kept and rewritten in place, `false`
    /// when it was filtered out. Python instead returns the mutated read (or
    /// tuple of reads), or `None`; rewriting in place and returning a verdict
    /// avoids cloning every record and makes the discard case a `bool` rather
    /// than an `Option` the caller must destructure.
    ///
    /// # Errors
    ///
    /// - [`PortError::Value`] when a read is shorter than its string mask,
    ///   with Python's message. Python added this check deliberately (issue
    ///   #424) so the failure is a clear error rather than an out-of-bounds
    ///   surprise; regex patterns get no such check, because a too-short read
    ///   simply fails to match.
    /// - Whatever [`get_barcodes`](Self::get_barcodes),
    ///   [`filter_quality`](Self::filter_quality) or
    ///   [`mask_quality`](Self::mask_quality) return.
    pub fn process(&mut self, read1: &mut Record, read2: Option<&mut Record>) -> Result<bool> {
        if self.config.method == ExtractMethod::String {
            if let Some(mask) = &self.mask {
                if read1.seq.chars().count() < mask.length {
                    return Err(PortError::value(format!(
                        "Read sequence: {} is shorter than pattern: {}",
                        read1.seq,
                        self.config.pattern.as_deref().unwrap_or("")
                    )));
                }
            }
            if let (Some(r2), Some(mask2)) = (read2.as_deref(), &self.mask2) {
                if r2.seq.chars().count() < mask2.length {
                    return Err(PortError::value(format!(
                        "Read2 sequence: {} is shorter than pattern2: {}",
                        r2.seq,
                        self.config.pattern2.as_deref().unwrap_or("")
                    )));
                }
            }
        }

        self.count("Input Reads");

        let Some(values) = self.get_barcodes(read1, read2.as_deref())? else {
            return Ok(false);
        };
        let BarcodeValues {
            mut cell,
            mut umi,
            umi_quals,
            new_seq,
            new_quals,
            new_seq2,
            new_quals2,
        } = values;

        if self.filter_quality(&umi_quals)? {
            return Ok(false);
        }
        umi = self.mask_quality(&umi, &umi_quals)?;

        if self.config.filter_umi_barcode {
            let Some(filtered) = self.filter_umi_barcode(&umi) else {
                return Ok(false);
            };
            umi = filtered;
        }
        if self.config.filter_cell_barcode {
            let Some(filtered) = self.filter_cell_barcode(&cell) else {
                return Ok(false);
            };
            cell = filtered;
        }

        self.count("Reads output");
        let sep = &self.config.umi_separator;

        if self.config.either_read {
            let id1 = add_barcodes_to_identifier(&read1.identifier, &umi, &cell, sep);
            read1.identifier.clone_from(&id1);
            if let Some(r2) = read2 {
                // Python computes `new_identifier2` here and then assigns
                // `new_identifier` — a copy-paste slip. Both are built from the
                // same `umi` and `cell`, so the value is identical either way;
                // mirrored as written. See 05_design.md § Intentional deviations.
                r2.identifier = id1;

                if new_seq2.is_empty() && new_quals2.is_empty() {
                    read1.seq = new_seq;
                    read1.quals = new_quals;
                } else if new_seq.is_empty() && new_quals.is_empty() {
                    r2.seq = new_seq2;
                    r2.quals = new_quals2;
                } else {
                    read1.seq = new_seq;
                    read1.quals = new_quals;
                    r2.seq = new_seq2;
                    r2.quals = new_quals2;
                }
            }
            return Ok(true);
        }

        read1.identifier = add_barcodes_to_identifier(&read1.identifier, &umi, &cell, sep);
        if self.config.pattern.is_some() {
            read1.seq = new_seq;
            read1.quals = new_quals;
        }
        if let Some(r2) = read2 {
            r2.identifier = add_barcodes_to_identifier(&r2.identifier, &umi, &cell, sep);
            if self.config.pattern2.is_some() {
                r2.seq = new_seq2;
                r2.quals = new_quals2;
            }
        }
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn string_extractor(pattern: &str) -> ExtractFilterAndUpdate {
        ExtractFilterAndUpdate::new(ExtractConfig {
            pattern: Some(pattern.to_string()),
            ..ExtractConfig::new()
        })
        .expect("config is valid")
    }

    fn regex_extractor(pattern: &str) -> ExtractFilterAndUpdate {
        ExtractFilterAndUpdate::new(ExtractConfig {
            method: ExtractMethod::Regex,
            pattern: Some(pattern.to_string()),
            ..ExtractConfig::new()
        })
        .expect("config is valid")
    }

    #[test]
    fn identifier_rewriting_matches_python() {
        assert_eq!(
            add_barcodes_to_identifier("r1 desc", "AAA", "", "_"),
            "r1_AAA desc"
        );
        assert_eq!(
            add_barcodes_to_identifier("r1 desc", "AAA", "CCC", "_"),
            "r1_CCC_AAA desc"
        );
        assert_eq!(add_barcodes_to_identifier("r1", "AAA", "", "-"), "r1-AAA");
        assert_eq!(
            add_barcodes_to_identifier("r1 a b", "U", "", "_"),
            "r1_U a b"
        );
    }

    #[test]
    fn seq_and_quals_split_matches_python() {
        let umi: BTreeSet<usize> = [0, 1].into_iter().collect();
        let cell: BTreeSet<usize> = [2, 3].into_iter().collect();
        let discard: BTreeSet<usize> = [4].into_iter().collect();

        let out = extract_seq_and_quals("ACGTACGT", "12345678", &umi, &cell, &discard, false);
        assert_eq!(out.new_seq, "CGT");
        assert_eq!(out.new_quals, "678");
        assert_eq!(out.umi_quals, "12");
        assert_eq!(out.cell_quals, "34");
    }

    #[test]
    fn retain_umi_keeps_umi_bases_on_the_read() {
        let umi: BTreeSet<usize> = [0, 1].into_iter().collect();
        let cell: BTreeSet<usize> = [2, 3].into_iter().collect();
        let discard: BTreeSet<usize> = [4].into_iter().collect();

        let out = extract_seq_and_quals("ACGTACGT", "12345678", &umi, &cell, &discard, true);
        assert_eq!(out.new_seq, "ACCGT");
        assert_eq!(out.new_quals, "12678");
        // Python's retain_umi branch appends to umi_quals for every retained
        // position, not only the UMI ones; see the function's documentation.
        assert_eq!(out.umi_quals, "12678");
        assert_eq!(out.cell_quals, "34");
    }

    #[test]
    fn quality_helpers_match_python() {
        assert_eq!(
            get_below_threshold("I!I", QualityEncoding::Phred33, 30),
            vec![false, true, false]
        );
        assert!(umi_below_threshold("I!I", QualityEncoding::Phred33, 30));
        assert!(!umi_below_threshold("III", QualityEncoding::Phred33, 30));
        assert!(!umi_below_threshold("", QualityEncoding::Phred33, 30));
        assert_eq!(mask_umi("ACG", "I!I", QualityEncoding::Phred33, 30), "ANG");
        assert_eq!(mask_umi("ACG", "III", QualityEncoding::Phred33, 30), "ACG");
    }

    #[test]
    fn quality_offsets_differ_per_encoding() {
        // 'I' is 73; phred33 -> 40, phred64 -> 9, solexa -> 14.
        assert!(!umi_below_threshold("I", QualityEncoding::Phred33, 30));
        assert!(umi_below_threshold("I", QualityEncoding::Phred64, 30));
        assert!(umi_below_threshold("I", QualityEncoding::Solexa, 30));
    }

    #[test]
    fn string_extraction_matches_python() {
        let mut e = string_extractor("NNNXXXXNN");
        let mut read = Record::new("read1 d", "ACGTTTTGG", "IIIIIIIII");
        assert!(e.process(&mut read, None).unwrap());
        assert_eq!(read.to_string(), "@read1_ACGGG d\nTTTT\n+\nIIII");
        assert_eq!(e.read_counts()["Input Reads"], 1);
        assert_eq!(e.read_counts()["Reads output"], 1);
    }

    #[test]
    fn three_prime_string_extraction_matches_python() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            pattern: Some("NNN".to_string()),
            prime3: true,
            ..ExtractConfig::new()
        })
        .unwrap();
        let mut read = Record::new("r", "ACGTTTT", "IIIIIII");
        assert!(e.process(&mut read, None).unwrap());
        assert_eq!(read.to_string(), "@r_TTT\nACGT\n+\nIIII");
    }

    #[test]
    fn cell_barcode_extraction_matches_python() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            pattern: Some("CCCCCCNNNNNNNNNN".to_string()),
            extract_cell: true,
            ..ExtractConfig::new()
        })
        .unwrap();
        let mut read = Record::new("r", "CCCCCCUUUUUUUUUUAAAA", "IIIIIIIIIIIIIIIIIIII");
        assert!(e.process(&mut read, None).unwrap());
        assert_eq!(read.identifier, "r_CCCCCC_UUUUUUUUUU");
        assert_eq!(read.seq, "AAAA");
    }

    #[test]
    fn regex_extraction_matches_python() {
        let mut e = regex_extractor("^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})");
        let mut read = Record::new("read1 d", "ACGTTTTGG", "IIIIIIIII");
        assert!(e.process(&mut read, None).unwrap());
        assert_eq!(read.to_string(), "@read1_ACGGG d\nTTTT\n+\nIIII");
    }

    #[test]
    fn regex_retain_umi_keeps_the_full_sequence() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            method: ExtractMethod::Regex,
            pattern: Some("^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})".to_string()),
            retain_umi: true,
            ..ExtractConfig::new()
        })
        .unwrap();
        let mut read = Record::new("r", "ACGTTTTGG", "IIIIIIIII");
        assert!(e.process(&mut read, None).unwrap());
        assert_eq!(read.to_string(), "@r_ACGGG\nACGTTTTGG\n+\nIIIIIIIII");
    }

    #[test]
    fn regex_that_does_not_match_discards_the_read() {
        let mut e = regex_extractor("^(?P<umi_1>ACGT)GGGG");
        let mut read = Record::new("r", "TTTTTTTT", "IIIIIIII");
        assert!(!e.process(&mut read, None).unwrap());
        assert_eq!(e.read_counts()["regex does not match read1"], 1);
        assert!(!e.read_counts().contains_key("Reads output"));
    }

    #[test]
    fn read_shorter_than_the_string_pattern_is_an_error() {
        let mut e = string_extractor("NNNXXXXNN");
        let mut read = Record::new("r", "AC", "II");
        let err = e.process(&mut read, None).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m == "Read sequence: AC is shorter than pattern: NNNXXXXNN"));
    }

    #[test]
    fn read2_shorter_than_pattern2_is_an_error() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            pattern: Some("NN".to_string()),
            pattern2: Some("NNNN".to_string()),
            ..ExtractConfig::new()
        })
        .unwrap();
        let mut r1 = Record::new("r", "ACGT", "IIII");
        let mut r2 = Record::new("r", "AC", "II");
        let err = e.process(&mut r1, Some(&mut r2)).unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.starts_with("Read2 sequence: AC is shorter than pattern2:")));
    }

    #[test]
    fn invalid_string_pattern_is_rejected() {
        let err = ExtractFilterAndUpdate::new(ExtractConfig {
            pattern: Some("NNZ".to_string()),
            ..ExtractConfig::new()
        })
        .unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m == "barcode pattern (NNZ) should only contain N/X/C characters"));
    }

    #[test]
    fn invalid_pattern2_names_the_second_pattern() {
        let err = ExtractFilterAndUpdate::new(ExtractConfig {
            pattern: Some("NN".to_string()),
            pattern2: Some("NQ".to_string()),
            ..ExtractConfig::new()
        })
        .unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m == "barcode pattern2 (NQ) should only contain N/X/C characters"));
    }

    #[test]
    fn extract_barcodes_needs_a_role_to_extract() {
        let read = Record::new("r", "ACGT", "IIII");
        let p = Pattern::parse("(?P<umi_1>.{2})").unwrap();
        let m = p.match_at_start(&read.seq).unwrap();
        let err = extract_barcodes(&read, &m, ExtractBarcodeOpts::default()).unwrap_err();
        assert!(matches!(err, PortError::Value(ref msg)
            if msg.contains("extract_cell and/or extract_umi")));
    }

    #[test]
    fn quality_filter_discards_low_quality_umis() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            pattern: Some("NNN".to_string()),
            quality_encoding: Some(QualityEncoding::Phred33),
            quality_filter_threshold: Some(30),
            ..ExtractConfig::new()
        })
        .unwrap();
        let mut good = Record::new("r", "ACGT", "IIII");
        assert!(e.process(&mut good, None).unwrap());
        let mut bad = Record::new("r", "ACGT", "!!!!");
        assert!(!e.process(&mut bad, None).unwrap());
        assert_eq!(e.read_counts()["filtered: umi quality"], 1);
    }

    #[test]
    fn quality_mask_replaces_low_quality_bases() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            pattern: Some("NNN".to_string()),
            quality_encoding: Some(QualityEncoding::Phred33),
            quality_filter_mask: Some(30),
            ..ExtractConfig::new()
        })
        .unwrap();
        let mut read = Record::new("r", "ACGT", "I!IJ");
        assert!(e.process(&mut read, None).unwrap());
        assert_eq!(read.identifier, "r_ANG");
        assert_eq!(e.read_counts()["UMI masked"], 1);
    }

    #[test]
    fn quality_filter_without_an_encoding_is_an_error() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            pattern: Some("NNN".to_string()),
            quality_filter_threshold: Some(30),
            ..ExtractConfig::new()
        })
        .unwrap();
        let mut read = Record::new("r", "ACGT", "IIII");
        assert!(e.process(&mut read, None).is_err());
    }

    #[test]
    fn cell_barcode_filtering_uses_whitelist_then_correction() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            pattern: Some("CCCNNN".to_string()),
            extract_cell: true,
            filter_cell_barcode: true,
            ..ExtractConfig::new()
        })
        .unwrap();
        e.cell_whitelist = Some(["AAA".to_string()].into_iter().collect());
        e.false_to_true_map = Some(
            [("AAT".to_string(), "AAA".to_string())]
                .into_iter()
                .collect(),
        );

        let mut good = Record::new("r", "AAACGT", "IIIIII");
        assert!(e.process(&mut good, None).unwrap());
        assert_eq!(good.identifier, "r_AAA_CGT");

        let mut correctable = Record::new("r", "AATCGT", "IIIIII");
        assert!(e.process(&mut correctable, None).unwrap());
        assert_eq!(correctable.identifier, "r_AAA_CGT");
        assert_eq!(e.read_counts()["False cell barcode. Error-corrected"], 1);

        let mut bad = Record::new("r", "GGGCGT", "IIIIII");
        assert!(!e.process(&mut bad, None).unwrap());
        assert_eq!(e.read_counts()["Filtered cell barcode. Not correctable"], 1);
    }

    #[test]
    fn blacklist_is_checked_before_and_after_correction() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            pattern: Some("CCCNNN".to_string()),
            extract_cell: true,
            filter_cell_barcode: true,
            ..ExtractConfig::new()
        })
        .unwrap();
        e.cell_whitelist = Some(["AAA".to_string()].into_iter().collect());
        e.false_to_true_map = Some(
            [("AAT".to_string(), "AAA".to_string())]
                .into_iter()
                .collect(),
        );
        e.cell_blacklist = Some(["AAA".to_string(), "TTT".to_string()].into_iter().collect());

        let mut direct = Record::new("r", "TTTCGT", "IIIIII");
        assert!(!e.process(&mut direct, None).unwrap());
        assert_eq!(e.read_counts()["Cell barcode in blacklist"], 1);

        // AAT is not itself blacklisted, but corrects into a blacklisted one.
        let mut corrected = Record::new("r", "AATCGT", "IIIIII");
        assert!(!e.process(&mut corrected, None).unwrap());
        assert_eq!(
            e.read_counts()["Cell barcode corrected to barcode blacklist"],
            1
        );
    }

    #[test]
    fn umi_filtering_counts_each_outcome_separately() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            pattern: Some("NNN".to_string()),
            filter_umi_barcode: true,
            ..ExtractConfig::new()
        })
        .unwrap();
        e.umi_whitelist = Some(["ACG".to_string()].into_iter().collect());
        e.umi_false_to_true_map = Some(
            [
                ("ACA".to_string(), Some("ACG".to_string())),
                ("TTT".to_string(), None),
            ]
            .into_iter()
            .collect(),
        );

        let mut exact = Record::new("r", "ACGT", "IIII");
        assert!(e.process(&mut exact, None).unwrap());
        assert_eq!(e.umi_whitelist_counts["ACG"].no_error, 1);

        let mut correctable = Record::new("r", "ACAT", "IIII");
        assert!(e.process(&mut correctable, None).unwrap());
        assert_eq!(correctable.identifier, "r_ACG");
        assert_eq!(e.umi_whitelist_counts["ACG"].error, 1);

        let mut ambiguous = Record::new("r", "TTTT", "IIII");
        assert!(!e.process(&mut ambiguous, None).unwrap());
        assert_eq!(
            e.read_counts()["False UMI barcode. Not correctable - within \
                             threshold to more than one whitelisted UMI"],
            1
        );

        let mut unknown = Record::new("r", "GGGT", "IIII");
        assert!(!e.process(&mut unknown, None).unwrap());
        assert_eq!(
            e.read_counts()["False UMI barcode. Not correctable - not within \
                             threshold to whitelisted UMI"],
            1
        );
    }

    #[test]
    fn umi_filtering_without_a_map_just_filters() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            pattern: Some("NNN".to_string()),
            filter_umi_barcode: true,
            ..ExtractConfig::new()
        })
        .unwrap();
        e.umi_whitelist = Some(["ACG".to_string()].into_iter().collect());
        let mut read = Record::new("r", "TTTT", "IIII");
        assert!(!e.process(&mut read, None).unwrap());
        assert_eq!(e.read_counts()["Filtered umi barcode"], 1);
    }

    #[test]
    fn paired_string_extraction_concatenates_both_barcodes() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            pattern: Some("NN".to_string()),
            pattern2: Some("NN".to_string()),
            ..ExtractConfig::new()
        })
        .unwrap();
        let mut r1 = Record::new("r a", "AACCC", "IIIII");
        let mut r2 = Record::new("r b", "GGTTT", "JJJJJ");
        assert!(e.process(&mut r1, Some(&mut r2)).unwrap());
        assert_eq!(r1.identifier, "r_AAGG a");
        assert_eq!(r2.identifier, "r_AAGG b");
        assert_eq!(r1.seq, "CCC");
        assert_eq!(r2.seq, "TTT");
    }

    #[test]
    fn either_read_discard_resolution_drops_double_matches() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            method: ExtractMethod::Regex,
            pattern: Some("^(?P<umi_1>.{3}).*".to_string()),
            pattern2: Some("^(?P<umi_1>.{3}).*".to_string()),
            either_read: true,
            either_read_resolve: EitherReadResolve::Discard,
            ..ExtractConfig::new()
        })
        .unwrap();
        let mut r1 = Record::new("r", "ACGTT", "IIIII");
        let mut r2 = Record::new("r", "GGGTT", "IIIII");
        assert!(!e.process(&mut r1, Some(&mut r2)).unwrap());
        assert_eq!(e.read_counts()["regex matches both. discarded"], 1);
    }

    #[test]
    fn either_read_quality_resolution_picks_the_better_umi() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            method: ExtractMethod::Regex,
            pattern: Some("^(?P<umi_1>.{3}).*".to_string()),
            pattern2: Some("^(?P<umi_1>.{3}).*".to_string()),
            either_read: true,
            either_read_resolve: EitherReadResolve::Quality,
            quality_encoding: Some(QualityEncoding::Phred33),
            ..ExtractConfig::new()
        })
        .unwrap();
        let mut r1 = Record::new("r", "ACGTT", "!!!II");
        let mut r2 = Record::new("r", "GGGTT", "JJJII");
        assert!(e.process(&mut r1, Some(&mut r2)).unwrap());
        assert_eq!(r1.identifier, "r_GGG");
        assert_eq!(e.read_counts()["regex matches both. read2 used"], 1);
    }

    #[test]
    fn either_read_single_match_is_reported() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            method: ExtractMethod::Regex,
            pattern: Some("^(?P<umi_1>ACG).*".to_string()),
            pattern2: Some("^(?P<umi_1>TTT).*".to_string()),
            either_read: true,
            ..ExtractConfig::new()
        })
        .unwrap();
        let mut r1 = Record::new("r", "ACGTT", "IIIII");
        let mut r2 = Record::new("r", "GGGTT", "IIIII");
        assert!(e.process(&mut r1, Some(&mut r2)).unwrap());
        assert_eq!(e.read_counts()["regex only matches read1"], 1);
        assert_eq!(r1.identifier, "r_ACG");
    }

    #[test]
    fn get_cell_barcode_does_not_rewrite_the_read() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            pattern: Some("CCCNNN".to_string()),
            extract_cell: true,
            ..ExtractConfig::new()
        })
        .unwrap();
        let read = Record::new("r", "AAACGT", "IIIIII");
        assert_eq!(
            e.get_cell_barcode(&read, None).unwrap().as_deref(),
            Some("AAA")
        );
        assert_eq!(read.seq, "AAACGT");
    }

    #[test]
    fn get_cell_barcode_regex_returns_none_when_unmatched() {
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            method: ExtractMethod::Regex,
            pattern: Some("^(?P<cell_1>ACG)".to_string()),
            extract_cell: true,
            ..ExtractConfig::new()
        })
        .unwrap();
        let hit = Record::new("r", "ACGTT", "IIIII");
        assert_eq!(
            e.get_cell_barcode(&hit, None).unwrap().as_deref(),
            Some("ACG")
        );
        let miss = Record::new("r", "TTTTT", "IIIII");
        assert_eq!(e.get_cell_barcode(&miss, None).unwrap(), None);
    }

    #[test]
    fn extract_method_and_resolve_parsing() {
        assert_eq!(
            "string".parse::<ExtractMethod>().unwrap(),
            ExtractMethod::String
        );
        assert_eq!(
            "regex".parse::<ExtractMethod>().unwrap(),
            ExtractMethod::Regex
        );
        assert!("other".parse::<ExtractMethod>().is_err());
        assert_eq!(
            "quality".parse::<EitherReadResolve>().unwrap(),
            EitherReadResolve::Quality
        );
        let err = "other".parse::<EitherReadResolve>().unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m == "unexpected value for either_read_resolve"));
    }

    #[test]
    fn multi_group_cell_barcode_is_assembled_in_name_order() {
        // The inDrop shape: a cell barcode split by a linker, reassembled from
        // cell_1 then cell_3 — by name, not by position.
        let mut e = ExtractFilterAndUpdate::new(ExtractConfig {
            method: ExtractMethod::Regex,
            pattern: Some(
                "(?P<cell_1>.{4})(?P<discard_2>GGGG)(?P<cell_3>.{4})(?P<umi_1>.{3}).*".to_string(),
            ),
            extract_cell: true,
            ..ExtractConfig::new()
        })
        .unwrap();
        let mut read = Record::new("r", "AAAAGGGGCCCCTTTAC", "IIIIIIIIIIIIIIIII");
        assert!(e.process(&mut read, None).unwrap());
        assert_eq!(read.identifier, "r_AAAACCCC_TTT");
        // The linker is discarded, the barcodes are removed, the rest remains.
        assert_eq!(read.seq, "AC");
    }
}
