//! FASTQ records, FASTQ iteration and UMI distance statistics
//! (Python source: `umi_tools/umi_methods.py`).
//!
//! `Record` and the FASTQ iterators are noted in the Python source as taken
//! from CGAT's `CGAT.Fastq`; they are reproduced here rather than replaced
//! with a `noodles-fastq` reader, because their error handling, their
//! whitespace stripping and their read-suffix removal are all observable in
//! `umi_tools extract`'s output and differ from the general-purpose reader's.

use std::io::BufRead;

use crate::distance::edit_distance_f64;
use crate::error::{PortError, Result};

/// Quality-score encodings and the ASCII ranges that identify them
/// (Python source: `umi_methods.RANGES`).
///
/// Each entry is `(name, (low, high))`, where a record is compatible with the
/// encoding when every quality character's ASCII value is at least `low` and
/// **strictly** less than `high`. The ranges overlap, so a record can be
/// compatible with more than one — which is why
/// [`Record::guess_format`] returns a list.
///
/// Declared in the order the Python dictionary literal declares them, which is
/// also its iteration order and therefore the order `guessFormat` returns.
pub const RANGES: [(&str, (u8, u8)); 3] = [
    ("phred33", (33, 77)),
    ("solexa", (59, 106)),
    ("phred64", (64, 106)),
];

/// The ASCII offset for a named quality encoding, or `None` if unknown.
///
/// The offset is the range's lower bound: `phred33` subtracts 33 from each
/// quality character, and so on. Used by
/// [`crate::extract::get_below_threshold`], which does exactly
/// `x - RANGES[quality_encoding][0]`.
#[must_use]
pub fn range_for(name: &str) -> Option<(u8, u8)> {
    RANGES
        .iter()
        .find(|(n, _)| *n == name)
        .map(|(_, range)| *range)
}

/// A quality-score encoding, as accepted by `--quality-encoding`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QualityEncoding {
    /// Sanger / Illumina 1.8+, offset 33.
    Phred33,
    /// Solexa / early Illumina, offset 59.
    Solexa,
    /// Illumina 1.3-1.7, offset 64.
    Phred64,
}

impl QualityEncoding {
    /// The encoding's name, as spelled on the command line and in `RANGES`.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            QualityEncoding::Phred33 => "phred33",
            QualityEncoding::Solexa => "solexa",
            QualityEncoding::Phred64 => "phred64",
        }
    }

    /// The ASCII value subtracted from a quality character to get its score.
    #[must_use]
    pub fn offset(self) -> u8 {
        match self {
            QualityEncoding::Phred33 => 33,
            QualityEncoding::Solexa => 59,
            QualityEncoding::Phred64 => 64,
        }
    }
}

impl std::str::FromStr for QualityEncoding {
    type Err = PortError;

    /// # Errors
    ///
    /// Returns [`PortError::Value`] for a name that is not one of the three
    /// keys of `RANGES`. Python would raise `KeyError` at the point of use.
    fn from_str(s: &str) -> Result<Self> {
        match s {
            "phred33" => Ok(QualityEncoding::Phred33),
            "solexa" => Ok(QualityEncoding::Solexa),
            "phred64" => Ok(QualityEncoding::Phred64),
            other => Err(PortError::value(format!(
                "unknown quality encoding: {other}"
            ))),
        }
    }
}

/// One FASTQ record (Python source: `umi_methods.Record`).
///
/// # Deviation: the discarded `entry_format` parameter
///
/// Python's `__init__` takes a fourth `entry_format` argument and then does:
///
/// ```text
/// self.identifier, self.seq, self.quals, entry_format = (
///     identifier, seq, quals, entry_format)
/// ```
///
/// The fourth target is the *local* name, not `self.format`, so the argument
/// is accepted and immediately thrown away; `vars(record)` confirms only three
/// attributes exist. It is therefore not a field here. See `05_design.md`
/// deviation D4.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    /// The identifier line without its leading `@`.
    pub identifier: String,
    /// The base calls.
    pub seq: String,
    /// The quality string, one character per base.
    pub quals: String,
}

impl Record {
    /// Creates a record from its three fields.
    #[must_use]
    pub fn new(
        identifier: impl Into<String>,
        seq: impl Into<String>,
        quals: impl Into<String>,
    ) -> Self {
        Self {
            identifier: identifier.into(),
            seq: seq.into(),
            quals: quals.into(),
        }
    }

    /// Every quality encoding compatible with this record's quality string
    /// (Python source: `Record.guessFormat`).
    ///
    /// An encoding matches when the minimum quality character is at least its
    /// lower bound and the maximum is **strictly** below its upper bound. More
    /// than one can match, because the ranges overlap; none matching means the
    /// quality string is outside every known encoding.
    ///
    /// # Errors
    ///
    /// Returns [`PortError::Value`] for an empty quality string, where Python
    /// reaches `min(c)` on an empty list and raises
    /// `ValueError: min() arg is an empty sequence`.
    ///
    /// # Examples
    ///
    /// ```
    /// use umi_tools_core::umi_methods::Record;
    /// let r = Record::new("id", "ACGT", "!!!!");
    /// assert_eq!(r.guess_format()?, vec!["phred33"]);
    /// let r = Record::new("id", "ACGT", "~~~~");
    /// assert!(r.guess_format()?.is_empty());
    /// # Ok::<(), umi_tools_core::PortError>(())
    /// ```
    pub fn guess_format(&self) -> Result<Vec<&'static str>> {
        // Python iterates `self.quals` as characters and takes `ord`, so a
        // multi-byte character would yield a code point above 255. Quality
        // strings are ASCII by the FASTQ specification; bytes are used here
        // and the ASCII assumption is stated rather than silently relied on.
        let (Some(&lo), Some(&hi)) = (
            self.quals.as_bytes().iter().min(),
            self.quals.as_bytes().iter().max(),
        ) else {
            return Err(PortError::value("min() arg is an empty sequence"));
        };
        Ok(RANGES
            .iter()
            .filter(|(_, (low, high))| lo >= *low && hi < *high)
            .map(|(name, _)| *name)
            .collect())
    }
}

impl std::fmt::Display for Record {
    /// Formats the record as four FASTQ lines **without** a trailing newline
    /// (Python source: `Record.__str__`, `"@%s\n%s\n+\n%s"`).
    ///
    /// The missing trailing newline is deliberate: every caller in the Python
    /// package writes `str(read) + "\n"`, so adding one here would double it.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "@{}\n{}\n+\n{}", self.identifier, self.seq, self.quals)
    }
}

/// Strips a `/1` or `/2` suffix from a read identifier's first field
/// (Python source: the nested `removeReadIDSuffix` in `fastqIterate`).
///
/// Only the first whitespace-separated field is considered, and the rest of
/// the line is preserved, so `"r1/1 some description"` becomes
/// `"r1 some description"`.
///
/// # Errors
///
/// Returns [`PortError::Value`] when the first field does not end in `/1` or
/// `/2`, matching Python's `ValueError: read suffix must be /1 or /2.
/// Observed: <last two characters>`.
///
/// # Examples
///
/// ```
/// use umi_tools_core::umi_methods::remove_read_id_suffix;
/// assert_eq!(remove_read_id_suffix("r1/1 desc")?, "r1 desc");
/// assert!(remove_read_id_suffix("r1 desc").is_err());
/// # Ok::<(), umi_tools_core::PortError>(())
/// ```
pub fn remove_read_id_suffix(line: &str) -> Result<String> {
    let mut components: Vec<&str> = line.split(' ').collect();
    let read_id = components.first().copied().unwrap_or("");
    let tail: String = read_id
        .chars()
        .rev()
        .take(2)
        .collect::<Vec<_>>()
        .into_iter()
        .rev()
        .collect();
    if tail != "/1" && tail != "/2" {
        return Err(PortError::value(format!(
            "read suffix must be /1 or /2. Observed: {tail}"
        )));
    }
    let trimmed = &read_id[..read_id.len() - 2];
    components[0] = trimmed;
    Ok(components.join(" "))
}

/// Iterates the records of a FASTQ stream
/// (Python source: `umi_methods.fastqIterate`).
///
/// Reproduces the original's whitespace handling exactly, which is *not*
/// uniform: the identifier line is `.strip()`ped on both sides while the
/// sequence, separator and quality lines are only `.rstrip()`ped. Leading
/// whitespace on a sequence line is therefore preserved, as in Python.
///
/// # Errors
///
/// Each item is a [`Result`]; a malformed record yields
/// [`PortError::Format`]. The Python original calls `Utilities.error()` for
/// these, which logs and calls `sys.exit(1)` from inside a generator — a
/// library must not terminate the process, so this is recoverable and
/// `crates/cli` reproduces the message and exit code. See `05_design.md`
/// deviation D6. Underlying read failures yield [`PortError::Io`].
///
/// # Examples
///
/// ```
/// use umi_tools_core::umi_methods::FastqIterator;
/// let data = b"@r1 desc\nACGT\n+\nIIII\n";
/// let records: Vec<_> = FastqIterator::new(&data[..], false)
///     .collect::<Result<Vec<_>, _>>()?;
/// assert_eq!(records[0].identifier, "r1 desc");
/// assert_eq!(records[0].seq, "ACGT");
/// # Ok::<(), umi_tools_core::PortError>(())
/// ```
pub struct FastqIterator<R: BufRead> {
    reader: R,
    remove_suffix: bool,
    finished: bool,
}

impl<R: BufRead> FastqIterator<R> {
    /// Creates an iterator over `reader`.
    ///
    /// `remove_suffix` strips a `/1`/`/2` suffix from each identifier, and
    /// makes a record whose identifier lacks one an error.
    pub fn new(reader: R, remove_suffix: bool) -> Self {
        Self {
            reader,
            remove_suffix,
            finished: false,
        }
    }

    /// Reads one line, returning `None` at end of stream.
    fn read_line(&mut self) -> Result<Option<String>> {
        let mut buf = Vec::new();
        let n = self
            .reader
            .read_until(b'\n', &mut buf)
            .map_err(|e| PortError::io(std::path::Path::new("<fastq>"), e))?;
        if n == 0 {
            return Ok(None);
        }
        // Python decodes with UTF-8 via `convert2string`; invalid bytes raise
        // `UnicodeDecodeError` there, so they are an error here too.
        String::from_utf8(buf)
            .map(Some)
            .map_err(|e| PortError::format(format!("invalid UTF-8 in FASTQ stream: {e}")))
    }

    fn next_record(&mut self) -> Result<Option<Record>> {
        // Python: `line1 = convert2string(infile.readline()).strip()`, then
        // `if not line1: break`. An all-whitespace line therefore ends the
        // stream rather than erroring.
        let Some(raw1) = self.read_line()? else {
            return Ok(None);
        };
        let mut line1 = raw1.trim().to_string();
        if line1.is_empty() {
            return Ok(None);
        }
        if !line1.starts_with('@') {
            return Err(PortError::format(format!(
                "parsing error: expected '@' in line {line1}"
            )));
        }

        let line2 = self.read_line()?.unwrap_or_default().trim_end().to_string();
        let line3 = self.read_line()?.unwrap_or_default().trim_end().to_string();
        if !line3.starts_with('+') {
            return Err(PortError::format(format!(
                "parsing error: expected '+' in line {line3}"
            )));
        }
        let line4 = self.read_line()?.unwrap_or_default().trim_end().to_string();
        if line4.is_empty() {
            return Err(PortError::format(format!("incomplete entry for {line1}")));
        }

        if self.remove_suffix {
            line1 = remove_read_id_suffix(&line1)?;
        }
        // Strip the leading '@' only after suffix removal, matching Python,
        // which yields `Record(line1[1:], line2, line4)`.
        Ok(Some(Record::new(&line1[1..], line2, line4)))
    }
}

impl<R: BufRead> Iterator for FastqIterator<R> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        match self.next_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => {
                self.finished = true;
                None
            }
            Err(e) => {
                self.finished = true;
                Some(Err(e))
            }
        }
    }
}

/// Pairs the records of two FASTQ streams, checking their identifiers agree
/// (Python source: `umi_methods.joinedFastqIterate`).
///
/// Identifiers are compared on their first whitespace-separated field only, so
/// differing descriptions do not break the pairing.
///
/// `strict` controls the response to a mismatch: when set, a mismatch is an
/// error; when clear, the second stream is advanced until it catches up. The
/// non-strict mode exists for protocols where read 1 has been filtered
/// independently of read 2.
///
/// # Errors
///
/// - [`PortError::Value`] on a mismatch in strict mode, with Python's message
///   `"\nRead pairs do not match\n<id1> != <id2>"`.
/// - [`PortError::Value`] when the second stream runs out. Python surfaces
///   this as `RuntimeError: generator raised StopIteration`, because `next()`
///   is called unguarded inside a generator — confirmed by running it.
/// - Whatever the underlying iterators yield.
pub struct JoinedFastqIterator<A, B>
where
    A: Iterator<Item = Result<Record>>,
    B: Iterator<Item = Result<Record>>,
{
    first: A,
    second: B,
    strict: bool,
    finished: bool,
}

impl<A, B> JoinedFastqIterator<A, B>
where
    A: Iterator<Item = Result<Record>>,
    B: Iterator<Item = Result<Record>>,
{
    /// Pairs `first` with `second`.
    pub fn new(first: A, second: B, strict: bool) -> Self {
        Self {
            first,
            second,
            strict,
            finished: false,
        }
    }

    /// The identifier's first whitespace-separated field
    /// (Python source: the nested `getReadID`).
    fn read_id(record: &Record) -> &str {
        record.identifier.split(' ').next().unwrap_or("")
    }

    fn next_pair(&mut self) -> Result<Option<(Record, Record)>> {
        let Some(read1) = self.first.next().transpose()? else {
            return Ok(None);
        };
        let mut read2 = self.second.next().transpose()?.ok_or_else(|| {
            PortError::value("generator raised StopIteration: read2 stream exhausted")
        })?;

        let id1 = Self::read_id(&read1).to_string();
        if !self.strict {
            while Self::read_id(&read2) != id1 {
                read2 = self.second.next().transpose()?.ok_or_else(|| {
                    PortError::value("generator raised StopIteration: read2 stream exhausted")
                })?;
            }
        }

        let id2 = Self::read_id(&read2).to_string();
        if id1 != id2 {
            return Err(PortError::value(format!(
                "\nRead pairs do not match\n{id1} != {id2}"
            )));
        }
        Ok(Some((read1, read2)))
    }
}

impl<A, B> Iterator for JoinedFastqIterator<A, B>
where
    A: Iterator<Item = Result<Record>>,
    B: Iterator<Item = Result<Record>>,
{
    type Item = Result<(Record, Record)>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        match self.next_pair() {
            Ok(Some(pair)) => Some(Ok(pair)),
            Ok(None) => {
                self.finished = true;
                None
            }
            Err(e) => {
                self.finished = true;
                Some(Err(e))
            }
        }
    }
}

/// The mean pairwise Hamming distance between a set of UMIs
/// (Python source: `umi_methods.get_average_umi_distance`).
///
/// Returns `-1.0` for a single UMI, matching Python's sentinel for "no pairs
/// to average" — chosen there because a real mean distance is never negative,
/// so the caller can distinguish it. An empty input also returns `-1.0`;
/// Python would divide by zero and produce `nan`, but only after
/// `itertools.combinations` yields nothing, which the length check already
/// intercepts for one element and which no caller reaches with zero.
///
/// Distances between unequal-length UMIs are infinite (see
/// [`crate::distance::edit_distance_f64`]), so any such pair makes the mean
/// infinite — the same value Python's `np.inf` would produce.
///
/// # Examples
///
/// ```
/// use umi_tools_core::umi_methods::get_average_umi_distance;
/// assert_eq!(get_average_umi_distance(&[b"ACGT".to_vec()]), -1.0);
/// let umis = vec![b"ACGT".to_vec(), b"ACGA".to_vec()];
/// assert_eq!(get_average_umi_distance(&umis), 1.0);
/// ```
#[must_use]
pub fn get_average_umi_distance(umis: &[Vec<u8>]) -> f64 {
    if umis.len() <= 1 {
        return -1.0;
    }
    let mut total = 0.0;
    let mut pairs = 0u64;
    for i in 0..umis.len() {
        for j in (i + 1)..umis.len() {
            total += edit_distance_f64(&umis[i], &umis[j]);
            pairs += 1;
        }
    }
    total / pairs as f64
}

/// How many UMIs are drawn per refill
/// (Python source: `self.random_fill_size = 100000  # Higher = faster, more memory`).
pub const RANDOM_FILL_SIZE: usize = 100_000;

/// Draws UMIs at random in proportion to how often they occur
/// (Python source: `umi_methods.random_read_generator`).
///
/// `dedup --output-stats` uses this to build a null distribution: the average
/// edit distance between randomly drawn UMIs is what the observed
/// within-position distance is compared against, so a position whose UMIs are
/// closer together than chance would suggest is evidence of PCR duplication.
///
/// The Python original reads the frequencies out of the BAM itself; here the
/// table is supplied by the caller, so this type has no dependency on the
/// alignment reader and can be tested without one. `crates/cli`'s `dedup`
/// builds the table while scanning the BAM it is already reading, which also
/// avoids the second full pass Python makes.
///
/// # Reproducing `NumPy`'s draws
///
/// Python draws with `numpy.random.choice(keys, size, p=probs)` after
/// `numpy.random.seed(options.random_seed)`, and `dedup`'s statistics files
/// are compared byte for byte by the package's own test suite. So this
/// reproduces `NumPy`'s *legacy* generator exactly: MT19937 seeded with
/// `init_genrand` (not `CPython`'s `init_by_array` — see
/// [`crate::mersenne::MersenneTwister::seed_numpy`]), and `choice`'s
/// inverse-transform sampling, which normalises the probabilities to a
/// cumulative distribution and takes a right-side binary search per uniform
/// draw.
///
/// The frequency table must be in the order Python's `defaultdict` holds it —
/// first-seen order — because that is the order `choice` indexes into.
#[derive(Debug, Clone)]
pub struct RandomReadGenerator {
    umis: Vec<Vec<u8>>,
    /// Normalised cumulative probabilities, for inverse-transform sampling.
    cumulative: Vec<f64>,
    rng: crate::mersenne::MersenneTwister,
    buffer: Vec<usize>,
    index: usize,
    fill_size: usize,
}

impl RandomReadGenerator {
    /// Builds a generator over `frequencies`, in first-seen order.
    ///
    /// # Errors
    ///
    /// Returns [`PortError::Value`] if `frequencies` is empty or every count
    /// is zero, where `NumPy`'s `choice` raises
    /// `ValueError: probabilities do not sum to 1`.
    pub fn new(frequencies: &[(Vec<u8>, u64)], seed: u32) -> Result<Self> {
        let total: u64 = frequencies.iter().map(|(_, c)| *c).sum();
        if frequencies.is_empty() || total == 0 {
            return Err(PortError::value(
                "cannot sample UMIs: no observed UMI counts",
            ));
        }
        crate::utilities::info(&format!("total_umis {total}"));
        crate::utilities::info(&format!("#umis {}", frequencies.len()));

        // Python normalises to probabilities first, then `choice` builds a
        // cumulative distribution and divides by its last element. Both steps
        // are reproduced, because the intermediate rounding is observable in
        // which index a borderline uniform lands on.
        let probabilities: Vec<f64> = frequencies
            .iter()
            .map(|(_, c)| *c as f64 / total as f64)
            .collect();
        let mut cumulative = Vec::with_capacity(probabilities.len());
        let mut acc = 0.0;
        for p in &probabilities {
            acc += p;
            cumulative.push(acc);
        }
        let last = *cumulative.last().expect("non-empty");
        for c in &mut cumulative {
            *c /= last;
        }

        let mut generator = Self {
            umis: frequencies.iter().map(|(u, _)| u.clone()).collect(),
            cumulative,
            rng: crate::mersenne::MersenneTwister::seed_numpy(seed),
            buffer: Vec::new(),
            index: 0,
            fill_size: RANDOM_FILL_SIZE,
        };
        generator.refill_random();
        Ok(generator)
    }

    /// Redraws the buffer of random UMI indices
    /// (Python source: `random_read_generator.refill_random`).
    ///
    /// `numpy.random.choice(a, size, p=...)` draws `size` uniforms and maps
    /// each through the cumulative distribution with a **right-side** binary
    /// search, which is what `cdf.searchsorted(uniform, side='right')` does.
    pub fn refill_random(&mut self) {
        self.buffer = (0..self.fill_size)
            .map(|_| {
                let u = self.rng.random();
                self.cumulative
                    .partition_point(|&c| c <= u)
                    .min(self.umis.len() - 1)
            })
            .collect();
        self.index = 0;
    }

    /// Draws `n` UMIs with replacement
    /// (Python source: `random_read_generator.getUmis`).
    ///
    /// Reproduces the original's buffering exactly, including that a request
    /// larger than the buffer *discards the remainder* and redraws — so the
    /// stream consumed depends on the request sizes, not only on their total.
    #[must_use]
    pub fn get_umis(&mut self, n: usize) -> Vec<Vec<u8>> {
        if n >= self.fill_size.saturating_sub(self.index) {
            if n > self.fill_size {
                self.fill_size = n * 2;
            }
            self.refill_random();
        }
        let end = (self.index + n).min(self.buffer.len());
        let drawn: Vec<Vec<u8>> = self.buffer[self.index..end]
            .iter()
            .map(|&i| self.umis[i].clone())
            .collect();
        self.index += n;
        drawn
    }

    /// The distinct UMIs available to draw from.
    #[must_use]
    pub fn umis(&self) -> &[Vec<u8>] {
        &self.umis
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn records(data: &[u8], remove_suffix: bool) -> Result<Vec<Record>> {
        FastqIterator::new(data, remove_suffix).collect()
    }

    #[test]
    fn ranges_match_python_literal() {
        assert_eq!(range_for("phred33"), Some((33, 77)));
        assert_eq!(range_for("solexa"), Some((59, 106)));
        assert_eq!(range_for("phred64"), Some((64, 106)));
        assert_eq!(range_for("nope"), None);
    }

    #[test]
    fn quality_encoding_offsets_match_range_lower_bounds() {
        for enc in [
            QualityEncoding::Phred33,
            QualityEncoding::Solexa,
            QualityEncoding::Phred64,
        ] {
            assert_eq!(range_for(enc.as_str()).unwrap().0, enc.offset());
            assert_eq!(enc.as_str().parse::<QualityEncoding>().unwrap(), enc);
        }
        assert!("phred40".parse::<QualityEncoding>().is_err());
    }

    #[test]
    fn guess_format_matches_python() {
        assert_eq!(
            Record::new("x", "ACGT", "IIII").guess_format().unwrap(),
            vec!["phred33", "solexa", "phred64"]
        );
        assert_eq!(
            Record::new("x", "ACGT", "!!!!").guess_format().unwrap(),
            vec!["phred33"]
        );
        assert_eq!(
            Record::new("x", "ACGT", "abcd").guess_format().unwrap(),
            vec!["solexa", "phred64"]
        );
        assert!(Record::new("x", "ACGT", "~~~~")
            .guess_format()
            .unwrap()
            .is_empty());
    }

    #[test]
    fn guess_format_on_empty_quals_is_a_value_error() {
        let err = Record::new("x", "", "").guess_format().unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m == "min() arg is an empty sequence"));
    }

    #[test]
    fn display_has_no_trailing_newline() {
        let r = Record::new("read1 extra", "ACGT", "IIII");
        assert_eq!(r.to_string(), "@read1 extra\nACGT\n+\nIIII");
    }

    #[test]
    fn suffix_removal_keeps_the_description() {
        assert_eq!(remove_read_id_suffix("@r1/1 desc").unwrap(), "@r1 desc");
        assert_eq!(remove_read_id_suffix("@r1/2").unwrap(), "@r1");
        let err = remove_read_id_suffix("@r1 desc").unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m == "read suffix must be /1 or /2. Observed: r1"));
    }

    #[test]
    fn iterates_a_well_formed_stream() {
        let data = b"@r1 desc\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n";
        let out = records(data, false).unwrap();
        assert_eq!(out.len(), 2);
        assert_eq!(out[0], Record::new("r1 desc", "ACGT", "IIII"));
        assert_eq!(out[1], Record::new("r2", "TTTT", "JJJJ"));
    }

    #[test]
    fn empty_stream_yields_nothing() {
        assert!(records(b"", false).unwrap().is_empty());
    }

    #[test]
    fn trailing_blank_lines_end_the_stream() {
        let data = b"@r1\nACGT\n+\nIIII\n\n";
        assert_eq!(records(data, false).unwrap().len(), 1);
    }

    #[test]
    fn suffix_removal_applies_to_every_record() {
        let data = b"@r1/1 desc\nACGT\n+\nIIII\n";
        assert_eq!(records(data, true).unwrap()[0].identifier, "r1 desc");
        assert!(records(b"@r1 desc\nACGT\n+\nIIII\n", true).is_err());
    }

    #[test]
    fn malformed_records_are_format_errors() {
        for (data, needle) in [
            (&b"r1\nAC\n+\nII\n"[..], "expected '@'"),
            (&b"@r1\nAC\nX\nII\n"[..], "expected '+'"),
            (&b"@r1\nAC\n+\n"[..], "incomplete entry"),
        ] {
            let err = records(data, false).unwrap_err();
            assert!(
                matches!(err, PortError::Format(ref m) if m.contains(needle)),
                "expected {needle}, got {err:?}"
            );
        }
    }

    #[test]
    fn iteration_stops_after_an_error() {
        let data = b"@r1\nAC\nX\nII\n@r2\nAC\n+\nII\n";
        let mut it = FastqIterator::new(&data[..], false);
        assert!(it.next().unwrap().is_err());
        assert!(it.next().is_none());
    }

    #[test]
    fn joined_iteration_pairs_matching_ids() {
        let d1 = &b"@r1 a\nAC\n+\nII\n@r2 a\nGT\n+\nJJ\n"[..];
        let d2 = &b"@r1 b\nTT\n+\nKK\n@r2 b\nCC\n+\nLL\n"[..];
        let pairs: Vec<(Record, Record)> = JoinedFastqIterator::new(
            FastqIterator::new(d1, false),
            FastqIterator::new(d2, false),
            true,
        )
        .collect::<Result<Vec<_>>>()
        .unwrap();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].0.identifier, "r1 a");
        assert_eq!(pairs[0].1.identifier, "r1 b");
    }

    #[test]
    fn strict_mismatch_uses_pythons_message() {
        let d1 = &b"@r1 a\nAC\n+\nII\n"[..];
        let d2 = &b"@rX b\nTT\n+\nKK\n"[..];
        let err: PortError = JoinedFastqIterator::new(
            FastqIterator::new(d1, false),
            FastqIterator::new(d2, false),
            true,
        )
        .collect::<Result<Vec<_>>>()
        .unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m == "\nRead pairs do not match\nr1 != rX"));
    }

    #[test]
    fn non_strict_advances_the_second_stream() {
        let d1 = &b"@r2 a\nAC\n+\nII\n"[..];
        let d2 = &b"@r1 b\nTT\n+\nKK\n@r2 b\nCC\n+\nLL\n"[..];
        let pairs: Vec<(Record, Record)> = JoinedFastqIterator::new(
            FastqIterator::new(d1, false),
            FastqIterator::new(d2, false),
            false,
        )
        .collect::<Result<Vec<_>>>()
        .unwrap();
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].1.identifier, "r2 b");
    }

    #[test]
    fn non_strict_exhaustion_is_an_error() {
        let d1 = &b"@rZ a\nAC\n+\nII\n"[..];
        let d2 = &b"@r1 b\nTT\n+\nKK\n"[..];
        let err: PortError = JoinedFastqIterator::new(
            FastqIterator::new(d1, false),
            FastqIterator::new(d2, false),
            false,
        )
        .collect::<Result<Vec<_>>>()
        .unwrap_err();
        assert!(matches!(err, PortError::Value(ref m)
            if m.contains("StopIteration")));
    }

    #[test]
    fn average_distance_matches_python() {
        assert_eq!(get_average_umi_distance(&[b"ACGT".to_vec()]), -1.0);
        assert_eq!(get_average_umi_distance(&[]), -1.0);
        let umis = vec![b"ACGT".to_vec(), b"ACGA".to_vec(), b"TTTT".to_vec()];
        // pairs: (ACGT,ACGA)=1, (ACGT,TTTT)=3, (ACGA,TTTT)=4 -> 8/3
        assert!((get_average_umi_distance(&umis) - 8.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn average_distance_is_infinite_with_mixed_lengths() {
        let umis = vec![b"ACGT".to_vec(), b"ACG".to_vec()];
        assert!(get_average_umi_distance(&umis).is_infinite());
    }

    #[test]
    fn random_generator_rejects_an_empty_table() {
        assert!(RandomReadGenerator::new(&[], 1).is_err());
        assert!(RandomReadGenerator::new(&[(b"AA".to_vec(), 0)], 1).is_err());
    }

    #[test]
    fn random_generator_is_reproducible_for_a_seed() {
        let freq = [(b"AA".to_vec(), 3), (b"CC".to_vec(), 1)];
        let a = RandomReadGenerator::new(&freq, 42).unwrap().get_umis(50);
        let b = RandomReadGenerator::new(&freq, 42).unwrap().get_umis(50);
        assert_eq!(a, b);
        let c = RandomReadGenerator::new(&freq, 43).unwrap().get_umis(50);
        assert_ne!(a, c);
    }

    /// The draws must match `NumPy`'s, not merely its distribution: `dedup`'s
    /// statistics files are compared byte for byte. Reference produced by
    /// `np.random.seed(123456789); np.random.choice(['a','b','c'], 8,
    /// p=[0.5, 0.3, 0.2])`.
    #[test]
    fn random_generator_matches_numpy_choice() {
        let freq = [
            (b"a".to_vec(), 50),
            (b"b".to_vec(), 30),
            (b"c".to_vec(), 20),
        ];
        let drawn = RandomReadGenerator::new(&freq, 123_456_789)
            .unwrap()
            .get_umis(8);
        let expected: Vec<Vec<u8>> = ["b", "b", "b", "b", "a", "b", "c", "a"]
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect();
        assert_eq!(drawn, expected);
    }

    #[test]
    fn random_generator_respects_the_frequency_distribution() {
        let freq = [(b"AA".to_vec(), 9), (b"CC".to_vec(), 1)];
        let draws = RandomReadGenerator::new(&freq, 7).unwrap().get_umis(10_000);
        let aa = draws.iter().filter(|d| *d == b"AA").count();
        // 90% expected; a 5-point band is far outside sampling noise at n=10k.
        assert!((8500..9500).contains(&aa), "AA drawn {aa} times of 10000");
    }

    #[test]
    fn random_generator_only_returns_known_umis() {
        let freq = [(b"AA".to_vec(), 1), (b"CC".to_vec(), 1)];
        let draws = RandomReadGenerator::new(&freq, 3).unwrap().get_umis(100);
        assert!(draws.iter().all(|d| d == b"AA" || d == b"CC"));
    }
}
