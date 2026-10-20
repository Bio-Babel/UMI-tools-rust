//! Opening, writing and sorting SAM, BAM and CRAM files
//! (Python source: `Utilities.open_input_alignments`,
//! `Utilities.open_output_alignments` and `Utilities.sort_output`).
//!
//! The Python original delegates all three to `pysam`, which wraps htslib;
//! `sort_output` in particular shells out to `pysam.sort`, the `samtools sort`
//! entry point. This module uses `noodles` instead, so the port has no C
//! dependency and no samtools on `PATH` — and it performs the coordinate sort
//! itself rather than shelling out.
//!
//! Records are streamed, never read in whole, and sorted output goes through
//! [`CoordinateSorter`], which spills to disk as `samtools sort` does.

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::Header;

use crate::error::{PortError, Result};
pub use crate::record::BamRecord;
use crate::utilities::AlignmentFormat;

/// An opened alignment file, ready to yield records.
///
/// The three formats need three different readers with three different type
/// parameters, so they are an enum rather than a `Box<dyn>` — which also keeps
/// the CRAM reader's reference-sequence repository attached to the reader that
/// needs it.
pub enum AlignmentReader {
    /// Text SAM.
    Sam(Box<noodles_sam::io::Reader<BufReader<File>>>),
    /// BGZF-compressed BAM.
    Bam(Box<noodles_bam::io::Reader<noodles_bgzf::io::Reader<File>>>),
    /// Reference-compressed CRAM.
    Cram(Box<noodles_cram::io::Reader<File>>),
}

/// An alignment file and the header read from it.
pub struct OpenedAlignment {
    /// The reader.
    pub reader: AlignmentReader,
    /// The header, needed to resolve reference-sequence ids to names.
    pub header: Header,
    header_text: Option<Vec<u8>>,
    path: PathBuf,
    first_records: FirstRecords,
    stream: StreamState,
}

/// Each reference's first record, from the BAM index, read on first use.
enum FirstRecords {
    Unread,
    Absent,
    Read(Vec<Option<noodles_bgzf::VirtualPosition>>),
}

#[derive(Default)]
struct StreamState {
    buf: RecordBuf,
    cram: Option<CramState>,
}

struct CramState {
    repository: noodles_fasta::repository::Repository,
    container: noodles_cram::io::reader::Container,
    pending: std::vec::IntoIter<BamRecord>,
}

impl OpenedAlignment {
    /// Reads every record into memory, decoded (for `prepare-for-rsem`).
    ///
    /// # Errors
    ///
    /// Returns [`PortError::Alignment`] if a record cannot be decoded.
    pub fn read_all(&mut self) -> Result<Vec<RecordBuf>> {
        let header = self.header.clone();
        let mut out = Vec::new();
        match &mut self.reader {
            AlignmentReader::Sam(reader) => {
                for record in reader.record_bufs(&header) {
                    out.push(record.map_err(|e| PortError::alignment(e.to_string()))?);
                }
            }
            AlignmentReader::Bam(reader) => {
                for record in reader.record_bufs(&header) {
                    out.push(record.map_err(|e| PortError::alignment(e.to_string()))?);
                }
            }
            AlignmentReader::Cram(reader) => {
                // The CRAM reader yields lazy `Record`s rather than
                // `RecordBuf`s, so each is converted here; the rest of the
                // crate only ever sees the owned form.
                for record in reader.records(&header) {
                    let record = record.map_err(|e| PortError::alignment(e.to_string()))?;
                    out.push(
                        RecordBuf::try_from_alignment_record(&header, &record)
                            .map_err(|e| PortError::alignment(e.to_string()))?,
                    );
                }
            }
        }
        Ok(out)
    }

    /// Reads the next record, in file order; `None` at the end of the file.
    pub fn read_record(&mut self) -> Result<Option<BamRecord>> {
        self.read_record_io()
            .map_err(|e| PortError::alignment(e.to_string()))
    }

    fn read_record_io(&mut self) -> std::io::Result<Option<BamRecord>> {
        match &mut self.reader {
            AlignmentReader::Bam(reader) => BamRecord::read_from(reader.get_mut()),
            AlignmentReader::Sam(reader) => {
                if reader.read_record_buf(&self.header, &mut self.stream.buf)? == 0 {
                    return Ok(None);
                }
                BamRecord::from_alignment_record(&self.header, &self.stream.buf).map(Some)
            }
            AlignmentReader::Cram(reader) => {
                let Some(cram) = self.stream.cram.as_mut() else {
                    return Ok(None);
                };
                loop {
                    if let Some(record) = cram.pending.next() {
                        return Ok(Some(record));
                    }
                    // As `noodles_cram::io::reader::Records`, without borrowing the reader.
                    if reader.read_container(&mut cram.container)? == 0 {
                        return Ok(None);
                    }
                    let compression_header = cram.container.compression_header()?;
                    let mut decoded = Vec::new();
                    for slice in cram.container.slices() {
                        let slice = slice?;
                        let (core, external) = slice.decode_blocks()?;
                        for record in slice.records(
                            cram.repository.clone(),
                            &self.header,
                            &compression_header,
                            &core,
                            &external,
                        )? {
                            decoded.push(BamRecord::from_alignment_record(&self.header, &record)?);
                        }
                    }
                    cram.pending = decoded.into_iter();
                }
            }
        }
    }

    /// Streams the records (Python source: `fetch()`, which skips unplaced
    /// records, or `fetch(until_eof=True)`).
    pub fn records(&mut self, until_eof: bool) -> Records<'_> {
        Records {
            opened: self,
            until_eof,
        }
    }

    /// The header text as the file spells it (`None` for CRAM), which
    /// Python's `template=infile` output carries over verbatim.
    #[must_use]
    pub fn header_text(&self) -> Option<&[u8]> {
        self.header_text.as_deref()
    }

    /// Whether the file has a BAM index (`.bai` or `.csi`).
    pub fn has_index(&mut self) -> Result<bool> {
        Ok(self.first_records()?.is_some())
    }

    fn first_records(&mut self) -> Result<Option<&[Option<noodles_bgzf::VirtualPosition>]>> {
        if matches!(self.first_records, FirstRecords::Unread) {
            let found = match self.reader {
                AlignmentReader::Bam(_) => read_first_records(&self.path)?,
                AlignmentReader::Sam(_) | AlignmentReader::Cram(_) => None,
            };
            self.first_records = found.map_or(FirstRecords::Absent, FirstRecords::Read);
        }
        Ok(match &self.first_records {
            FirstRecords::Read(starts) => Some(starts),
            FirstRecords::Unread | FirstRecords::Absent => None,
        })
    }

    /// Streams the records placed on one reference (Python source:
    /// `fetch(reference=...)`): seeks with an index, else reads forward.
    pub fn fetch_reference(&mut self, id: usize) -> Result<ReferenceRecords<'_>> {
        let scan = self.start_reference(id)?;
        Ok(ReferenceRecords { opened: self, scan })
    }

    pub(crate) fn start_reference(&mut self, id: usize) -> Result<ReferenceScan> {
        let mut scan = ReferenceScan {
            id,
            started: false,
            done: false,
        };
        if let Some(starts) = self.first_records()? {
            match starts.get(id).copied().flatten() {
                None => scan.done = true,
                Some(start) => {
                    let AlignmentReader::Bam(reader) = &mut self.reader else {
                        unreachable!("only BAM indices are read");
                    };
                    reader
                        .get_mut()
                        .seek(start)
                        .map_err(|e| PortError::alignment(e.to_string()))?;
                }
            }
        }
        Ok(scan)
    }

    /// The reference-sequence names, in header order.
    ///
    /// Equivalent to `pysam`'s `AlignmentFile.references`, which
    /// [`crate::sam::get_metacontig2contig`] needs.
    #[must_use]
    pub fn reference_names(&self) -> Vec<String> {
        self.header
            .reference_sequences()
            .keys()
            .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
            .collect()
    }

    /// The name of the reference sequence a record aligns to.
    #[must_use]
    pub fn reference_name(&self, record: &RecordBuf) -> Option<String> {
        let id = record.reference_sequence_id()?;
        self.header
            .reference_sequences()
            .get_index(id)
            .map(|(name, _)| String::from_utf8_lossy(name.as_ref()).into_owned())
    }
}

/// See [`OpenedAlignment::records`].
pub struct Records<'a> {
    opened: &'a mut OpenedAlignment,
    until_eof: bool,
}

impl Iterator for Records<'_> {
    type Item = Result<BamRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.opened.read_record().transpose()? {
                Ok(record) if !self.until_eof && record.reference_id().is_none() => {}
                result => return Some(result),
            }
        }
    }
}

/// See [`OpenedAlignment::fetch_reference`].
pub struct ReferenceRecords<'a> {
    opened: &'a mut OpenedAlignment,
    scan: ReferenceScan,
}

impl Iterator for ReferenceRecords<'_> {
    type Item = Result<BamRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.scan.next(self.opened)
    }
}

pub(crate) struct ReferenceScan {
    id: usize,
    started: bool,
    done: bool,
}

impl ReferenceScan {
    pub(crate) fn next(&mut self, opened: &mut OpenedAlignment) -> Option<Result<BamRecord>> {
        while !self.done {
            let record = match opened.read_record() {
                Ok(Some(record)) => record,
                Ok(None) => break,
                Err(e) => return Some(Err(e)),
            };
            match record.reference_id() {
                Some(id) if id == self.id => {
                    self.started = true;
                    return Some(Ok(record));
                }
                // Sorted input: once past this reference, no more of its records.
                Some(id) if id < self.id && !self.started => {}
                _ => self.done = true,
            }
        }
        self.done = true;
        None
    }
}

/// Looks for the index where htslib does: `path.bai`, `path.csi`, then the `.bai` beside it.
fn read_first_records(path: &Path) -> Result<Option<Vec<Option<noodles_bgzf::VirtualPosition>>>> {
    let with = |suffix: &str| {
        let mut p = path.as_os_str().to_owned();
        p.push(suffix);
        PathBuf::from(p)
    };
    let bad = |p: &Path, e: std::io::Error| PortError::alignment(format!("{}: {e}", p.display()));
    let bai = with(".bai");
    if bai.is_file() {
        let index = noodles_bam::bai::fs::read(&bai).map_err(|e| bad(&bai, e))?;
        return Ok(Some(first_chunk_starts(&index)));
    }
    let csi = with(".csi");
    if csi.is_file() {
        let index = noodles_csi::fs::read(&csi).map_err(|e| bad(&csi, e))?;
        return Ok(Some(first_chunk_starts(&index)));
    }
    let beside = path.with_extension("bai");
    if beside.is_file() && beside != path {
        let index = noodles_bam::bai::fs::read(&beside).map_err(|e| bad(&beside, e))?;
        return Ok(Some(first_chunk_starts(&index)));
    }
    Ok(None)
}

fn first_chunk_starts<I>(
    index: &noodles_csi::binning_index::Index<I>,
) -> Vec<Option<noodles_bgzf::VirtualPosition>>
where
    I: noodles_csi::binning_index::index::reference_sequence::Index,
{
    use noodles_csi::binning_index::index::reference_sequence::{bin::Chunk, Bin};
    index
        .reference_sequences()
        .iter()
        .map(|reference| {
            reference
                .bins()
                .values()
                .flat_map(Bin::chunks)
                .map(Chunk::start)
                .min()
        })
        .collect()
}

/// The input alignment, enough to open it again for the side passes Python
/// also reopens it for.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentSource {
    /// The file.
    pub path: PathBuf,
    /// Its format.
    pub format: AlignmentFormat,
    /// The reference, required for CRAM.
    pub reference: Option<PathBuf>,
}

impl AlignmentSource {
    /// Opens the file.
    pub fn open(&self) -> Result<OpenedAlignment> {
        open_input_alignments(&self.path, self.format, self.reference.as_deref())
    }
}

/// Opens an alignment file for reading
/// (Python source: `Utilities.open_input_alignments`).
///
/// `reference_filename` is required for CRAM, which cannot be decoded without
/// the reference the sequences were compressed against — the same requirement
/// `pysam` has.
///
/// # Errors
///
/// - [`PortError::Io`] if the file cannot be opened.
/// - [`PortError::Alignment`] if the header cannot be parsed, or if CRAM was
///   requested without a reference.
pub fn open_input_alignments(
    path: &Path,
    format: AlignmentFormat,
    reference_filename: Option<&Path>,
) -> Result<OpenedAlignment> {
    let file = File::open(path).map_err(|e| PortError::io(path, e))?;
    match format {
        AlignmentFormat::Sam => {
            let mut input = BufReader::new(file);
            let (header, text) = read_sam_header(&mut input).map_err(|e| bad_header(path, &e))?;
            Ok(OpenedAlignment {
                reader: AlignmentReader::Sam(Box::new(noodles_sam::io::Reader::new(input))),
                header,
                header_text: Some(text),
                path: path.to_path_buf(),
                first_records: FirstRecords::Unread,
                stream: StreamState::default(),
            })
        }
        AlignmentFormat::Bam => {
            let mut reader = noodles_bam::io::Reader::new(file);
            let (header, text) = read_bam_header(&mut reader).map_err(|e| bad_header(path, &e))?;
            Ok(OpenedAlignment {
                reader: AlignmentReader::Bam(Box::new(reader)),
                header,
                header_text: Some(text),
                path: path.to_path_buf(),
                first_records: FirstRecords::Unread,
                stream: StreamState::default(),
            })
        }
        AlignmentFormat::Cram => {
            let reference = reference_filename.ok_or_else(|| {
                PortError::alignment("CRAM input requires a reference sequence (--reference-fasta)")
            })?;
            let repository = fasta_repository(reference)?;
            let mut reader = noodles_cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repository.clone())
                .build_from_reader(file);
            let header = reader
                .read_header()
                .map_err(|e| PortError::alignment(format!("{}: {e}", path.display())))?;
            Ok(OpenedAlignment {
                reader: AlignmentReader::Cram(Box::new(reader)),
                header,
                header_text: None,
                path: path.to_path_buf(),
                first_records: FirstRecords::Unread,
                stream: StreamState {
                    cram: Some(CramState {
                        repository,
                        container: noodles_cram::io::reader::Container::default(),
                        pending: Vec::new().into_iter(),
                    }),
                    ..StreamState::default()
                },
            })
        }
    }
}

fn bad_header(path: &Path, e: &std::io::Error) -> PortError {
    PortError::alignment(format!("{}: {e}", path.display()))
}

fn parse_header_text(text: &[u8]) -> std::io::Result<Header> {
    let mut parser = noodles_sam::header::Parser::default();
    for line in text.split(|&b| b == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.is_empty() {
            continue;
        }
        parser
            .parse_partial(line)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    }
    Ok(parser.finish())
}

fn read_sam_header<R: std::io::BufRead>(input: &mut R) -> std::io::Result<(Header, Vec<u8>)> {
    let mut text = Vec::new();
    while input.fill_buf()?.first() == Some(&b'@') {
        input.read_until(b'\n', &mut text)?;
    }
    Ok((parse_header_text(&text)?, text))
}

/// Reads a BAM header, keeping its text; the binary reference list stands in
/// for missing `@SQ` lines, as in `noodles`.
fn read_bam_header<R: std::io::Read>(
    reader: &mut noodles_bam::io::Reader<R>,
) -> std::io::Result<(Header, Vec<u8>)> {
    use std::io::Read as _;

    let mut header_reader = reader.header_reader();
    if header_reader.read_magic_number()? != *b"BAM\x01" {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "invalid BAM magic number",
        ));
    }
    let mut text = Vec::new();
    header_reader
        .raw_sam_header_reader()?
        .read_to_end(&mut text)?;
    if let Some(nul) = text.iter().position(|&b| b == 0) {
        text.truncate(nul);
    }
    let binary = header_reader.read_reference_sequences()?;

    let mut header = parse_header_text(&text)?;
    if header.reference_sequences().is_empty() {
        *header.reference_sequences_mut() = binary;
    } else if header.reference_sequences().len() != binary.len()
        || header
            .reference_sequences()
            .iter()
            .zip(&binary)
            .any(|((a, a_map), (b, b_map))| a != b || a_map.length() != b_map.length())
    {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "SAM header and binary reference sequence dictionaries mismatch",
        ));
    }
    Ok((header, text))
}

/// Whether a `--reference-filename` value names a remote location rather than
/// a local file.
///
/// htslib — and so `pysam`, and so the Python original — will fetch a
/// reference over HTTP(S) or FTP. This port reads references from the local
/// filesystem only, so a URL is detected here and reported with an actionable
/// message instead of surfacing as a puzzling "No such file or directory".
/// Recorded as deviation D12 in `port_reports/umi_tools/05_design.md`.
#[must_use]
pub fn is_remote_reference(path: &Path) -> bool {
    let text = path.to_string_lossy();
    ["http://", "https://", "ftp://", "ftps://", "s3://", "gs://"]
        .iter()
        .any(|scheme| text.starts_with(scheme))
}

fn fasta_repository(path: &Path) -> Result<noodles_fasta::repository::Repository> {
    use noodles_fasta::repository::adapters::IndexedReader;

    if is_remote_reference(path) {
        return Err(PortError::alignment(format!(
            "cannot read the reference from {}: remote references are not supported. \
             htslib fetches these for the original; this port reads local files only. \
             Download the reference (and its .fai index) and pass the local path to \
             --reference-filename.",
            path.display()
        )));
    }

    let reader = noodles_fasta::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .map_err(|e| PortError::io(path, e))?;
    Ok(noodles_fasta::repository::Repository::new(
        IndexedReader::new(reader),
    ))
}

/// Where an alignment writer sends its bytes.
///
/// Python writes alignments to `options.stdout`, which is `sys.stdout` unless
/// `-S/--stdout` names a file — and the package's own test suite relies on
/// that, redirecting the command's standard output with the shell. So the
/// writer has to accept either, which means one concrete sink type rather
/// than a `File`.
pub enum OutputSink {
    /// A file on disk.
    File(BufWriter<File>),
    /// The process's standard output.
    Stdout(BufWriter<std::io::Stdout>),
}

impl Write for OutputSink {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            OutputSink::File(w) => w.write(buf),
            OutputSink::Stdout(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            OutputSink::File(w) => w.flush(),
            OutputSink::Stdout(w) => w.flush(),
        }
    }
}

/// An alignment file open for writing.
pub enum AlignmentWriter {
    /// Text SAM.
    Sam(Box<noodles_sam::io::Writer<OutputSink>>),
    /// BGZF-compressed BAM.
    Bam(Box<noodles_bam::io::Writer<noodles_bgzf::io::Writer<OutputSink>>>),
    /// Reference-compressed CRAM.
    Cram(Box<noodles_cram::io::Writer<OutputSink>>),
}

/// An alignment writer plus the header it was opened with.
pub struct OpenedAlignmentWriter {
    writer: AlignmentWriter,
    header: Header,
}

impl OpenedAlignmentWriter {
    /// Writes one record.
    ///
    /// # Errors
    ///
    /// Returns [`PortError::Alignment`] if the record cannot be encoded.
    pub fn write(&mut self, record: &RecordBuf) -> Result<()> {
        let result = match &mut self.writer {
            AlignmentWriter::Sam(w) => w.write_alignment_record(&self.header, record),
            AlignmentWriter::Bam(w) => w.write_alignment_record(&self.header, record),
            AlignmentWriter::Cram(w) => {
                check_cram_writable(record)?;
                w.write_alignment_record(&self.header, record)
            }
        };
        result.map_err(|e| PortError::alignment(e.to_string()))
    }

    /// Writes one [`BamRecord`]; to BAM its bytes pass through unchanged, as in htslib.
    pub fn write_bam_record(&mut self, record: &BamRecord) -> Result<()> {
        let result = match &mut self.writer {
            AlignmentWriter::Bam(w) => record.write_to(w.get_mut()),
            AlignmentWriter::Sam(w) => w.write_alignment_record(&self.header, &record.view()),
            AlignmentWriter::Cram(w) => {
                let view = record.view();
                if view.sequence().is_empty() && !view.cigar().is_empty() {
                    return Err(cram_unwritable(record.name()));
                }
                w.write_alignment_record(&self.header, &view)
            }
        };
        result.map_err(|e| PortError::alignment(e.to_string()))
    }

    /// Flushes and closes the file.
    ///
    /// Consuming `self` makes writing after close a compile error, where
    /// `pysam`'s `close()` leaves the object usable and failing at runtime.
    ///
    /// # Errors
    ///
    /// Returns [`PortError::Alignment`] if the trailer cannot be written.
    pub fn finish(self) -> Result<()> {
        match self.writer {
            AlignmentWriter::Sam(mut w) => w
                .get_mut()
                .flush()
                .map_err(|e| PortError::alignment(e.to_string())),
            AlignmentWriter::Bam(w) => w
                .into_inner()
                .finish()
                .and_then(|mut sink| sink.flush())
                .map_err(|e| PortError::alignment(e.to_string())),
            AlignmentWriter::Cram(mut w) => w
                .try_finish(&self.header)
                .map_err(|e| PortError::alignment(e.to_string())),
        }
    }
}

/// Rejects records the CRAM writer cannot encode, before it panics on them.
///
/// A SAM record may legitimately carry `*` for its sequence — the package's own
/// `chr19_1mb.cram` fixture is entirely such records — but `noodles-cram`'s
/// writer assumes a sequence is present when it converts the CIGAR into CRAM
/// features, and indexes past the end of an empty slice. That is an upstream
/// panic on ordinary user input, which this crate must not propagate: a
/// malformed-for-this-writer record is a recoverable error, not a bug in the
/// caller.
///
/// Recorded as deviation D13 in `port_reports/umi_tools/05_design.md`: CRAM
/// *output* of sequence-less records is unsupported, while CRAM input and
/// SAM/BAM output of the same records work normally.
fn check_cram_writable(record: &RecordBuf) -> Result<()> {
    if record.sequence().is_empty() && !record.cigar().as_ref().is_empty() {
        return Err(cram_unwritable(record.name().map(AsRef::as_ref)));
    }
    Ok(())
}

fn cram_unwritable(name: Option<&[u8]>) -> PortError {
    PortError::alignment(format!(
        "cannot write {} to CRAM: the record has a CIGAR but no sequence, \
         which this CRAM writer cannot encode. Write SAM or BAM instead, \
         or use --no-sort-output with a non-CRAM output format.",
        name.map_or_else(
            || "<unnamed record>".to_string(),
            |n| String::from_utf8_lossy(n).into_owned()
        )
    ))
}

/// Opens an alignment file for writing, using `header` as the template
/// (Python source: `Utilities.open_output_alignments`, whose `template=infile`
/// argument does the same job).
///
/// SAM and BAM headers are written as `header_text` when given; CRAM always from `header`.
///
/// # Errors
///
/// - [`PortError::Io`] if the file cannot be created.
/// - [`PortError::Alignment`] if the header cannot be written, or if CRAM was
///   requested without a reference.
pub fn open_output_alignments(
    path: &Path,
    format: AlignmentFormat,
    header: &Header,
    header_text: Option<&[u8]>,
    reference_filename: Option<&Path>,
) -> Result<OpenedAlignmentWriter> {
    let file = File::create(path).map_err(|e| PortError::io(path, e))?;
    open_output_alignments_to(
        OutputSink::File(BufWriter::new(file)),
        format,
        header,
        header_text,
        reference_filename,
    )
}

/// Opens an alignment writer over standard output
/// (Python source: the `options.stdout is sys.stdout` case of
/// `Utilities.open_output_alignments`).
///
/// # Errors
///
/// - [`PortError::Alignment`] if the header cannot be written, or if CRAM was
///   requested without a reference.
pub fn open_output_alignments_to_stdout(
    format: AlignmentFormat,
    header: &Header,
    header_text: Option<&[u8]>,
    reference_filename: Option<&Path>,
) -> Result<OpenedAlignmentWriter> {
    open_output_alignments_to(
        OutputSink::Stdout(BufWriter::new(std::io::stdout())),
        format,
        header,
        header_text,
        reference_filename,
    )
}

/// Opens an alignment writer over an arbitrary sink.
///
/// # Errors
///
/// - [`PortError::Alignment`] if the header cannot be written, or if CRAM was
///   requested without a reference.
pub fn open_output_alignments_to(
    inner: OutputSink,
    format: AlignmentFormat,
    header: &Header,
    header_text: Option<&[u8]>,
    reference_filename: Option<&Path>,
) -> Result<OpenedAlignmentWriter> {
    let failed = |e: std::io::Error| PortError::alignment(e.to_string());
    let writer = match format {
        AlignmentFormat::Sam => {
            let mut w = noodles_sam::io::Writer::new(inner);
            match header_text {
                Some(text) => w.get_mut().write_all(text).map_err(failed)?,
                None => w.write_header(header).map_err(failed)?,
            }
            AlignmentWriter::Sam(Box::new(w))
        }
        AlignmentFormat::Bam => {
            let mut w = noodles_bam::io::Writer::new(inner);
            match header_text {
                Some(text) => write_bam_header(w.get_mut(), text, header).map_err(failed)?,
                None => w.write_header(header).map_err(failed)?,
            }
            AlignmentWriter::Bam(Box::new(w))
        }
        AlignmentFormat::Cram => {
            let reference = reference_filename.ok_or_else(|| {
                PortError::alignment(
                    "CRAM output requires a reference sequence (--reference-fasta)",
                )
            })?;
            let repository = fasta_repository(reference)?;
            let mut w = noodles_cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repository)
                .build_from_writer(inner);
            w.write_header(header)
                .map_err(|e| PortError::alignment(e.to_string()))?;
            AlignmentWriter::Cram(Box::new(w))
        }
    };

    Ok(OpenedAlignmentWriter {
        writer,
        header: header.clone(),
    })
}

fn write_bam_header<W: Write>(writer: &mut W, text: &[u8], header: &Header) -> std::io::Result<()> {
    let to_u32 = |n: usize| {
        u32::try_from(n).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))
    };
    writer.write_all(b"BAM\x01")?;
    writer.write_all(&to_u32(text.len())?.to_le_bytes())?;
    writer.write_all(text)?;
    let references = header.reference_sequences();
    writer.write_all(&to_u32(references.len())?.to_le_bytes())?;
    for (name, map) in references {
        let name: &[u8] = name.as_ref();
        writer.write_all(&to_u32(name.len() + 1)?.to_le_bytes())?;
        writer.write_all(name)?;
        writer.write_all(&[0])?;
        writer.write_all(&to_u32(usize::from(map.length()))?.to_le_bytes())?;
    }
    Ok(())
}

/// Marks header text as coordinate-sorted as `samtools sort` does: `SO`
/// replaced or appended, `GO`/`SS` removed, `@HD VN:1.6 SO:coordinate` if absent.
#[must_use]
pub fn stamp_coordinate_sort_order_text(text: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(text.len() + 32);
    let mut stamped = false;
    for line in text.split_inclusive(|&b| b == b'\n') {
        let body = line.strip_suffix(b"\n").unwrap_or(line);
        let body = body.strip_suffix(b"\r").unwrap_or(body);
        if stamped || !(body == b"@HD" || body.starts_with(b"@HD\t")) {
            out.extend_from_slice(line);
            continue;
        }
        stamped = true;
        let mut fields: Vec<&[u8]> = body
            .split(|&b| b == b'\t')
            .filter(|f| !f.starts_with(b"GO:") && !f.starts_with(b"SS:"))
            .collect();
        match fields.iter().position(|f| f.starts_with(b"SO:")) {
            Some(i) => fields[i] = b"SO:coordinate",
            None => fields.push(b"SO:coordinate"),
        }
        out.extend_from_slice(&fields.join(&b'\t'));
        out.push(b'\n');
    }
    if !stamped {
        let mut with_hd = b"@HD\tVN:1.6\tSO:coordinate\n".to_vec();
        with_hd.extend_from_slice(&out);
        return with_hd;
    }
    out
}

/// Stamp `@HD VN:1.6 SO:coordinate` onto a header that is about to be written
/// as coordinate-sorted output.
///
/// The Python original sorts by shelling out to `pysam.sort`, i.e. `samtools
/// sort`, which rewrites the `@HD` line of whatever it emits: it sets the
/// format version and declares `SO:coordinate`. Several of the test fixtures
/// carry no `@HD` line at all, so copying the input header verbatim drops the
/// sort-order declaration that the Python output always has — a difference
/// that matters, because downstream tools read `SO` to decide whether a file
/// can be indexed or must be re-sorted.
///
/// An existing `@HD` is preserved apart from its `SO` field, matching
/// `samtools sort`, which overwrites `SO` but leaves other header tags alone.
fn stamp_coordinate_sort_order(header: &mut Header) {
    use noodles_sam::header::record::value::map::{
        self,
        header::{sort_order, tag, Version},
    };
    use noodles_sam::header::record::value::Map;

    let hd = header
        .header_mut()
        .get_or_insert_with(|| Map::<map::Header>::new(Version::new(1, 6)));
    hd.other_fields_mut()
        .insert(tag::SORT_ORDER, sort_order::COORDINATE.into());
    hd.other_fields_mut().shift_remove(&tag::GROUP_ORDER);
    hd.other_fields_mut().shift_remove(&tag::SUBSORT_ORDER);
}

/// Bytes buffered before a sorted run is spilled: a third of samtools' 768 MiB,
/// since here the buffer coexists with the bundles in flight.
pub const SORT_MEMORY: usize = 256 << 20;

/// Per-record memory beyond the encoded bytes.
const RECORD_OVERHEAD: usize = 32;

/// Sorts records in samtools' coordinate order (Python source: `pysam.sort`)
/// within a memory budget, spilling sorted runs to temporary files.
pub struct CoordinateSorter {
    buffer: Vec<BamRecord>,
    buffered: usize,
    limit: usize,
    tmpdir: Option<PathBuf>,
    runs: Vec<tempfile::NamedTempFile>,
}

impl CoordinateSorter {
    /// A sorter that spills to `tmpdir` (system default if `None`) past `limit` bytes.
    #[must_use]
    pub fn new(tmpdir: Option<&Path>, limit: usize) -> Self {
        Self {
            buffer: Vec::new(),
            buffered: 0,
            limit,
            tmpdir: tmpdir.map(Path::to_path_buf),
            runs: Vec::new(),
        }
    }

    /// Adds a record.
    pub fn push(&mut self, record: BamRecord) -> Result<()> {
        self.buffered += record.heap_size() + RECORD_OVERHEAD;
        self.buffer.push(record);
        if self.buffered >= self.limit {
            self.spill()?;
        }
        Ok(())
    }

    fn spill(&mut self) -> Result<()> {
        let dir = self.tmpdir.clone().unwrap_or_else(std::env::temp_dir);
        let io_err = |e| PortError::io(&dir, e);
        let run = tempfile::Builder::new()
            .prefix("umi_tools.sort.")
            .suffix(".bam")
            .tempfile_in(&dir)
            .map_err(io_err)?;
        let mut writer = noodles_bgzf::io::writer::Builder::default()
            .set_compression_level(noodles_bgzf::io::writer::CompressionLevel::FAST)
            .build_from_writer(BufWriter::new(run.reopen().map_err(io_err)?));

        // Stable, so equal keys keep their push order.
        self.buffer.sort_by_key(BamRecord::coordinate_key);
        for record in self.buffer.drain(..) {
            record.write_to(&mut writer).map_err(io_err)?;
        }
        writer.finish().map_err(io_err)?.flush().map_err(io_err)?;
        self.buffered = 0;
        self.runs.push(run);
        Ok(())
    }

    /// Writes every record to `writer` in coordinate order, then closes it.
    pub fn finish(mut self, mut writer: OpenedAlignmentWriter) -> Result<()> {
        self.buffer.sort_by_key(BamRecord::coordinate_key);
        if self.runs.is_empty() {
            for record in &self.buffer {
                writer.write_bam_record(record)?;
            }
            return writer.finish();
        }

        // The in-memory remainder is the newest run; ties go to the earlier run.
        let mut sources: Vec<RunReader> = Vec::with_capacity(self.runs.len() + 1);
        for run in &self.runs {
            let file = run.reopen().map_err(|e| PortError::io(run.path(), e))?;
            sources.push(RunReader::File(noodles_bgzf::io::Reader::new(file)));
        }
        sources.push(RunReader::Memory(
            std::mem::take(&mut self.buffer).into_iter(),
        ));

        let mut heads: Vec<Option<BamRecord>> = Vec::with_capacity(sources.len());
        let mut heap = BinaryHeap::with_capacity(sources.len());
        for (i, source) in sources.iter_mut().enumerate() {
            let head = source.next()?;
            if let Some(record) = &head {
                heap.push(Reverse((record.coordinate_key(), i)));
            }
            heads.push(head);
        }
        while let Some(Reverse((_, i))) = heap.pop() {
            if let Some(record) = heads[i].take() {
                writer.write_bam_record(&record)?;
            }
            heads[i] = sources[i].next()?;
            if let Some(record) = &heads[i] {
                heap.push(Reverse((record.coordinate_key(), i)));
            }
        }
        writer.finish()
    }
}

enum RunReader {
    File(noodles_bgzf::io::Reader<File>),
    Memory(std::vec::IntoIter<BamRecord>),
}

impl RunReader {
    fn next(&mut self) -> Result<Option<BamRecord>> {
        match self {
            RunReader::File(reader) => BamRecord::read_from(reader)
                .map_err(|e| PortError::alignment(format!("reading a sort run: {e}"))),
            RunReader::Memory(records) => Ok(records.next()),
        }
    }
}

/// Where `dedup` and `group` write: sorted unless `--no-sort-output`.
pub enum AlignmentSink {
    /// `--no-sort-output`: straight to the output, in the order written.
    Unsorted(OpenedAlignmentWriter),
    /// Coordinate-sorted on [`finish`](Self::finish).
    Sorted {
        /// The records so far.
        sorter: CoordinateSorter,
        /// The eventual output.
        writer: OpenedAlignmentWriter,
    },
}

impl AlignmentSink {
    /// Opens `path` (standard output if `None`); sorted output declares `SO:coordinate`.
    pub fn open(
        path: Option<&Path>,
        format: AlignmentFormat,
        header: &Header,
        header_text: Option<&[u8]>,
        reference_filename: Option<&Path>,
        sort: bool,
        tmpdir: Option<&Path>,
    ) -> Result<Self> {
        let mut header = header.clone();
        let mut text = header_text.map(<[u8]>::to_vec);
        if sort {
            stamp_coordinate_sort_order(&mut header);
            text = text.as_deref().map(stamp_coordinate_sort_order_text);
        }
        let text = text.as_deref();
        let writer = match path {
            Some(path) => open_output_alignments(path, format, &header, text, reference_filename)?,
            None => open_output_alignments_to_stdout(format, &header, text, reference_filename)?,
        };
        Ok(if sort {
            AlignmentSink::Sorted {
                sorter: CoordinateSorter::new(tmpdir, SORT_MEMORY),
                writer,
            }
        } else {
            AlignmentSink::Unsorted(writer)
        })
    }

    /// Writes one record.
    pub fn write(&mut self, record: BamRecord) -> Result<()> {
        match self {
            AlignmentSink::Unsorted(writer) => writer.write_bam_record(&record),
            AlignmentSink::Sorted { sorter, .. } => sorter.push(record),
        }
    }

    /// Sorts if required, and closes the output.
    pub fn finish(self) -> Result<()> {
        match self {
            AlignmentSink::Unsorted(writer) => writer.finish(),
            AlignmentSink::Sorted { sorter, writer } => sorter.finish(writer),
        }
    }
}

/// Sorts an alignment file by coordinate and deletes the input
/// (Python source: `Utilities.sort_output`).
///
/// Records are ordered as `samtools sort` orders them; see [`CoordinateSorter`].
///
/// The input file is removed once the sorted output is written, matching
/// Python's documented "Deletes the input file" behaviour.
///
/// # Errors
///
/// - [`PortError::Io`] if either file cannot be opened, or the input cannot be
///   removed.
/// - [`PortError::Alignment`] if a record cannot be read or written.
pub fn sort_output(
    input: &Path,
    sorted_out: &Path,
    input_format: AlignmentFormat,
    output_format: AlignmentFormat,
    reference_filename: Option<&Path>,
) -> Result<()> {
    sort_output_to(
        input,
        Some(sorted_out),
        input_format,
        output_format,
        reference_filename,
    )
}

/// Sorts an alignment file by coordinate, writing to a file or to standard
/// output, and deletes the input
/// (Python source: `Utilities.sort_output`, whose `sorted_out_name` is `"-"`
/// when the command's output goes to a pipe).
///
/// # Errors
///
/// - [`PortError::Io`] if either file cannot be opened, or the input cannot be
///   removed.
/// - [`PortError::Alignment`] if a record cannot be read or written.
pub fn sort_output_to(
    input: &Path,
    sorted_out: Option<&Path>,
    input_format: AlignmentFormat,
    output_format: AlignmentFormat,
    reference_filename: Option<&Path>,
) -> Result<()> {
    let mut opened = open_input_alignments(input, input_format, reference_filename)?;
    let header = opened.header.clone();
    let mut sink = AlignmentSink::open(
        sorted_out,
        output_format,
        &header,
        opened.header_text(),
        reference_filename,
        true,
        None,
    )?;
    for record in opened.records(true) {
        sink.write(record?)?;
    }
    drop(opened);
    sink.finish()?;

    std::fs::remove_file(input).map_err(|e| PortError::io(input, e))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_core::Position;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    use std::num::NonZeroUsize;

    fn header() -> Header {
        let mut builder = Header::builder();
        for name in ["chr1", "chr2"] {
            builder = builder.add_reference_sequence(
                name.as_bytes(),
                Map::<ReferenceSequence>::new(NonZeroUsize::new(1_000_000).unwrap()),
            );
        }
        builder.build()
    }

    fn record(name: &str, ref_id: Option<usize>, start: Option<usize>) -> RecordBuf {
        let mut r = RecordBuf::default();
        *r.flags_mut() = if ref_id.is_some() {
            Flags::empty()
        } else {
            Flags::UNMAPPED
        };
        *r.name_mut() = Some(name.into());
        *r.reference_sequence_id_mut() = ref_id;
        *r.alignment_start_mut() = start.and_then(Position::new);
        r
    }

    fn write_sam(path: &Path, records: &[RecordBuf]) {
        let mut w =
            open_output_alignments(path, AlignmentFormat::Sam, &header(), None, None).unwrap();
        for r in records {
            w.write(r).unwrap();
        }
        w.finish().unwrap();
    }

    #[test]
    fn sam_round_trips() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.sam");
        let records = vec![
            record("r1", Some(0), Some(100)),
            record("r2", Some(1), Some(50)),
        ];
        write_sam(&path, &records);

        let mut opened = open_input_alignments(&path, AlignmentFormat::Sam, None).unwrap();
        let read_back = opened.read_all().unwrap();
        assert_eq!(read_back.len(), 2);
        assert_eq!(read_back[0].name().unwrap().to_vec(), b"r1".to_vec());
        assert_eq!(opened.reference_names(), vec!["chr1", "chr2"]);
        assert_eq!(
            opened.reference_name(&read_back[1]).as_deref(),
            Some("chr2")
        );
    }

    #[test]
    fn bam_round_trips() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.bam");
        let records = vec![record("r1", Some(0), Some(100))];
        let mut w =
            open_output_alignments(&path, AlignmentFormat::Bam, &header(), None, None).unwrap();
        for r in &records {
            w.write(r).unwrap();
        }
        w.finish().unwrap();

        let mut opened = open_input_alignments(&path, AlignmentFormat::Bam, None).unwrap();
        let read_back = opened.read_all().unwrap();
        assert_eq!(read_back.len(), 1);
        assert_eq!(read_back[0].alignment_start(), Position::new(100));
    }

    #[test]
    fn sorting_orders_by_reference_then_position_and_removes_the_input() {
        let dir = tempfile::tempdir().unwrap();
        let unsorted = dir.path().join("unsorted.sam");
        let sorted = dir.path().join("sorted.sam");

        write_sam(
            &unsorted,
            &[
                record("c", Some(1), Some(10)),
                record("a", Some(0), Some(500)),
                record("b", Some(0), Some(100)),
            ],
        );
        sort_output(
            &unsorted,
            &sorted,
            AlignmentFormat::Sam,
            AlignmentFormat::Sam,
            None,
        )
        .unwrap();

        assert!(!unsorted.exists(), "the input should be deleted");
        let mut opened = open_input_alignments(&sorted, AlignmentFormat::Sam, None).unwrap();
        let names: Vec<Vec<u8>> = opened
            .read_all()
            .unwrap()
            .iter()
            .map(|r| r.name().unwrap().to_vec())
            .collect();
        assert_eq!(names, vec![b"b".to_vec(), b"a".to_vec(), b"c".to_vec()]);
    }

    #[test]
    fn unmapped_records_sort_last() {
        let dir = tempfile::tempdir().unwrap();
        let unsorted = dir.path().join("u.sam");
        let sorted = dir.path().join("s.sam");
        write_sam(
            &unsorted,
            &[
                record("unmapped", None, None),
                record("mapped", Some(0), Some(10)),
            ],
        );
        sort_output(
            &unsorted,
            &sorted,
            AlignmentFormat::Sam,
            AlignmentFormat::Sam,
            None,
        )
        .unwrap();

        let mut opened = open_input_alignments(&sorted, AlignmentFormat::Sam, None).unwrap();
        let names: Vec<Vec<u8>> = opened
            .read_all()
            .unwrap()
            .iter()
            .map(|r| r.name().unwrap().to_vec())
            .collect();
        assert_eq!(names, vec![b"mapped".to_vec(), b"unmapped".to_vec()]);
    }

    /// `samtools sort` — which the Python original shells out to via
    /// `pysam.sort` — stamps `@HD VN:1.6 SO:coordinate` on its output. Test
    /// fixtures frequently carry no `@HD` line at all, so copying the input
    /// header verbatim silently drops the sort-order declaration that every
    /// Python-produced file has.
    #[test]
    fn sorted_output_declares_coordinate_sort_order() {
        let dir = tempfile::tempdir().unwrap();
        let unsorted = dir.path().join("u.sam");
        let sorted = dir.path().join("s.sam");

        // header() builds a header with no @HD line, like the BAM/CRAM fixtures.
        assert!(header().header().is_none());
        write_sam(&unsorted, &[record("r1", Some(0), Some(10))]);

        sort_output(
            &unsorted,
            &sorted,
            AlignmentFormat::Sam,
            AlignmentFormat::Sam,
            None,
        )
        .unwrap();

        let text = std::fs::read_to_string(&sorted).unwrap();
        assert_eq!(
            text.lines().next(),
            Some("@HD\tVN:1.6\tSO:coordinate"),
            "sorted output must declare its sort order; got:\n{text}"
        );
    }

    fn names_in(path: &Path) -> Vec<String> {
        let mut opened = open_input_alignments(path, AlignmentFormat::Sam, None).unwrap();
        opened
            .records(true)
            .map(|r| String::from_utf8_lossy(r.unwrap().name().unwrap()).into_owned())
            .collect()
    }

    #[test]
    fn spilled_runs_merge_stably_in_samtools_order() {
        let dir = tempfile::tempdir().unwrap();
        let spill = tempfile::tempdir().unwrap();
        let out = dir.path().join("sorted.sam");

        let mut input = Vec::new();
        for (i, (tid, start, reverse)) in [
            (Some(1), Some(10), false),
            (Some(0), Some(50), true),
            (None, None, false),
            (Some(0), Some(50), false),
            (Some(0), Some(50), true),
            (Some(0), Some(7), false),
            (Some(1), Some(10), false),
            (Some(0), Some(50), false),
            (None, None, false),
        ]
        .into_iter()
        .enumerate()
        {
            let mut r = record(&format!("r{i}"), tid, start);
            if reverse {
                *r.flags_mut() = Flags::REVERSE_COMPLEMENTED;
            }
            input.push(BamRecord::from_alignment_record(&header(), &r).unwrap());
        }
        let mut expected: Vec<(u64, String)> = input
            .iter()
            .map(|r| {
                (
                    r.coordinate_key(),
                    String::from_utf8_lossy(r.name().unwrap()).into_owned(),
                )
            })
            .collect();
        expected.sort_by_key(|(key, _)| *key);

        // A budget of one byte spills after every record.
        let mut sorter = CoordinateSorter::new(Some(spill.path()), 1);
        for r in input {
            sorter.push(r).unwrap();
        }
        assert!(sorter.runs.len() > 1, "the test must exercise the merge");
        let writer =
            open_output_alignments(&out, AlignmentFormat::Sam, &header(), None, None).unwrap();
        sorter.finish(writer).unwrap();

        let expected: Vec<String> = expected.into_iter().map(|(_, n)| n).collect();
        assert_eq!(
            expected,
            ["r5", "r3", "r7", "r1", "r4", "r0", "r6", "r2", "r8"]
        );
        assert_eq!(names_in(&out), expected);
        assert_eq!(std::fs::read_dir(spill.path()).unwrap().count(), 0);
    }

    #[test]
    fn an_unsorted_sink_keeps_the_order_written() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.sam");
        let mut sink = AlignmentSink::open(
            Some(&out),
            AlignmentFormat::Sam,
            &header(),
            None,
            None,
            false,
            None,
        )
        .unwrap();
        for (name, start) in [("b", 500), ("a", 100)] {
            let r = record(name, Some(0), Some(start));
            sink.write(BamRecord::from_alignment_record(&header(), &r).unwrap())
                .unwrap();
        }
        sink.finish().unwrap();
        assert_eq!(names_in(&out), vec!["b", "a"]);
        let text = std::fs::read_to_string(&out).unwrap();
        assert!(!text.contains("SO:coordinate"));
    }

    #[test]
    fn records_stream_placed_reads_unless_until_eof() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.sam");
        write_sam(
            &path,
            &[
                record("placed", Some(0), Some(10)),
                record("unplaced", None, None),
            ],
        );
        let mut opened = open_input_alignments(&path, AlignmentFormat::Sam, None).unwrap();
        assert_eq!(opened.records(false).count(), 1);
        let mut opened = open_input_alignments(&path, AlignmentFormat::Sam, None).unwrap();
        assert_eq!(opened.records(true).count(), 2);
    }

    #[test]
    fn header_text_passes_through_to_sam_and_bam() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.sam");
        let text = "@HD\tVN:1.4\n@SQ\tSN:chr1\tLN:1000\tAS:x\n@SQ\tSN:chr2\tLN:1000\n\
                    @RG\tSM:s\tID:rg1\n@PG\tPN:bwa\tID:bwa\tVN:1\n@CO\tfree\ttext\n";
        std::fs::write(
            &input,
            format!("{text}r1\t0\tchr1\t10\t30\t4M\t*\t0\t0\tACGT\tIIII\n"),
        )
        .unwrap();
        let opened = open_input_alignments(&input, AlignmentFormat::Sam, None).unwrap();
        assert_eq!(opened.header_text(), Some(text.as_bytes()));

        for format in [AlignmentFormat::Sam, AlignmentFormat::Bam] {
            let out = dir.path().join(format!("out.{}", format.as_str()));
            open_output_alignments(&out, format, &opened.header, opened.header_text(), None)
                .unwrap()
                .finish()
                .unwrap();
            let back = open_input_alignments(&out, format, None).unwrap();
            assert_eq!(back.header_text(), Some(text.as_bytes()), "{format:?}");
            assert_eq!(back.reference_names(), vec!["chr1", "chr2"]);
        }
    }

    #[test]
    fn sorted_header_text_is_stamped_as_samtools_stamps_it() {
        let stamp =
            |t: &str| String::from_utf8(stamp_coordinate_sort_order_text(t.as_bytes())).unwrap();
        assert_eq!(
            stamp("@HD\tVN:1.4\tSO:unsorted\tGO:query\tSS:x\n@SQ\tSN:c\tLN:1\n"),
            "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:c\tLN:1\n"
        );
        assert_eq!(
            stamp("@HD\tVN:1.5\n@CO\tx\n"),
            "@HD\tVN:1.5\tSO:coordinate\n@CO\tx\n"
        );
        assert_eq!(
            stamp("@SQ\tSN:c\tLN:1\n"),
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:1\n"
        );
    }

    fn indexed_bam(dir: &Path, indexed: bool) -> PathBuf {
        let path = dir.join("x.bam");
        // Sorted, so the header says so, which indexing requires.
        let mut sink = AlignmentSink::open(
            Some(&path),
            AlignmentFormat::Bam,
            &header(),
            None,
            None,
            true,
            None,
        )
        .unwrap();
        for (name, tid, start) in [
            ("a", Some(0), Some(10)),
            ("b", Some(0), Some(20)),
            ("c", Some(1), Some(5)),
            ("d", Some(1), Some(9)),
            ("u", None, None),
        ] {
            let r = record(name, tid, start);
            sink.write(BamRecord::from_alignment_record(&header(), &r).unwrap())
                .unwrap();
        }
        sink.finish().unwrap();
        if indexed {
            let index = noodles_bam::fs::index(&path).unwrap();
            noodles_bam::bai::fs::write(dir.join("x.bam.bai"), &index).unwrap();
        }
        path
    }

    #[test]
    fn a_reference_is_fetched_with_or_without_an_index() {
        for indexed in [true, false] {
            let dir = tempfile::tempdir().unwrap();
            let path = indexed_bam(dir.path(), indexed);
            let names = |id: usize| -> Vec<String> {
                let mut opened = open_input_alignments(&path, AlignmentFormat::Bam, None).unwrap();
                assert_eq!(opened.has_index().unwrap(), indexed);
                opened
                    .fetch_reference(id)
                    .unwrap()
                    .map(|r| String::from_utf8_lossy(r.unwrap().name().unwrap()).into_owned())
                    .collect()
            };
            assert_eq!(names(1), vec!["c", "d"], "indexed: {indexed}");
            assert_eq!(names(0), vec!["a", "b"], "indexed: {indexed}");
        }
    }

    #[test]
    fn sam_to_bam_conversion_preserves_records() {
        let dir = tempfile::tempdir().unwrap();
        let unsorted = dir.path().join("u.sam");
        let sorted = dir.path().join("s.bam");
        write_sam(&unsorted, &[record("r1", Some(0), Some(10))]);
        sort_output(
            &unsorted,
            &sorted,
            AlignmentFormat::Sam,
            AlignmentFormat::Bam,
            None,
        )
        .unwrap();

        let mut opened = open_input_alignments(&sorted, AlignmentFormat::Bam, None).unwrap();
        assert_eq!(opened.read_all().unwrap().len(), 1);
    }

    /// A record with a CIGAR but no sequence must be refused, not crash the
    /// process — `noodles-cram`'s writer panics on it.
    #[test]
    fn cram_rejects_sequence_less_records_instead_of_panicking() {
        use noodles_sam::alignment::record::cigar::op::{Kind, Op};
        use noodles_sam::alignment::record_buf::Cigar as CigarBuf;

        let mut r = record("r1", Some(0), Some(10));
        *r.cigar_mut() = [Op::new(Kind::Match, 33)].into_iter().collect::<CigarBuf>();
        assert!(r.sequence().is_empty());
        assert!(matches!(
            check_cram_writable(&r),
            Err(PortError::Alignment(_))
        ));

        // The same record is fine for SAM and BAM.
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.sam");
        let mut w =
            open_output_alignments(&path, AlignmentFormat::Sam, &header(), None, None).unwrap();
        assert!(w.write(&r).is_ok());
        w.finish().unwrap();
    }

    #[test]
    fn cram_without_a_reference_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.cram");
        assert!(matches!(
            open_output_alignments(&path, AlignmentFormat::Cram, &header(), None, None),
            Err(PortError::Alignment(_))
        ));
    }

    #[test]
    fn remote_references_are_reported_clearly() {
        for url in [
            "https://example.com/ref.fa",
            "http://example.com/ref.fa",
            "ftp://example.com/ref.fa",
            "s3://bucket/ref.fa",
        ] {
            assert!(
                is_remote_reference(Path::new(url)),
                "{url} should be remote"
            );
        }
        assert!(!is_remote_reference(Path::new("/local/ref.fa")));
        assert!(!is_remote_reference(Path::new("ref.fa")));

        let err = fasta_repository(Path::new("https://example.com/ref.fa")).unwrap_err();
        assert!(matches!(err, PortError::Alignment(ref m)
            if m.contains("remote references are not supported")));
    }

    #[test]
    fn opening_a_missing_file_is_an_io_error() {
        assert!(matches!(
            open_input_alignments(Path::new("/nonexistent.bam"), AlignmentFormat::Bam, None),
            Err(PortError::Io { .. })
        ));
    }
}
