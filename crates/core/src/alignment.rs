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
//! # The sort is in-memory
//!
//! `samtools sort` spills to temporary files for inputs larger than memory;
//! [`sort_output`] holds the records in memory. For the deduplicated output of
//! a UMI-tools run that is a much smaller file than the input, but it is a
//! real difference from the original for a very large output, and it is
//! recorded in `port_reports/umi_tools/08_validation.md`.

use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::Header;

use crate::error::{PortError, Result};
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
}

impl OpenedAlignment {
    /// Reads every record into memory.
    ///
    /// Python iterates lazily via `pysam`'s `fetch()`. Reading eagerly here
    /// keeps the reader enum from leaking into every consumer's type; the
    /// callers that stream (`dedup`, `group`, `count`) hand records straight
    /// to [`crate::bundles::GetBundles`], which bounds its own memory by
    /// flushing as it goes.
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
            let mut reader = noodles_sam::io::Reader::new(BufReader::new(file));
            let header = reader
                .read_header()
                .map_err(|e| PortError::alignment(format!("{}: {e}", path.display())))?;
            Ok(OpenedAlignment {
                reader: AlignmentReader::Sam(Box::new(reader)),
                header,
            })
        }
        AlignmentFormat::Bam => {
            let mut reader = noodles_bam::io::Reader::new(file);
            let header = reader
                .read_header()
                .map_err(|e| PortError::alignment(format!("{}: {e}", path.display())))?;
            Ok(OpenedAlignment {
                reader: AlignmentReader::Bam(Box::new(reader)),
                header,
            })
        }
        AlignmentFormat::Cram => {
            let reference = reference_filename.ok_or_else(|| {
                PortError::alignment("CRAM input requires a reference sequence (--reference-fasta)")
            })?;
            let repository = fasta_repository(reference)?;
            let mut reader = noodles_cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repository)
                .build_from_reader(file);
            let header = reader
                .read_header()
                .map_err(|e| PortError::alignment(format!("{}: {e}", path.display())))?;
            Ok(OpenedAlignment {
                reader: AlignmentReader::Cram(Box::new(reader)),
                header,
            })
        }
    }
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
        return Err(PortError::alignment(format!(
            "cannot write {} to CRAM: the record has a CIGAR but no sequence, \
             which this CRAM writer cannot encode. Write SAM or BAM instead, \
             or use --no-sort-output with a non-CRAM output format.",
            record.name().map_or_else(
                || "<unnamed record>".to_string(),
                |n| String::from_utf8_lossy(n.as_ref()).into_owned()
            )
        )));
    }
    Ok(())
}

/// Opens an alignment file for writing, using `header` as the template
/// (Python source: `Utilities.open_output_alignments`, whose `template=infile`
/// argument does the same job).
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
    reference_filename: Option<&Path>,
) -> Result<OpenedAlignmentWriter> {
    let file = File::create(path).map_err(|e| PortError::io(path, e))?;
    open_output_alignments_to(
        OutputSink::File(BufWriter::new(file)),
        format,
        header,
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
    reference_filename: Option<&Path>,
) -> Result<OpenedAlignmentWriter> {
    open_output_alignments_to(
        OutputSink::Stdout(BufWriter::new(std::io::stdout())),
        format,
        header,
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
    reference_filename: Option<&Path>,
) -> Result<OpenedAlignmentWriter> {
    let writer = match format {
        AlignmentFormat::Sam => {
            let mut w = noodles_sam::io::Writer::new(inner);
            w.write_header(header)
                .map_err(|e| PortError::alignment(e.to_string()))?;
            AlignmentWriter::Sam(Box::new(w))
        }
        AlignmentFormat::Bam => {
            let mut w = noodles_bam::io::Writer::new(inner);
            w.write_header(header)
                .map_err(|e| PortError::alignment(e.to_string()))?;
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

/// Filters records the way `pysam`'s `AlignmentFile.fetch()` does
/// (Python source: the `infile.fetch()` / `infile.fetch(until_eof=...)` calls
/// in `dedup.main`, `group.main` and `count.main`).
///
/// `fetch()` on an indexed file walks the index, so it yields only *placed*
/// records — the unplaced block at the end of a coordinate-sorted BAM is
/// skipped entirely. `fetch(until_eof=True)` reads straight through and yields
/// everything. `group --unmapped-reads=output` uses the second form; the other
/// commands use the first.
///
/// Reading the whole file and filtering here is equivalent, and keeps the
/// distinction explicit rather than buried in a reader flag.
#[must_use]
pub fn fetch(records: &[RecordBuf], until_eof: bool) -> Vec<&RecordBuf> {
    records
        .iter()
        .filter(|r| until_eof || r.reference_sequence_id().is_some())
        .collect()
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
}

/// Sorts an alignment file by coordinate and deletes the input
/// (Python source: `Utilities.sort_output`).
///
/// Records are ordered by reference-sequence id then alignment start, which is
/// `samtools sort`'s default and what `pysam.sort` performs. Unmapped records
/// sort last, as they do in a coordinate-sorted BAM.
///
/// The input file is removed once the sorted output is written, matching
/// Python's documented "Deletes the input file" behaviour — the input is the
/// temporary file [`crate::utilities::output_names_and_formats`] created.
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
/// The stdout case is not a corner: `dedup` and `group` sort their output
/// **even when it is a pipe**, because `output_names_and_formats` routes the
/// unsorted records through a temporary file whenever `--no-sort-output` is
/// absent. Skipping the sort there produces the right reads in the wrong
/// order, which the package's own md5-comparing test suite rejects.
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
    let mut records = opened.read_all()?;

    records.sort_by(|a, b| {
        let key = |r: &RecordBuf| {
            (
                // Unmapped records have no reference id and sort last.
                r.reference_sequence_id().unwrap_or(usize::MAX),
                r.alignment_start().map_or(usize::MAX, usize::from),
            )
        };
        key(a).cmp(&key(b))
    });

    let mut header = opened.header.clone();
    drop(opened);
    stamp_coordinate_sort_order(&mut header);

    let mut writer = match sorted_out {
        Some(path) => open_output_alignments(path, output_format, &header, reference_filename)?,
        None => open_output_alignments_to_stdout(output_format, &header, reference_filename)?,
    };
    for record in &records {
        writer.write(record)?;
    }
    writer.finish()?;

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
        let mut w = open_output_alignments(path, AlignmentFormat::Sam, &header(), None).unwrap();
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
        let mut w = open_output_alignments(&path, AlignmentFormat::Bam, &header(), None).unwrap();
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
        let mut w = open_output_alignments(&path, AlignmentFormat::Sam, &header(), None).unwrap();
        assert!(w.write(&r).is_ok());
        w.finish().unwrap();
    }

    #[test]
    fn cram_without_a_reference_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.cram");
        assert!(matches!(
            open_output_alignments(&path, AlignmentFormat::Cram, &header(), None),
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
