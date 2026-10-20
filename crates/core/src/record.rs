//! One alignment record, held as its BAM encoding (Python source:
//! `pysam.AlignedSegment`), a fraction of a decoded `RecordBuf`'s size.

use std::fmt;
use std::io::{self, Read, Write};

use noodles_bam::RecordRef;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::Header;

/// Offsets into the fixed-length part of a BAM record (SAM specification §4.2).
const REFERENCE_ID: usize = 0;
const POSITION: usize = 4;
const NAME_LENGTH: usize = 8;
const MAPPING_QUALITY: usize = 9;
const CIGAR_OP_COUNT: usize = 12;
const FLAGS: usize = 14;
const READ_LENGTH: usize = 16;
const MATE_REFERENCE_ID: usize = 20;
const MATE_POSITION: usize = 24;
const TEMPLATE_LENGTH: usize = 28;
const FIXED_LENGTH: usize = 32;

/// An alignment record in its BAM encoding: the bytes after `block_size`.
#[derive(Clone, PartialEq, Eq)]
pub struct BamRecord(Box<[u8]>);

/// A value for [`BamRecord::set_tag`]: the two kinds `group` writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TagValue<'a> {
    /// An integer, stored in the smallest BAM type that holds it.
    Int(i64),
    /// A string (`Z`).
    String(&'a [u8]),
}

impl BamRecord {
    /// Wraps an encoded record, checking (as `noodles-bam` does) that its
    /// variable-length fields fit.
    pub fn from_bytes(buf: Vec<u8>) -> io::Result<Self> {
        if buf.len() < FIXED_LENGTH {
            return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
        }
        let name_len = usize::from(buf[NAME_LENGTH]);
        let cigar_ops = usize::from(u16_at(&buf, CIGAR_OP_COUNT));
        let bases = usize::try_from(u32_at(&buf, READ_LENGTH))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let end = FIXED_LENGTH + name_len + 4 * cigar_ops + bases.div_ceil(2) + bases;
        if buf.len() < end {
            return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
        }
        Ok(Self(buf.into_boxed_slice()))
    }

    /// Encodes any alignment record the way a BAM writer would.
    pub fn from_alignment_record(
        header: &Header,
        record: &dyn noodles_sam::alignment::Record,
    ) -> io::Result<Self> {
        use noodles_sam::alignment::io::Write as _;

        // Uncompressed `block_size` + record; only the record is kept.
        let mut writer = noodles_bam::io::Writer::from(Vec::new());
        writer.write_alignment_record(header, record)?;
        let mut buf = writer.into_inner();
        buf.drain(..4);
        Self::from_bytes(buf)
    }

    /// Reads one length-prefixed record, or `None` at a clean end of input.
    pub fn read_from<R: Read>(reader: &mut R) -> io::Result<Option<Self>> {
        let mut size = [0u8; 4];
        let mut filled = 0;
        while filled < size.len() {
            match reader.read(&mut size[filled..]) {
                Ok(0) if filled == 0 => return Ok(None),
                Ok(0) => return Err(io::Error::from(io::ErrorKind::UnexpectedEof)),
                Ok(n) => filled += n,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }
        let size = usize::try_from(u32::from_le_bytes(size))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let mut buf = vec![0u8; size];
        reader.read_exact(&mut buf)?;
        Self::from_bytes(buf).map(Some)
    }

    /// Writes the record with its length prefix.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        let size = u32::try_from(self.0.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        writer.write_all(&size.to_le_bytes())?;
        writer.write_all(&self.0)
    }

    /// A field-by-field view.
    #[must_use]
    pub fn view(&self) -> RecordRef<'_> {
        // `from_bytes` guaranteed the fixed-length part is present.
        RecordRef::new(&self.0).expect("validated BAM record")
    }

    /// Decodes the record into an owned, mutable form.
    pub fn to_record_buf(&self, header: &Header) -> io::Result<RecordBuf> {
        RecordBuf::try_from_alignment_record(header, &self.view())
    }

    /// The encoded bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        &self.0
    }

    /// The bytes this record occupies on the heap, for memory budgets.
    #[must_use]
    pub fn heap_size(&self) -> usize {
        self.0.len()
    }

    /// The SAM flags.
    #[must_use]
    pub fn flags(&self) -> Flags {
        Flags::from(u16_at(&self.0, FLAGS))
    }

    /// The reference-sequence index, `None` when unplaced
    /// (`pysam`: `reference_id`, which is -1 there).
    #[must_use]
    pub fn reference_id(&self) -> Option<usize> {
        usize::try_from(i32_at(&self.0, REFERENCE_ID)).ok()
    }

    /// The 0-based leftmost coordinate, `None` when unset
    /// (`pysam`: `reference_start`).
    #[must_use]
    pub fn position(&self) -> Option<i64> {
        let pos = i32_at(&self.0, POSITION);
        (pos >= 0).then_some(i64::from(pos))
    }

    /// The mate's reference-sequence index (`pysam`: `next_reference_id`).
    #[must_use]
    pub fn mate_reference_id(&self) -> Option<usize> {
        usize::try_from(i32_at(&self.0, MATE_REFERENCE_ID)).ok()
    }

    /// The mate's 0-based leftmost coordinate (`pysam`: `next_reference_start`).
    #[must_use]
    pub fn mate_position(&self) -> Option<i64> {
        let pos = i32_at(&self.0, MATE_POSITION);
        (pos >= 0).then_some(i64::from(pos))
    }

    /// The observed template length (`pysam`: `template_length`).
    #[must_use]
    pub fn template_length(&self) -> i32 {
        i32_at(&self.0, TEMPLATE_LENGTH)
    }

    /// The raw mapping-quality byte, 255 included (`pysam`: `mapping_quality`).
    #[must_use]
    pub fn mapping_quality(&self) -> u8 {
        self.0[MAPPING_QUALITY]
    }

    /// The read name, `None` for `*` (`pysam`: `query_name`).
    #[must_use]
    pub fn name(&self) -> Option<&[u8]> {
        self.view().name().map(AsRef::as_ref)
    }

    /// The number of bases (`pysam`: `query_length`).
    #[must_use]
    pub fn sequence_length(&self) -> usize {
        u32_at(&self.0, READ_LENGTH) as usize
    }

    /// The `samtools sort` key: reference (unplaced last), `pos + 1`, strand.
    #[must_use]
    pub fn coordinate_key(&self) -> u64 {
        // tid -1 sorts after every real reference, as in samtools.
        let tid = u64::try_from(i32_at(&self.0, REFERENCE_ID)).unwrap_or(i32::MAX as u64);
        // 31 bits of tid, 32 of pos + 1 (0..=2^31) and 1 of strand.
        let pos = (i64::from(i32_at(&self.0, POSITION)) + 1) as u64;
        let reverse = u64::from(self.flags().is_reverse_complemented());
        (tid << 33) | (pos << 1) | reverse
    }

    /// Sets a tag as `pysam`'s `set_tag` does: any old one removed, the new
    /// one appended, integers in the smallest type that holds them.
    pub fn set_tag(&mut self, tag: [u8; 2], value: TagValue<'_>) -> io::Result<()> {
        let data_start = self.data_start();
        let mut buf = Vec::with_capacity(self.0.len() + 16);
        buf.extend_from_slice(&self.0[..data_start]);

        let mut src = &self.0[data_start..];
        while !src.is_empty() {
            let len = field_length(src)?;
            if src[..2] != tag {
                buf.extend_from_slice(&src[..len]);
            }
            src = &src[len..];
        }

        buf.extend_from_slice(&tag);
        match value {
            TagValue::String(s) => {
                buf.push(b'Z');
                buf.extend_from_slice(s);
                buf.push(0);
            }
            TagValue::Int(n) => {
                let out_of_range =
                    || io::Error::new(io::ErrorKind::InvalidInput, "tag value out of range");
                if n >= 0 {
                    if let Ok(v) = u8::try_from(n) {
                        buf.push(b'C');
                        buf.push(v);
                    } else if let Ok(v) = u16::try_from(n) {
                        buf.push(b'S');
                        buf.extend_from_slice(&v.to_le_bytes());
                    } else {
                        let v = u32::try_from(n).map_err(|_| out_of_range())?;
                        buf.push(b'I');
                        buf.extend_from_slice(&v.to_le_bytes());
                    }
                } else if let Ok(v) = i8::try_from(n) {
                    buf.push(b'c');
                    buf.extend_from_slice(&v.to_le_bytes());
                } else if let Ok(v) = i16::try_from(n) {
                    buf.push(b's');
                    buf.extend_from_slice(&v.to_le_bytes());
                } else {
                    let v = i32::try_from(n).map_err(|_| out_of_range())?;
                    buf.push(b'i');
                    buf.extend_from_slice(&v.to_le_bytes());
                }
            }
        }
        self.0 = buf.into_boxed_slice();
        Ok(())
    }

    fn data_start(&self) -> usize {
        let name_len = usize::from(self.0[NAME_LENGTH]);
        let cigar_ops = usize::from(u16_at(&self.0, CIGAR_OP_COUNT));
        let bases = self.sequence_length();
        FIXED_LENGTH + name_len + 4 * cigar_ops + bases.div_ceil(2) + bases
    }
}

impl fmt::Debug for BamRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.view().fmt(f)
    }
}

fn field_length(src: &[u8]) -> io::Result<usize> {
    let bad = || io::Error::new(io::ErrorKind::InvalidData, "malformed tag field");
    let ty = *src.get(2).ok_or_else(bad)?;
    let value_len = match ty {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'Z' | b'H' => src[3..].iter().position(|&b| b == 0).ok_or_else(bad)? + 1,
        b'B' => {
            let width = match *src.get(3).ok_or_else(bad)? {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => return Err(bad()),
            };
            let count = src.get(4..8).ok_or_else(bad)?;
            let count = u32::from_le_bytes(count.try_into().map_err(|_| bad())?) as usize;
            5 + width * count
        }
        _ => return Err(bad()),
    };
    let len = 3 + value_len;
    if len > src.len() {
        return Err(bad());
    }
    Ok(len)
}

fn u16_at(buf: &[u8], i: usize) -> u16 {
    u16::from_le_bytes([buf[i], buf[i + 1]])
}

fn u32_at(buf: &[u8], i: usize) -> u32 {
    u32::from_le_bytes([buf[i], buf[i + 1], buf[i + 2], buf[i + 3]])
}

fn i32_at(buf: &[u8], i: usize) -> i32 {
    i32::from_le_bytes([buf[i], buf[i + 1], buf[i + 2], buf[i + 3]])
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_core::Position;
    use noodles_sam::alignment::record::cigar::op::{Kind, Op};
    use noodles_sam::alignment::record::data::field::Tag;
    use noodles_sam::alignment::record_buf::data::field::Value;
    use noodles_sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
    use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    use std::num::NonZeroUsize;

    fn header() -> Header {
        Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(10_000).unwrap()),
            )
            .build()
    }

    fn record_buf() -> RecordBuf {
        let mut r = RecordBuf::default();
        *r.name_mut() = Some("read1".into());
        *r.flags_mut() = Flags::REVERSE_COMPLEMENTED;
        *r.reference_sequence_id_mut() = Some(0);
        *r.alignment_start_mut() = Position::new(101);
        *r.mapping_quality_mut() = noodles_sam::alignment::record::MappingQuality::new(37);
        *r.cigar_mut() = [Op::new(Kind::Match, 4)].into_iter().collect::<Cigar>();
        *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
        *r.quality_scores_mut() = QualityScores::from(vec![30, 30, 30, 30]);
        r.data_mut()
            .insert(Tag::from(*b"BX"), Value::String("OLD".into()));
        r.data_mut().insert(Tag::from(*b"NH"), Value::UInt8(1));
        r
    }

    #[test]
    fn fields_read_back_as_pysam_reports_them() {
        let r = BamRecord::from_alignment_record(&header(), &record_buf()).unwrap();
        assert_eq!(r.name(), Some(&b"read1"[..]));
        assert_eq!(r.reference_id(), Some(0));
        assert_eq!(r.position(), Some(100));
        assert_eq!(r.mapping_quality(), 37);
        assert!(r.flags().is_reverse_complemented());
        assert_eq!(r.mate_reference_id(), None);
        assert_eq!(r.mate_position(), None);
        assert_eq!(r.sequence_length(), 4);
        assert_eq!(r.to_record_buf(&header()).unwrap(), record_buf());
    }

    #[test]
    fn a_missing_mapping_quality_reads_back_as_255() {
        let mut rb = record_buf();
        *rb.mapping_quality_mut() = None;
        let r = BamRecord::from_alignment_record(&header(), &rb).unwrap();
        assert_eq!(r.mapping_quality(), 255);
    }

    #[test]
    fn records_round_trip_through_their_length_prefix() {
        let r = BamRecord::from_alignment_record(&header(), &record_buf()).unwrap();
        let mut buf = Vec::new();
        r.write_to(&mut buf).unwrap();
        r.write_to(&mut buf).unwrap();
        let mut src = &buf[..];
        assert_eq!(BamRecord::read_from(&mut src).unwrap().as_ref(), Some(&r));
        assert_eq!(BamRecord::read_from(&mut src).unwrap().as_ref(), Some(&r));
        assert!(BamRecord::read_from(&mut src).unwrap().is_none());
        assert!(BamRecord::read_from(&mut &buf[..2]).is_err());
    }

    #[test]
    fn truncated_records_are_rejected() {
        let r = BamRecord::from_alignment_record(&header(), &record_buf()).unwrap();
        let bytes = r.as_bytes();
        assert!(BamRecord::from_bytes(bytes[..bytes.len() - 20].to_vec()).is_err());
        assert!(BamRecord::from_bytes(vec![0; 10]).is_err());
    }

    #[test]
    fn set_tag_matches_pysam() {
        let mut r = BamRecord::from_alignment_record(&header(), &record_buf()).unwrap();
        r.set_tag(*b"UG", TagValue::Int(5)).unwrap();
        r.set_tag(*b"BX", TagValue::String(b"NEW")).unwrap();
        // BX moved to the end, as pysam's bam_aux_del + bam_aux_append does.
        let data = r.view().data();
        let tags: Vec<[u8; 2]> = data.iter().map(|f| *f.unwrap().0.as_ref()).collect();
        assert_eq!(tags, vec![*b"NH", *b"UG", *b"BX"]);
        let rb = r.to_record_buf(&header()).unwrap();
        assert_eq!(rb.data().get(&Tag::from(*b"UG")), Some(&Value::UInt8(5)));
        assert_eq!(
            rb.data().get(&Tag::from(*b"BX")),
            Some(&Value::String("NEW".into()))
        );
    }

    #[test]
    fn integer_tags_use_the_smallest_type() {
        let mut r = BamRecord::from_alignment_record(&header(), &record_buf()).unwrap();
        for (n, expected) in [
            (0, Value::UInt8(0)),
            (255, Value::UInt8(255)),
            (256, Value::UInt16(256)),
            (70_000, Value::UInt32(70_000)),
            (-1, Value::Int8(-1)),
            (-200, Value::Int16(-200)),
            (-40_000, Value::Int32(-40_000)),
        ] {
            r.set_tag(*b"UG", TagValue::Int(n)).unwrap();
            let rb = r.to_record_buf(&header()).unwrap();
            assert_eq!(rb.data().get(&Tag::from(*b"UG")), Some(&expected), "{n}");
        }
        assert!(r.set_tag(*b"UG", TagValue::Int(1 << 40)).is_err());
    }

    #[test]
    fn the_sort_key_orders_like_samtools() {
        let at = |tid: Option<usize>, start: Option<usize>, reverse: bool| {
            let mut rb = RecordBuf::default();
            *rb.flags_mut() = if reverse {
                Flags::REVERSE_COMPLEMENTED
            } else {
                Flags::empty()
            };
            *rb.reference_sequence_id_mut() = tid;
            *rb.alignment_start_mut() = start.and_then(Position::new);
            BamRecord::from_alignment_record(&header(), &rb)
                .unwrap()
                .coordinate_key()
        };
        assert!(at(Some(0), Some(5), false) < at(Some(0), Some(5), true));
        assert!(at(Some(0), Some(5), true) < at(Some(0), Some(6), false));
        assert!(at(Some(0), None, false) < at(Some(0), Some(1), false));
        assert!(at(Some(0), Some(9_000), true) < at(None, None, false));
    }
}
