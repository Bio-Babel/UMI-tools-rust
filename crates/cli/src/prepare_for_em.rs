//! `umi_tools prepare-for-em` — make output usable by RSEM and Salmon
//! (Python source: `umi_tools/prepare_for_em.py`).
//!
//! `dedup` and `group` output one read 2 per selected read 1, which is correct
//! but leaves some read 1s pointing at a mate that is not in the file. RSEM and
//! Salmon reject that. This command rewrites the file so every read 1 has
//! exactly one read 2 immediately after it, and the two point at each other.
//!
//! The input must be **name-sorted** (`samtools sort -n` or `samtools
//! collate`), because the grouping is a streaming group-by on the read name.

use std::collections::BTreeMap;

use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::RecordBuf;

use umi_tools_core::alignment::open_output_alignments;
use umi_tools_core::utilities::{self, determine_format};
use umi_tools_core::{PortError, Result};

use crate::documentation;
use crate::options::{self, OptionGroups};

/// The subcommand's short description.
pub const ABOUT: &str = "prepare_for_em - make output from dedup or group compatible with EM tools";

/// The tags copied from read 1 to read 2 by default
/// (Python source: `default="UG,BX"`).
pub const DEFAULT_TAGS: &str = "UG,BX";

/// Groups consecutive records sharing a read name
/// (Python source: `prepare_for_em.chunk_bam`).
///
/// A streaming group-by, so the input must be name-sorted; an unsorted input
/// silently splits a template across several groups, exactly as in Python.
///
/// # Examples
///
/// ```
/// use noodles_sam::alignment::RecordBuf;
/// use umi_tools_cli::prepare_for_em::chunk_by_name;
///
/// let mut a = RecordBuf::default();
/// *a.name_mut() = Some("r1".into());
/// let mut b = RecordBuf::default();
/// *b.name_mut() = Some("r2".into());
///
/// let chunks = chunk_by_name(&[a.clone(), a, b]);
/// assert_eq!(chunks.len(), 2);
/// assert_eq!(chunks[0].len(), 2);
/// ```
#[must_use]
pub fn chunk_by_name(records: &[RecordBuf]) -> Vec<Vec<RecordBuf>> {
    let mut chunks: Vec<Vec<RecordBuf>> = Vec::new();
    let mut current: Vec<RecordBuf> = Vec::new();
    let mut last_name: Option<Vec<u8>> = None;

    for record in records {
        let name = record.name().map(|n| n.to_vec());
        if last_name.is_some() && last_name != name && !current.is_empty() {
            chunks.push(std::mem::take(&mut current));
        }
        last_name = name;
        current.push(record.clone());
    }
    // Python yields the final buffer unconditionally, even when empty; an
    // empty trailing group carries no records, so it is simply not emitted.
    if !current.is_empty() {
        chunks.push(current);
    }
    chunks
}

/// Copies the named tags from `source` onto `target`
/// (Python source: `prepare_for_em.copy_tags`).
///
/// A tag absent from `source` is skipped rather than cleared, matching
/// Python's `except KeyError: pass`.
#[must_use]
pub fn copy_tags(tags: &[[u8; 2]], source: &RecordBuf, mut target: RecordBuf) -> RecordBuf {
    for tag in tags {
        if let Some(value) = source.data().get(&Tag::from(*tag)) {
            target.data_mut().insert(Tag::from(*tag), value.clone());
        }
    }
    target
}

/// Picks `read`'s mate from the candidates filed under `mate_key`
/// (Python source: `prepare_for_em.pick_mate`).
///
/// Prefers a candidate that names `read` as *its* mate; failing that, takes
/// the first candidate at the position. Python's comment explains the
/// fallback: it happens when UMI-tools output the wrong read as something's
/// pair.
#[must_use]
pub fn pick_mate<'a>(read: &RecordBuf, candidates: &'a [RecordBuf]) -> Option<&'a RecordBuf> {
    let points_back = candidates.iter().find(|candidate| {
        candidate.mate_reference_sequence_id() == read.reference_sequence_id()
            && candidate.mate_alignment_start() == read.alignment_start()
    });
    points_back.or_else(|| candidates.first())
}

/// The key a record is filed under within a template
/// (Python source: `(read.reference_name, read.pos, not read.is_secondary)`).
type MateKey = (Option<usize>, Option<usize>, bool);

fn mate_key(record: &RecordBuf) -> MateKey {
    (
        record.reference_sequence_id(),
        record.alignment_start().map(usize::from),
        !record.flags().is_secondary(),
    )
}

/// Runs `umi_tools prepare-for-em` (Python source: `prepare_for_em.main`).
///
/// # Errors
///
/// Returns [`umi_tools_core::PortError::Value`] for an invalid option, and
/// [`umi_tools_core::PortError::Alignment`] or [`umi_tools_core::PortError::Io`] for input or output problems.
pub fn main(args: &[String], argv: &[String]) -> Result<i32> {
    if crate::wants_help(args) {
        // `optparse`'s layout, reproduced verbatim; see `crate::help`.
        print!("{}", crate::help::PREPARE_FOR_EM);
        return Ok(0);
    }

    let started = std::time::Instant::now();
    let groups = OptionGroups {
        group_dedup: false,
        umi_grouping: false,
        sam: false,
        sbcram: true,
        ..OptionGroups::default()
    };

    let command = options::build_parser("prepare-for-em", ABOUT, groups)
        .after_help(format!(
            "{}{}",
            documentation::GENERIC_DOCSTRING_SBCRAM_INPUT,
            documentation::GENERIC_DOCSTRING_SBCRAM_OUTPUT
        ))
        .arg(
            clap::Arg::new("tags")
                .overrides_with("tags")
                .long("tags")
                .num_args(1)
                .help("comma-separated list of tags to transfer from read1 to read2"),
        )
        .arg(
            clap::Arg::new("sam")
                .overrides_with("sam")
                .long("sam")
                .action(clap::ArgAction::SetTrue)
                .help("input and output SAM rather than BAM"),
        );

    let matches = match command.try_get_matches_from(
        std::iter::once("prepare-for-em".to_string()).chain(args.iter().cloned()),
    ) {
        Ok(m) => m,
        Err(e) => {
            e.print().ok();
            return Ok(i32::from(e.use_stderr()));
        }
    };

    let io = options::io_options(&matches)?;
    let run_id = options::start(&io, argv)?;

    let tags: Vec<[u8; 2]> = matches
        .get_one::<String>("tags")
        .map_or(DEFAULT_TAGS, String::as_str)
        .split(',')
        .filter(|t| !t.is_empty())
        .map(|t| {
            let b = t.as_bytes();
            if b.len() == 2 {
                Ok([b[0], b[1]])
            } else {
                Err(PortError::value(format!(
                    "SAM tags are two characters; got '{t}'"
                )))
            }
        })
        .collect::<Result<Vec<_>>>()?;

    let in_path = io.stdin.as_deref().ok_or_else(|| {
        PortError::value(
            "reading an alignment from standard input is not supported; use -I/--stdin",
        )
    })?;
    let sam_flag = matches.get_flag("sam");
    let in_format = determine_format(&in_path.to_string_lossy(), sam_flag, io.in_format);
    let mut opened = umi_tools_core::alignment::open_input_alignments(
        in_path,
        in_format,
        io.reference_filename.as_deref(),
    )?;
    let records = opened.read_all()?;
    let header = opened.header.clone();

    let out_name = io
        .stdout
        .as_ref()
        .map_or_else(|| "-".to_string(), |p| p.to_string_lossy().into_owned());
    let out_format = determine_format(&out_name, io.out_sam, io.out_format);
    let mut writer = match io.stdout.as_deref() {
        Some(path) => {
            open_output_alignments(path, out_format, &header, io.reference_filename.as_deref())?
        }
        None => umi_tools_core::alignment::open_output_alignments_to_stdout(
            out_format,
            &header,
            io.reference_filename.as_deref(),
        )?,
    };

    let mut pairs_output = 0u64;
    let mut no_mate = 0u64;
    let mut not_read12 = 0u64;

    for template in chunk_by_name(&records) {
        // Read 1s and read 2s indexed separately, as Python's
        // `current_template[read.is_read1]` does.
        let mut by_read1: BTreeMap<MateKey, Vec<RecordBuf>> = BTreeMap::new();
        let mut by_read2: BTreeMap<MateKey, Vec<RecordBuf>> = BTreeMap::new();
        for record in &template {
            let target = if is_read1(record) {
                &mut by_read1
            } else {
                &mut by_read2
            };
            target
                .entry(mate_key(record))
                .or_default()
                .push(record.clone());
        }

        let mut emitted: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();

        for record in &template {
            let opposite = if is_read1(record) {
                &by_read2
            } else {
                &by_read1
            };
            let base = (
                record.mate_reference_sequence_id(),
                record.mate_alignment_start().map(usize::from),
            );

            // Prefer a mate with the same primary/secondary status, then the
            // opposite — Python's two-step lookup.
            let mate = [
                record.flags().is_secondary(),
                !record.flags().is_secondary(),
            ]
            .into_iter()
            .find_map(|secondary| {
                opposite
                    .get(&(base.0, base.1, !secondary))
                    .and_then(|candidates| pick_mate(record, candidates))
            })
            .cloned();

            let Some(mut mate) = mate else {
                no_mate += 1;
                utilities::warn(&format!(
                    "Alignment {} has no mate -- skipped",
                    record
                        .name()
                        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                        .unwrap_or_default()
                ));
                continue;
            };

            let read = record.clone();
            // A secondary read's mate is made secondary too; the reverse is
            // not done, or the file would gain a second primary alignment.
            if read.flags().is_secondary() {
                *mate.flags_mut() |= noodles_sam::alignment::record::Flags::SECONDARY;
            }

            let (first, second) = if is_read1(&read) {
                let mate = copy_tags(&tags, &read, mate);
                (read, mate)
            } else if read.flags().is_segmented() && read.flags().is_last_segment() {
                let read = copy_tags(&tags, &mate, read);
                (mate, read)
            } else {
                not_read12 += 1;
                utilities::warn(&format!(
                    "Alignment {} is neither read1 nor read2 -- skipped",
                    read.name()
                        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                        .unwrap_or_default()
                ));
                continue;
            };

            // Each pair is reached twice — once from each side — so an
            // identity key keeps it from being written twice.
            let key = format!("{first:?}{second:?}");
            if emitted.insert(key) {
                writer.write(&first)?;
                writer.write(&second)?;
                pairs_output += 1;
            }
        }
    }
    writer.finish()?;

    utilities::info(&format!(
        "Total pairs output: {pairs_output}, Pairs skipped - no mates: {no_mate}, \
         Pairs skipped - not read1 or 2: {not_read12}"
    ));
    options::stop(&io, &run_id, started.elapsed().as_secs())?;
    Ok(0)
}

/// Whether a record is read 1 of a pair.
///
/// An unpaired record is neither read 1 nor read 2 in Python's terms, because
/// `pysam`'s `is_read1` requires the segmented flag.
fn is_read1(record: &RecordBuf) -> bool {
    let flags = record.flags();
    flags.is_segmented() && flags.is_first_segment()
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_core::Position;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record_buf::data::field::Value;

    fn record(name: &str, read1: bool, start: usize, mate_start: usize) -> RecordBuf {
        let mut r = RecordBuf::default();
        *r.flags_mut() = Flags::SEGMENTED
            | if read1 {
                Flags::FIRST_SEGMENT
            } else {
                Flags::LAST_SEGMENT
            };
        *r.name_mut() = Some(name.into());
        *r.reference_sequence_id_mut() = Some(0);
        *r.alignment_start_mut() = Position::new(start);
        *r.mate_reference_sequence_id_mut() = Some(0);
        *r.mate_alignment_start_mut() = Position::new(mate_start);
        r
    }

    #[test]
    fn about_names_the_command() {
        assert!(ABOUT.starts_with("prepare_for_em - "));
    }

    #[test]
    fn chunking_groups_consecutive_names() {
        let a = record("r1", true, 10, 20);
        let b = record("r1", false, 20, 10);
        let c = record("r2", true, 30, 40);
        let chunks = chunk_by_name(&[a, b, c]);
        assert_eq!(chunks.len(), 2);
        assert_eq!(chunks[0].len(), 2);
        assert_eq!(chunks[1].len(), 1);
    }

    #[test]
    fn chunking_an_empty_input_yields_nothing() {
        assert!(chunk_by_name(&[]).is_empty());
    }

    #[test]
    fn tags_are_copied_only_when_present() {
        let mut source = record("r1", true, 10, 20);
        source.data_mut().insert(Tag::from(*b"UG"), Value::Int32(7));
        let target = record("r1", false, 20, 10);

        let copied = copy_tags(&[*b"UG", *b"BX"], &source, target);
        assert!(copied.data().get(&Tag::from(*b"UG")).is_some());
        assert!(copied.data().get(&Tag::from(*b"BX")).is_none());
    }

    #[test]
    fn mate_selection_prefers_a_read_pointing_back() {
        let read = record("r1", true, 10, 20);
        // The first candidate does not point back; the second does.
        let wrong = record("r1", false, 20, 99);
        let right = record("r1", false, 20, 10);
        let candidates = vec![wrong, right.clone()];
        assert_eq!(
            pick_mate(&read, &candidates)
                .unwrap()
                .mate_alignment_start(),
            right.mate_alignment_start()
        );
    }

    #[test]
    fn mate_selection_falls_back_to_the_first_candidate() {
        let read = record("r1", true, 10, 20);
        let only = record("r1", false, 20, 99);
        let candidates = vec![only.clone()];
        assert_eq!(
            pick_mate(&read, &candidates).unwrap().alignment_start(),
            only.alignment_start()
        );
        assert!(pick_mate(&read, &[]).is_none());
    }

    #[test]
    fn default_tags_match_python() {
        assert_eq!(DEFAULT_TAGS, "UG,BX");
    }
}
