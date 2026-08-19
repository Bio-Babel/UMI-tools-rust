//! Logging, file access, temporary files and output-format helpers
//! (Python source: `umi_tools/Utilities.py`).
//!
//! `Utilities.py` mixes two unrelated concerns: a library half (logging,
//! `openFile`, temp files, format determination) and an `optparse` half. Only
//! the library half lives here; the option-parsing machinery is in
//! `crates/cli`, because a library has no business owning a command-line
//! parser. See `port_reports/umi_tools/05_design.md` § Crate layout.

use std::fmt::Write as _;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::{OnceLock, RwLock};

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::error::{PortError, Result};

/// Log verbosity, mirroring `Utilities.py`'s `-v/--verbose` integer scale
/// and its mapping onto Python `logging` levels in `Start`.
///
/// Python maps `loglevel == 0` to `logging.ERROR`, `== 1` to `logging.INFO`
/// and anything greater to `logging.DEBUG`; there is no level between INFO and
/// DEBUG, so values above 2 are indistinguishable.
///
/// The variants are ordered quietest-first (`Error < Info < Debug`), so a
/// record is emitted when its own level is **not greater** than the
/// configured one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum LogLevel {
    /// Only `error`/`critical` messages (`--verbose=0`).
    Error,
    /// Adds `info` and `warning` (`--verbose=1`, the default for most tools).
    Info,
    /// Adds `debug` (`--verbose` greater than 1).
    Debug,
}

impl LogLevel {
    /// Maps the integer `--verbose` value onto a level, exactly as
    /// `Utilities.Start` does.
    #[must_use]
    pub fn from_verbosity(loglevel: i32) -> Self {
        match loglevel {
            i32::MIN..=0 => LogLevel::Error,
            1 => LogLevel::Info,
            _ => LogLevel::Debug,
        }
    }
}

/// Where log records are written.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LogSink {
    /// Standard output (Python's default `stdlog`).
    Stdout,
    /// Standard error (`--log2stderr`).
    Stderr,
    /// A file opened in append mode (`-L/--log`); Python opens it with mode `"a"`.
    File(PathBuf),
}

/// Process-global logging configuration.
///
/// Mirrors `Utilities.DefaultOptions` plus the state `Start` mutates on it.
/// It is the only genuinely global mutable state in the port, and it exists
/// for the same reason it does in Python: `info()` is called from deep inside
/// library code and must honour `--log`/`--verbose` without threading a
/// logger through forty functions.
#[derive(Debug, Clone)]
pub struct GlobalOptions {
    /// Current verbosity.
    pub loglevel: LogLevel,
    /// Where log records go.
    pub sink: LogSink,
    /// Whether log lines are prefixed with `# `.
    ///
    /// Python chooses `'# %(asctime)s %(levelname)s %(message)s'` when
    /// `stdout` and `stdlog` are the same stream, so the log does not corrupt
    /// the data stream it shares, and the unprefixed form otherwise.
    pub comment_prefix: bool,
    /// gzip compression level for `.gz` outputs (Python default: 6).
    pub compresslevel: u32,
}

impl Default for GlobalOptions {
    fn default() -> Self {
        Self {
            loglevel: LogLevel::Info,
            sink: LogSink::Stdout,
            comment_prefix: true,
            compresslevel: 6,
        }
    }
}

fn global() -> &'static RwLock<GlobalOptions> {
    static G: OnceLock<RwLock<GlobalOptions>> = OnceLock::new();
    G.get_or_init(|| RwLock::new(GlobalOptions::default()))
}

/// Replaces the process-global logging configuration.
///
/// Called once by `crates/cli`'s `start`, matching the point at which
/// `Utilities.Start` calls `logging.basicConfig`.
pub fn set_global_options(options: GlobalOptions) {
    if let Ok(mut guard) = global().write() {
        *guard = options;
    }
}

/// Returns a copy of the current process-global logging configuration.
#[must_use]
pub fn global_options() -> GlobalOptions {
    global().read().map(|g| g.clone()).unwrap_or_default()
}

/// Indents continuation lines of a multi-line log record
/// (Python source: `Utilities.MultiLineFormatter.format`).
///
/// Python's implementation is:
///
/// ```text
/// s = logging.Formatter.format(self, record)
/// prefix = "#" if s.startswith("#") else ""
/// header, footer = s.split(record.message)
/// s = prefix + s.replace('\n', '\n%s' % prefix + ' ' * len(header))
/// ```
///
/// Two consequences that are easy to miss, both reproduced here because both
/// are visible in the log files `umi_tools` writes:
///
/// 1. The continuation indent is `len(header)` spaces *after* the `#`, so a
///    continuation line is one character longer than its header.
/// 2. `prefix` is prepended to the **whole** record, and the record already
///    begins with `#` whenever the format string does — so the first line of
///    every record comes out with a **doubled** `##`. Confirmed by running
///    `umi_tools extract` with the log going to stdout, which emits
///    `## 2026-08-18 20:28:33,938 INFO Starting barcode extraction`.
///
/// The doubling happens only when log and data share a stream; with
/// `-L/--log` pointing at a file the format string carries no `#`, `prefix` is
/// empty, and lines are unprefixed.
///
/// # Examples
///
/// ```
/// use umi_tools_core::utilities::format_multiline;
/// let header = "# 12:00:00 INFO ";
/// let out = format_multiline(header, "first\nsecond");
/// let indent = format!("#{}", " ".repeat(header.len()));
/// assert_eq!(out, format!("##{} first\n{indent}second", &header[1..header.len() - 1]));
/// ```
#[must_use]
pub fn format_multiline(header: &str, message: &str) -> String {
    let prefix = if header.starts_with('#') { "#" } else { "" };
    let pad = " ".repeat(header.len());
    let mut out = String::with_capacity(header.len() + message.len() + 2);
    out.push_str(prefix);
    for (i, line) in message.split('\n').enumerate() {
        if i > 0 {
            out.push('\n');
            out.push_str(prefix);
            out.push_str(&pad);
        } else {
            out.push_str(header);
        }
        out.push_str(line);
    }
    out
}

fn timestamp() -> String {
    // Python's logging emits `%(asctime)s` as `YYYY-MM-DD HH:MM:SS,mmm`.
    // Reproducing the wall-clock value is not parity-relevant (it differs on
    // every run in Python too, and the test harness strips it), but the shape
    // is, because downstream log readers split on it.
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default();
    let secs = now.as_secs();
    let millis = now.subsec_millis();
    let (days, rem) = (secs / 86_400, secs % 86_400);
    let (h, m, s) = (rem / 3600, (rem % 3600) / 60, rem % 60);
    let (y, mo, d) = civil_from_days(days as i64);
    format!("{y:04}-{mo:02}-{d:02} {h:02}:{m:02}:{s:02},{millis:03}")
}

/// Converts days since the Unix epoch to a `(year, month, day)` triple.
///
/// Howard Hinnant's `civil_from_days` algorithm, used so the logger can format
/// a timestamp without taking a date-time dependency for what is otherwise a
/// purely cosmetic field.
fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = u32::try_from(if mp < 10 { mp + 3 } else { mp - 9 }).unwrap_or(1);
    (if m <= 2 { y + 1 } else { y }, m, d)
}

/// Formats a wall-clock time the way Python's `time.asctime` does
/// (`"Tue Aug 18 20:28:18 2026"`), which is the format `getHeader` and
/// `getFooter` use — distinct from the `%(asctime)s` form the logger uses.
fn asctime(secs: u64) -> String {
    const DAYS: [&str; 7] = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"];
    const MONTHS: [&str; 12] = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ];
    let (days, rem) = (secs / 86_400, secs % 86_400);
    let (h, m, sec) = (rem / 3600, (rem % 3600) / 60, rem % 60);
    let (y, mo, d) = civil_from_days(days as i64);
    // 1 Jan 1970 was a Thursday, index 4 in `DAYS`.
    let weekday = DAYS[((days + 4) % 7) as usize];
    let month = MONTHS[(mo as usize).saturating_sub(1).min(11)];
    format!("{weekday} {month} {d:2} {h:02}:{m:02}:{sec:02} {y:04}")
}

/// Seconds since the Unix epoch, for the header and footer timestamps.
fn now_secs() -> u64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs()
}

fn emit(level: LogLevel, label: &str, message: &str) {
    let opts = global_options();
    // `LogLevel` is ordered by *verbosity*, not by severity: `Error` is the
    // quietest setting and `Debug` the loudest. A record is emitted when its
    // level is no louder than the configured one, so raising the setting adds
    // records rather than removing them. Getting this comparison backwards
    // silently swallows every error and warning — which is exactly what
    // happened until the end-to-end CLI suite caught it.
    if level > opts.loglevel {
        return;
    }
    let header = if opts.comment_prefix {
        format!("# {} {} ", timestamp(), label)
    } else {
        format!("{} {} ", timestamp(), label)
    };
    let line = format_multiline(&header, message);
    match &opts.sink {
        LogSink::Stdout => {
            let mut out = std::io::stdout().lock();
            let _ = writeln!(out, "{line}");
        }
        LogSink::Stderr => {
            let mut out = std::io::stderr().lock();
            let _ = writeln!(out, "{line}");
        }
        LogSink::File(path) => {
            if let Ok(f) = OpenOptions::new().create(true).append(true).open(path) {
                let mut out = BufWriter::new(f);
                let _ = writeln!(out, "{line}");
            }
        }
    }
}

/// Logs `message` at an explicit level (Python source: `Utilities.log`).
pub fn log(level: LogLevel, message: &str) {
    let label = match level {
        LogLevel::Error => "ERROR",
        LogLevel::Info => "INFO",
        LogLevel::Debug => "DEBUG",
    };
    emit(level, label, message);
}

/// Logs an informational message (Python source: `Utilities.info`).
pub fn info(message: &str) {
    emit(LogLevel::Info, "INFO", message);
}

/// Logs a warning (Python source: `Utilities.warning`).
pub fn warning(message: &str) {
    emit(LogLevel::Info, "WARNING", message);
}

/// Logs a warning. Alias of [`warning`], as in the Python original, which
/// defines both `warn` and `warning` with identical bodies.
pub fn warn(message: &str) {
    warning(message);
}

/// Logs a debug message (Python source: `Utilities.debug`).
pub fn debug(message: &str) {
    emit(LogLevel::Debug, "DEBUG", message);
}

/// Logs a critical message (Python source: `Utilities.critical`).
pub fn critical(message: &str) {
    emit(LogLevel::Error, "CRITICAL", message);
}

/// Logs an error and returns the fatal-error notice the Python original writes
/// to stderr before exiting.
///
/// `Utilities.error` calls `sys.exit(1)` after logging, which a library must
/// not do. This performs the logging half and hands the caller the message;
/// `crates/cli` writes it to stderr and exits 1, so the CLI's observable
/// behaviour is unchanged. See `05_design.md` deviation D6.
///
/// # Examples
///
/// ```
/// use umi_tools_core::utilities::error;
/// assert_eq!(error("bad input"), "UMI-tools failed with an error. Check the log file\n");
/// ```
#[must_use]
pub fn error(message: &str) -> &'static str {
    emit(LogLevel::Error, "ERROR", message);
    "UMI-tools failed with an error. Check the log file\n"
}

/// A reader over a plain or gzip-compressed file.
///
/// Python's `openFile` returns either a `gzip` file object or a builtin file
/// object depending on the extension; Rust needs one concrete type, so the two
/// cases are an enum rather than a `Box<dyn BufRead>` — callers can then name
/// the type in a struct field without a virtual call per byte.
pub enum InputFile {
    /// An uncompressed file.
    Plain(BufReader<File>),
    /// A gzip-compressed file (`.gz` or `.z`).
    Gzip(Box<BufReader<MultiGzDecoder<File>>>),
}

impl std::io::Read for InputFile {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            InputFile::Plain(r) => std::io::Read::read(r, buf),
            InputFile::Gzip(r) => std::io::Read::read(r, buf),
        }
    }
}

impl BufRead for InputFile {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        match self {
            InputFile::Plain(r) => r.fill_buf(),
            InputFile::Gzip(r) => r.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            InputFile::Plain(r) => r.consume(amt),
            InputFile::Gzip(r) => r.consume(amt),
        }
    }
}

/// A writer to a plain or gzip-compressed file.
pub enum OutputFile {
    /// An uncompressed file.
    Plain(BufWriter<File>),
    /// A gzip-compressed file (`.gz` or `.z`).
    Gzip(Box<GzEncoder<BufWriter<File>>>),
}

impl OutputFile {
    /// Finishes the gzip stream, if any, and flushes to disk.
    ///
    /// A `GzEncoder` must be finished to emit its trailer; dropping it flushes
    /// but does not guarantee the trailer is written before the file handle
    /// closes. Python's `gzip` file objects do this in `close()`.
    ///
    /// # Errors
    ///
    /// Returns [`PortError::Io`] if the final write fails.
    pub fn finish(self) -> Result<()> {
        match self {
            OutputFile::Plain(mut w) => w
                .flush()
                .map_err(|e| PortError::io(Path::new("<output>"), e)),
            OutputFile::Gzip(w) => w
                .finish()
                .and_then(|mut inner| inner.flush())
                .map_err(|e| PortError::io(Path::new("<output>"), e)),
        }
    }
}

impl Write for OutputFile {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            OutputFile::Plain(w) => w.write(buf),
            OutputFile::Gzip(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            OutputFile::Plain(w) => w.flush(),
            OutputFile::Gzip(w) => w.flush(),
        }
    }
}

/// Whether a path names a gzip-compressed file by extension.
///
/// Python's `openFile` tests `os.path.splitext(filename)[1].lower() in
/// (".gz", ".z")`, so the check is case-insensitive and looks only at the
/// final extension.
#[must_use]
pub fn is_compressed(path: &Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("gz") || e.eq_ignore_ascii_case("z"))
}

/// Opens a file for reading, transparently decompressing `.gz`/`.z`
/// (Python source: `Utilities.openFile` with `mode="r"`).
///
/// Uses a multi-member gzip decoder, matching Python's `gzip.open`, which also
/// reads concatenated members — relevant because sequencing FASTQs are
/// routinely produced by concatenating per-lane gzip streams.
///
/// # Errors
///
/// Returns [`PortError::Io`] if the file cannot be opened.
pub fn open_file(path: &Path) -> Result<InputFile> {
    let file = File::open(path).map_err(|e| PortError::io(path, e))?;
    Ok(if is_compressed(path) {
        InputFile::Gzip(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        InputFile::Plain(BufReader::new(file))
    })
}

/// Opens a file for writing, transparently compressing `.gz`/`.z`
/// (Python source: `Utilities.openFile` with `mode="w"`).
///
/// `create_dir` mirrors Python's argument of the same name: when set, missing
/// parent directories are created first. The compression level comes from the
/// process-global [`GlobalOptions::compresslevel`], as it does in Python,
/// where `openFile` reads `global_options.compresslevel`.
///
/// # Errors
///
/// Returns [`PortError::Io`] if the parent directory cannot be created or the
/// file cannot be opened.
pub fn create_file(path: &Path, create_dir: bool) -> Result<OutputFile> {
    if create_dir {
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() && !parent.exists() {
                std::fs::create_dir_all(parent).map_err(|e| PortError::io(parent, e))?;
            }
        }
    }
    let file = File::create(path).map_err(|e| PortError::io(path, e))?;
    Ok(if is_compressed(path) {
        let level = Compression::new(global_options().compresslevel.min(9));
        OutputFile::Gzip(Box::new(GzEncoder::new(BufWriter::new(file), level)))
    } else {
        OutputFile::Plain(BufWriter::new(file))
    })
}

/// Opens a file for appending (Python source: `Utilities.openFile` is used
/// with mode `"a"` by `Start` for the log file).
///
/// # Errors
///
/// Returns [`PortError::NotImplemented`] for a `.gz` path, matching Python,
/// whose `openFile` raises `NotImplementedError` for any mode other than `"r"`
/// or `"w"` on a compressed file. Returns [`PortError::Io`] if the file cannot
/// be opened.
pub fn append_file(path: &Path) -> Result<File> {
    if is_compressed(path) {
        return Err(PortError::NotImplemented(
            "mode 'a' not implemented".to_string(),
        ));
    }
    OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .map_err(|e| PortError::io(path, e))
}

/// Creates a temporary file and returns its handle and path
/// (Python source: `Utilities.getTempFile`).
///
/// Python uses `tempfile.NamedTemporaryFile(dir=dir, delete=False,
/// prefix="ctmp", suffix=suffix)`, so the file persists after the handle is
/// dropped and the caller owns deletion. That is reproduced here: the returned
/// path is *not* removed on drop.
///
/// # Errors
///
/// Returns [`PortError::Io`] if the file cannot be created or persisted.
pub fn get_temp_file(dir: Option<&Path>, suffix: &str) -> Result<(File, PathBuf)> {
    let fallback = Path::new(".");
    let mut builder = tempfile::Builder::new();
    builder.prefix("ctmp").suffix(suffix);
    let named = match dir {
        Some(d) => builder.tempfile_in(d),
        None => builder.tempfile(),
    }
    .map_err(|e| PortError::io(dir.unwrap_or(fallback), e))?;
    let (file, path) = named
        .keep()
        .map_err(|e| PortError::io(dir.unwrap_or(fallback), std::io::Error::other(e)))?;
    Ok((file, path))
}

/// Creates a temporary file and returns only its path
/// (Python source: `Utilities.getTempFilename`).
///
/// The file exists on return and the caller owns deleting it, exactly as in
/// Python, where `getTempFilename` closes the handle and returns `.name`.
///
/// # Errors
///
/// Returns [`PortError::Io`] if the file cannot be created.
pub fn get_temp_filename(dir: Option<&Path>, suffix: &str) -> Result<PathBuf> {
    get_temp_file(dir, suffix).map(|(_, path)| path)
}

/// A SAM/BAM/CRAM alignment file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentFormat {
    /// Text SAM.
    Sam,
    /// BGZF-compressed BAM.
    Bam,
    /// Reference-compressed CRAM.
    Cram,
}

impl AlignmentFormat {
    /// The lowercase format name, matching the Python original's strings.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            AlignmentFormat::Sam => "sam",
            AlignmentFormat::Bam => "bam",
            AlignmentFormat::Cram => "cram",
        }
    }

    /// Parses a `--in-format`/`--out-format` value.
    ///
    /// # Errors
    ///
    /// Returns [`PortError::Value`] for anything other than `sam`, `bam` or
    /// `cram`.
    pub fn parse(value: &str) -> Result<Self> {
        match value.to_ascii_lowercase().as_str() {
            "sam" => Ok(AlignmentFormat::Sam),
            "bam" => Ok(AlignmentFormat::Bam),
            "cram" => Ok(AlignmentFormat::Cram),
            other => Err(PortError::value(format!(
                "unknown alignment format: {other}"
            ))),
        }
    }
}

/// Chooses the alignment format for an input or output file
/// (Python source: `Utilities.determine_format`).
///
/// Precedence, exactly as in Python: the `--in-sam`/`--out-sam` flag wins,
/// then an explicit `--in-format`/`--out-format`, then the filename extension
/// (case-insensitively `.sam` or `.cram`), and BAM otherwise. An unrecognised
/// extension yields BAM rather than an error — including for `-`
/// (stdin/stdout).
///
/// # Examples
///
/// ```
/// use umi_tools_core::utilities::{determine_format, AlignmentFormat};
/// assert_eq!(determine_format("x.CRAM", false, None), AlignmentFormat::Cram);
/// assert_eq!(determine_format("x.cram", true, None), AlignmentFormat::Sam);
/// assert_eq!(determine_format("-", false, None), AlignmentFormat::Bam);
/// ```
#[must_use]
pub fn determine_format(
    filename: &str,
    sam: bool,
    explicit: Option<AlignmentFormat>,
) -> AlignmentFormat {
    if sam {
        return AlignmentFormat::Sam;
    }
    if let Some(fmt) = explicit {
        return fmt;
    }
    // Python lowercases the whole filename before testing the suffix, so the
    // comparison is case-insensitive; the extension is compared rather than
    // the whole name so `clippy::case_sensitive_file_extension_comparisons`
    // is satisfied by construction.
    let extension = filename.rsplit('.').next().unwrap_or("");
    if filename.contains('.') && extension.eq_ignore_ascii_case("sam") {
        AlignmentFormat::Sam
    } else if filename.contains('.') && extension.eq_ignore_ascii_case("cram") {
        AlignmentFormat::Cram
    } else {
        AlignmentFormat::Bam
    }
}

/// The unsorted-output and sorted-output names and formats chosen for a run
/// (Python source: `Utilities.output_names_and_formats`' 4-tuple return).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OutputNames {
    /// Name of the (possibly temporary) unsorted output file.
    pub out_name: String,
    /// Format of the unsorted output file.
    pub out_format: AlignmentFormat,
    /// Name of the sorted output file, absent when `--no-sort-output` is set.
    pub sorted_out_name: Option<String>,
    /// Format of the sorted output file, absent when `--no-sort-output` is set.
    pub sorted_out_format: Option<AlignmentFormat>,
}

/// Decides the intermediate and final output names and formats
/// (Python source: `Utilities.output_names_and_formats`).
///
/// When sorting is enabled the eventual output is written through a temporary
/// file whose format is SAM if the eventual format is SAM and BAM otherwise —
/// CRAM is never used for the intermediate, because a CRAM must be sorted to
/// be useful and the intermediate by definition is not.
///
/// `eventual_name` is `"-"` for standard output, as Python passes when
/// `options.stdout is sys.stdout`.
///
/// # Errors
///
/// Returns [`PortError::Io`] if the temporary file cannot be created.
pub fn output_names_and_formats(
    eventual_name: &str,
    sam: bool,
    explicit: Option<AlignmentFormat>,
    no_sort_out: bool,
    tmpdir: Option<&Path>,
) -> Result<OutputNames> {
    let eventual_format = determine_format(eventual_name, sam, explicit);
    if no_sort_out {
        Ok(OutputNames {
            out_name: eventual_name.to_string(),
            out_format: eventual_format,
            sorted_out_name: None,
            sorted_out_format: None,
        })
    } else {
        let tmp = get_temp_filename(tmpdir, "")?;
        Ok(OutputNames {
            out_name: tmp.to_string_lossy().into_owned(),
            out_format: if eventual_format == AlignmentFormat::Sam {
                AlignmentFormat::Sam
            } else {
                AlignmentFormat::Bam
            },
            sorted_out_name: Some(eventual_name.to_string()),
            sorted_out_format: Some(eventual_format),
        })
    }
}

/// Builds the run header written at the top of every log
/// (Python source: `Utilities.getHeader`).
///
/// Python's four lines are:
///
/// ```text
/// # UMI-tools version: 1.1.6
/// # output generated by extract --extract-method=string --stdin=t.fastq
/// # job started at Tue Aug 18 20:28:18 2026 on my-host -- <uuid>
/// # pid: 78988, system: Darwin 25.5.0 Darwin Kernel Version ... arm64
/// ```
///
/// The variable fields (time, host, pid, uuid, uname) differ on every run in
/// Python too, so they are not parity-relevant; the line structure is, because
/// the test harness strips `#`-prefixed lines by matching this shape.
///
/// `argv` is the argument vector **with the program name already removed**,
/// matching what Python has by this point: `umi_tools.py`'s `main` does
/// `del sys.argv[0]` before dispatching.
#[must_use]
pub fn get_header(argv: &[String], version: &str, run_id: &str, host: &str, uname: &str) -> String {
    let mut s = String::new();
    let _ = writeln!(s, "# UMI-tools version: {version}");
    let _ = writeln!(s, "# output generated by {}", argv.join(" "));
    let _ = writeln!(
        s,
        "# job started at {} on {host} -- {run_id}",
        asctime(now_secs())
    );
    let _ = write!(s, "# pid: {}, system: {uname}", std::process::id());
    s
}

/// Renders the resolved options as `#`-prefixed lines
/// (Python source: `Utilities.getParams`).
///
/// Python sorts the option dictionary by key and formats each as
/// `"# %-40s: %s"`, or emits `"# no parameters."` when there is nothing to
/// report. Both behaviours are reproduced, including the 40-column left
/// padding, because the log files the test suite compares contain these lines.
///
/// # Examples
///
/// ```
/// use umi_tools_core::utilities::get_params;
/// assert_eq!(get_params(&[]), "# no parameters.");
/// let out = get_params(&[("method".into(), "directional".into())]);
/// assert!(out.starts_with("# method"));
/// assert!(out.ends_with(": directional"));
/// ```
#[must_use]
pub fn get_params(params: &[(String, String)]) -> String {
    if params.is_empty() {
        return "# no parameters.".to_string();
    }
    let mut sorted = params.to_vec();
    sorted.sort_by(|a, b| a.0.cmp(&b.0));
    sorted
        .iter()
        .map(|(k, v)| format!("# {k:<40}: {v}"))
        .collect::<Vec<_>>()
        .join("\n")
}

/// Builds the run footer written when a tool finishes
/// (Python source: `Utilities.getFooter`).
///
/// Python's line is:
///
/// ```text
/// # job finished in 0 seconds at Tue Aug 18 20:28:18 2026 --  1.19  0.11  0.00  0.00 -- <uuid>
/// ```
///
/// The four floats are `os.times()[:4]` — user, system, child-user and
/// child-system CPU seconds — each formatted as `%5.2f`. `cpu_times` carries
/// them in the same order; a caller with no way to measure them can pass
/// zeros, which keeps the line's shape without inventing figures.
#[must_use]
pub fn get_footer(elapsed_secs: u64, run_id: &str, cpu_times: [f64; 4]) -> String {
    let times = cpu_times
        .iter()
        .map(|t| format!("{t:5.2}"))
        .collect::<Vec<_>>()
        .join(" ");
    format!(
        "# job finished in {elapsed_secs} seconds at {} -- {times} -- {run_id}",
        asctime(now_secs())
    )
}

/// Decodes bytes as UTF-8, mirroring Python's `bytes.decode()` default
/// `strict` error handling.
///
/// Barcodes travel through this crate as raw bytes (see `05_design.md`
/// § Ownership model); this is the single boundary at which they become text
/// for output, so invalid UTF-8 fails here rather than silently producing
/// replacement characters downstream.
///
/// # Errors
///
/// Returns [`PortError::Value`] for invalid UTF-8, as Python raises
/// `UnicodeDecodeError` (a `ValueError` subclass).
pub fn decode_utf8(bytes: &[u8]) -> Result<String> {
    std::str::from_utf8(bytes)
        .map(str::to_owned)
        .map_err(|e| PortError::value(format!("invalid UTF-8 in barcode: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read as _;

    /// Errors and warnings must survive every verbosity setting; the
    /// comparison in `emit` was once inverted, which silently swallowed them
    /// and turned every CLI failure into a bare "check the log file".
    #[test]
    fn severity_ordering_keeps_errors_at_every_verbosity() {
        assert!(LogLevel::Error <= LogLevel::Error);
        assert!(LogLevel::Error <= LogLevel::Info);
        assert!(LogLevel::Error <= LogLevel::Debug);
        assert!(LogLevel::Debug > LogLevel::Info);
        assert!(LogLevel::Info > LogLevel::Error);
    }

    #[test]
    fn loglevel_mapping_matches_python() {
        assert_eq!(LogLevel::from_verbosity(0), LogLevel::Error);
        assert_eq!(LogLevel::from_verbosity(1), LogLevel::Info);
        assert_eq!(LogLevel::from_verbosity(2), LogLevel::Debug);
        assert_eq!(LogLevel::from_verbosity(10), LogLevel::Debug);
        assert_eq!(LogLevel::from_verbosity(-1), LogLevel::Error);
    }

    /// Python doubles the `#` on the first line of every record; see the
    /// function's own documentation for why, and the live-output evidence.
    #[test]
    fn multiline_indents_continuations() {
        let header = "# 12:00:00 INFO ";
        let out = format_multiline(header, "a\nb\nc");
        let lines: Vec<&str> = out.split('\n').collect();
        assert_eq!(lines[0], format!("#{header}a"));
        let indent = format!("#{}", " ".repeat(header.len()));
        assert_eq!(lines[1], format!("{indent}b"));
        assert_eq!(lines[2], format!("{indent}c"));
        // Python's `prefix + ' ' * len(header)` makes every continuation line
        // exactly as long as the first.
        assert_eq!(lines[1].len(), lines[0].len());
        assert_eq!(lines[2].len(), lines[0].len());
    }

    #[test]
    fn multiline_without_hash_prefix_is_not_doubled() {
        assert_eq!(
            format_multiline("12:00 INFO ", "a\nb"),
            "12:00 INFO a\n           b"
        );
    }

    #[test]
    fn multiline_single_line_still_gets_the_doubled_prefix() {
        assert_eq!(format_multiline("# H ", "solo"), "## H solo");
        assert_eq!(format_multiline("H ", "solo"), "H solo");
    }

    #[test]
    fn compression_detected_case_insensitively() {
        assert!(is_compressed(Path::new("a.gz")));
        assert!(is_compressed(Path::new("a.GZ")));
        assert!(is_compressed(Path::new("a.Z")));
        assert!(!is_compressed(Path::new("a.gzip")));
        assert!(!is_compressed(Path::new("a.fastq")));
        assert!(!is_compressed(Path::new("a")));
    }

    #[test]
    fn determine_format_precedence() {
        assert_eq!(determine_format("x.cram", true, None), AlignmentFormat::Sam);
        assert_eq!(
            determine_format("x.cram", false, Some(AlignmentFormat::Bam)),
            AlignmentFormat::Bam
        );
        assert_eq!(determine_format("x.SAM", false, None), AlignmentFormat::Sam);
        assert_eq!(
            determine_format("x.CrAm", false, None),
            AlignmentFormat::Cram
        );
        assert_eq!(determine_format("x.bam", false, None), AlignmentFormat::Bam);
        assert_eq!(determine_format("-", false, None), AlignmentFormat::Bam);
        assert_eq!(determine_format("", false, None), AlignmentFormat::Bam);
    }

    #[test]
    fn alignment_format_parse_rejects_unknown() {
        assert_eq!(AlignmentFormat::parse("SAM").unwrap(), AlignmentFormat::Sam);
        assert!(AlignmentFormat::parse("fasta").is_err());
    }

    #[test]
    fn output_names_no_sort_keeps_eventual_name() {
        let n = output_names_and_formats("out.bam", false, None, true, None).unwrap();
        assert_eq!(n.out_name, "out.bam");
        assert_eq!(n.out_format, AlignmentFormat::Bam);
        assert!(n.sorted_out_name.is_none());
        assert!(n.sorted_out_format.is_none());
    }

    #[test]
    fn output_names_sorted_uses_bam_intermediate_for_cram() {
        let n = output_names_and_formats("out.cram", false, None, false, None).unwrap();
        assert_eq!(n.out_format, AlignmentFormat::Bam);
        assert_eq!(n.sorted_out_format, Some(AlignmentFormat::Cram));
        assert_eq!(n.sorted_out_name.as_deref(), Some("out.cram"));
        let _ = std::fs::remove_file(&n.out_name);
    }

    #[test]
    fn output_names_sorted_keeps_sam_intermediate_for_sam() {
        let n = output_names_and_formats("out.sam", false, None, false, None).unwrap();
        assert_eq!(n.out_format, AlignmentFormat::Sam);
        assert_eq!(n.sorted_out_format, Some(AlignmentFormat::Sam));
        let _ = std::fs::remove_file(&n.out_name);
    }

    #[test]
    fn get_params_matches_python_layout() {
        assert_eq!(get_params(&[]), "# no parameters.");
        let out = get_params(&[
            ("zeta".to_string(), "1".to_string()),
            ("alpha".to_string(), "2".to_string()),
        ]);
        let lines: Vec<&str> = out.split('\n').collect();
        assert!(lines[0].starts_with("# alpha"));
        assert!(lines[1].starts_with("# zeta"));
        assert_eq!(lines[0].len(), "# ".len() + 40 + ": 2".len());
    }

    #[test]
    fn gzip_roundtrip_through_open_and_create() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.gz");
        {
            let mut w = create_file(&path, false).unwrap();
            w.write_all(b"hello\nworld\n").unwrap();
            w.finish().unwrap();
        }
        let mut buf = String::new();
        open_file(&path).unwrap().read_to_string(&mut buf).unwrap();
        assert_eq!(buf, "hello\nworld\n");
    }

    #[test]
    fn plain_roundtrip_through_open_and_create() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.txt");
        {
            let mut w = create_file(&path, false).unwrap();
            w.write_all(b"abc").unwrap();
            w.finish().unwrap();
        }
        let mut buf = String::new();
        open_file(&path).unwrap().read_to_string(&mut buf).unwrap();
        assert_eq!(buf, "abc");
    }

    #[test]
    fn create_file_makes_missing_parents_when_asked() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("a/b/c.txt");
        assert!(create_file(&path, false).is_err());
        assert!(create_file(&path, true).is_ok());
    }

    #[test]
    fn append_to_gzip_is_not_implemented() {
        let dir = tempfile::tempdir().unwrap();
        let err = append_file(&dir.path().join("l.gz")).unwrap_err();
        assert!(matches!(err, PortError::NotImplemented(_)));
    }

    #[test]
    fn temp_file_persists_after_handle_drop() {
        let path = get_temp_filename(None, ".tmp").unwrap();
        assert!(path.exists());
        assert!(path
            .file_name()
            .unwrap()
            .to_string_lossy()
            .starts_with("ctmp"));
        assert!(path.to_string_lossy().ends_with(".tmp"));
        std::fs::remove_file(&path).unwrap();
    }

    #[test]
    fn decode_utf8_rejects_invalid_sequences() {
        assert_eq!(decode_utf8(b"ACGT").unwrap(), "ACGT");
        assert!(decode_utf8(&[0xff, 0xfe]).is_err());
    }

    #[test]
    fn error_returns_pythons_stderr_notice() {
        assert_eq!(
            error("boom"),
            "UMI-tools failed with an error. Check the log file\n"
        );
    }

    #[test]
    fn header_matches_pythons_four_line_shape() {
        let h = get_header(
            &["extract".into(), "--bc-pattern=NNNN".into()],
            "1.1.6",
            "abc-uuid",
            "my-host",
            "Darwin 25.5.0 arm64",
        );
        let lines: Vec<&str> = h.split('\n').collect();
        assert_eq!(lines.len(), 4);
        assert_eq!(lines[0], "# UMI-tools version: 1.1.6");
        assert_eq!(lines[1], "# output generated by extract --bc-pattern=NNNN");
        assert!(lines[2].starts_with("# job started at "));
        assert!(lines[2].ends_with(" on my-host -- abc-uuid"));
        assert!(lines[3].starts_with("# pid: "));
        assert!(lines[3].ends_with(", system: Darwin 25.5.0 arm64"));
    }

    #[test]
    fn footer_matches_pythons_shape_including_cpu_times() {
        let f = get_footer(0, "abc-uuid", [1.19, 0.11, 0.0, 0.0]);
        assert!(f.starts_with("# job finished in 0 seconds at "));
        assert!(f.ends_with(" --  1.19  0.11  0.00  0.00 -- abc-uuid"));
    }

    #[test]
    fn asctime_matches_pythons_format() {
        // Python's `time.asctime(time.localtime(...))` renders in local
        // time; this renders in UTC, because taking a timezone dependency for
        // a cosmetic log field is not worth it and the test harness strips
        // these lines. The format is identical.
        assert_eq!(asctime(1_787_113_698), "Wed Aug 19 04:28:18 2026");
        // Single-digit days are space-padded, as `time.asctime` does.
        assert_eq!(&asctime(0)[..10], "Thu Jan  1");
    }
}
