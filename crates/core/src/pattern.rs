//! Barcode-extraction patterns (Python dependency: the `PyPI` `regex` module,
//! as used by `--extract-method=regex`).
//!
//! `umi_tools extract` and `umi_tools whitelist` accept a regular expression
//! whose **named groups** say which bases are UMI, which are cell barcode and
//! which to throw away: a group named `umi_1` marks UMI bases, `cell_1` marks
//! cell-barcode bases, `discard_1` marks bases to drop.
//!
//! The Rust `regex` crate cannot be used directly for this, for two reasons:
//! it has no approximate-matching extension (`{s<=n}` and friends), and it
//! provides no per-group **span**, which is exactly what
//! [`crate::extract::extract_barcodes`] needs in order to know which
//! *positions* each group covered. This module therefore implements a small
//! backtracking matcher over the pattern subset the tool actually uses.
//!
//! # The supported grammar
//!
//! ```text
//! pattern   := alternation
//! alternation := sequence ('|' sequence)*
//! sequence  := repeat*
//! repeat    := atom quantifier?
//! atom      := '.' | literal | class | group | '^' | '$'
//! group     := '(?P<' name '>' alternation ')' | '(?:' alternation ')' | '(' alternation ')'
//! class     := '[' '^'? (char | char '-' char)* ']'
//! quantifier := '*' | '+' | '?' | '{' n '}' | '{' n ',' m? '}'
//!             | '{s<=' n '}' | '{e<=' n '}' | '{i<=' n '}' | '{d<=' n '}'
//! ```
//!
//! Escapes `\d \w \s \D \W \S \. \\ \( \)` and the like are supported.
//!
//! # Fuzzy quantifiers bind to one element, not to a run
//!
//! This is the single most misread part of Python's fuzzy syntax, and getting
//! it wrong changes what the inDrop protocol's patterns match. In
//!
//! ```text
//! (?P<discard_2>GAGTGATTGCTTGTGACGCCTT{s<=2})
//! ```
//!
//! the `{s<=2}` applies **only to the final `T`**, not to the 22-base literal.
//! Confirmed against the live `regex` module: substituting the last base still
//! matches, substituting the first does not, and substituting the last *two*
//! does not either — because only one character is covered.
//!
//! Per-element fuzzy semantics, also confirmed by running rather than reading:
//!
//! | Quantifier | On a single-character atom |
//! |---|---|
//! | `{s<=n}`, `n >= 1` | matches any one character |
//! | `{e<=n}`, `n >= 1` | matches zero, one, or two characters (deletion, substitution, insertion) |
//! | `{i<=n}`, `n >= 1` | matches the character, or the character preceded/followed by extras |
//! | `{d<=n}`, `n >= 1` | matches the character, or nothing |
//!
//! A fuzzy quantifier on a multi-character atom (a group) is rejected with
//! [`PortError::Format`] rather than silently mis-matched. No pattern in the
//! package's documentation or test suite uses that form; recorded as an
//! input-form narrowing in `port_reports/umi_tools/05_design.md`.

use std::collections::BTreeMap;

use crate::error::{PortError, Result};

/// A parsed extraction pattern.
///
/// Matching is anchored at the start of the subject, like Python's
/// `.match()`, and unanchored at the end unless the pattern ends with `$`.
#[derive(Debug, Clone)]
pub struct Pattern {
    source: String,
    root: Node,
    group_names: Vec<String>,
}

/// The span a group covered, as `[start, end)` byte offsets into the subject.
pub type Span = (usize, usize);

/// What a successful match captured.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Match {
    /// The subject offsets the whole pattern covered.
    pub span: Span,
    /// Named group to the span it covered, for groups that participated.
    pub groups: BTreeMap<String, Span>,
}

impl Match {
    /// The text a named group captured, if it participated in the match.
    #[must_use]
    pub fn group<'s>(&self, name: &str, subject: &'s str) -> Option<&'s str> {
        self.groups.get(name).map(|&(s, e)| &subject[s..e])
    }

    /// The span a named group covered, if it participated in the match.
    #[must_use]
    pub fn span(&self, name: &str) -> Option<Span> {
        self.groups.get(name).copied()
    }
}

#[derive(Debug, Clone)]
enum Node {
    /// Alternatives, tried left to right.
    Alternation(Vec<Node>),
    /// Concatenation.
    Sequence(Vec<Node>),
    /// A repeated element.
    Repeat {
        node: Box<Node>,
        min: usize,
        /// `None` means unbounded.
        max: Option<usize>,
        greedy: bool,
    },
    /// A per-element approximate match.
    Fuzzy { byte: u8, kind: FuzzyKind, n: u32 },
    /// A named capturing group.
    Group {
        name: Option<String>,
        node: Box<Node>,
    },
    /// One literal byte.
    Literal(u8),
    /// Any byte except `\n`, as Python's `.` without `re.DOTALL`.
    AnyByte,
    /// A character class.
    Class {
        negated: bool,
        ranges: Vec<(u8, u8)>,
    },
    /// Start of subject.
    Start,
    /// End of subject.
    End,
    /// Matches the empty string.
    Empty,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FuzzyKind {
    Substitution,
    Edit,
    Insertion,
    Deletion,
}

struct Parser<'a> {
    bytes: &'a [u8],
    pos: usize,
    names: Vec<String>,
}

impl<'a> Parser<'a> {
    fn new(pattern: &'a str) -> Self {
        Self {
            bytes: pattern.as_bytes(),
            pos: 0,
            names: Vec::new(),
        }
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }

    fn eat(&mut self, b: u8) -> bool {
        if self.peek() == Some(b) {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    fn err(&self, msg: &str) -> PortError {
        PortError::format(format!(
            "invalid barcode pattern at offset {}: {msg}",
            self.pos
        ))
    }

    fn parse_alternation(&mut self) -> Result<Node> {
        let mut branches = vec![self.parse_sequence()?];
        while self.eat(b'|') {
            branches.push(self.parse_sequence()?);
        }
        Ok(if branches.len() == 1 {
            branches.pop().expect("checked")
        } else {
            Node::Alternation(branches)
        })
    }

    fn parse_sequence(&mut self) -> Result<Node> {
        let mut items = Vec::new();
        while let Some(b) = self.peek() {
            if b == b'|' || b == b')' {
                break;
            }
            items.push(self.parse_repeat()?);
        }
        Ok(match items.len() {
            0 => Node::Empty,
            1 => items.pop().expect("checked"),
            _ => Node::Sequence(items),
        })
    }

    fn parse_repeat(&mut self) -> Result<Node> {
        let atom = self.parse_atom()?;
        let Some(b) = self.peek() else {
            return Ok(atom);
        };
        let (min, max) = match b {
            b'*' => {
                self.pos += 1;
                (0, None)
            }
            b'+' => {
                self.pos += 1;
                (1, None)
            }
            b'?' => {
                self.pos += 1;
                (0, Some(1))
            }
            b'{' => {
                if let Some(node) = self.try_parse_fuzzy(&atom)? {
                    return Ok(node);
                }
                self.parse_counted()?
            }
            _ => return Ok(atom),
        };
        // A trailing `?` makes the quantifier lazy.
        let greedy = !self.eat(b'?');
        Ok(Node::Repeat {
            node: Box::new(atom),
            min,
            max,
            greedy,
        })
    }

    /// Parses `{s<=n}` and friends, returning `None` if this `{` starts an
    /// ordinary counted quantifier instead.
    fn try_parse_fuzzy(&mut self, atom: &Node) -> Result<Option<Node>> {
        let start = self.pos;
        self.pos += 1; // consume '{'
        let kind = match self.peek() {
            Some(b's') => FuzzyKind::Substitution,
            Some(b'e') => FuzzyKind::Edit,
            Some(b'i') => FuzzyKind::Insertion,
            Some(b'd') => FuzzyKind::Deletion,
            _ => {
                self.pos = start;
                return Ok(None);
            }
        };
        self.pos += 1;
        if !(self.eat(b'<') && self.eat(b'=')) {
            self.pos = start;
            return Ok(None);
        }
        let mut n = 0u32;
        let mut digits = 0;
        while let Some(d) = self.peek() {
            if !d.is_ascii_digit() {
                break;
            }
            n = n * 10 + u32::from(d - b'0');
            digits += 1;
            self.pos += 1;
        }
        if digits == 0 || !self.eat(b'}') {
            self.pos = start;
            return Ok(None);
        }

        let Node::Literal(byte) = atom else {
            return Err(PortError::format(
                "approximate quantifiers ({s<=n}, {e<=n}, {i<=n}, {d<=n}) are \
                 supported only on a single literal character, which is the \
                 only form UMI-tools' own patterns use"
                    .to_string(),
            ));
        };
        Ok(Some(Node::Fuzzy {
            byte: *byte,
            kind,
            n,
        }))
    }

    fn parse_counted(&mut self) -> Result<(usize, Option<usize>)> {
        self.pos += 1; // consume '{'
        let mut min = 0usize;
        let mut digits = 0;
        while let Some(d) = self.peek() {
            if !d.is_ascii_digit() {
                break;
            }
            min = min * 10 + usize::from(d - b'0');
            digits += 1;
            self.pos += 1;
        }
        if digits == 0 {
            return Err(self.err("expected a repeat count"));
        }
        if self.eat(b'}') {
            return Ok((min, Some(min)));
        }
        if !self.eat(b',') {
            return Err(self.err("expected ',' or '}' in a repeat count"));
        }
        let mut max = 0usize;
        let mut max_digits = 0;
        while let Some(d) = self.peek() {
            if !d.is_ascii_digit() {
                break;
            }
            max = max * 10 + usize::from(d - b'0');
            max_digits += 1;
            self.pos += 1;
        }
        if !self.eat(b'}') {
            return Err(self.err("unterminated repeat count"));
        }
        Ok((min, if max_digits == 0 { None } else { Some(max) }))
    }

    fn parse_atom(&mut self) -> Result<Node> {
        let Some(b) = self.peek() else {
            return Ok(Node::Empty);
        };
        match b {
            b'(' => self.parse_group(),
            b'[' => self.parse_class(),
            b'.' => {
                self.pos += 1;
                Ok(Node::AnyByte)
            }
            b'^' => {
                self.pos += 1;
                Ok(Node::Start)
            }
            b'$' => {
                self.pos += 1;
                Ok(Node::End)
            }
            b'\\' => {
                self.pos += 1;
                let Some(e) = self.peek() else {
                    return Err(self.err("trailing backslash"));
                };
                self.pos += 1;
                Ok(escape_node(e))
            }
            other => {
                self.pos += 1;
                Ok(Node::Literal(other))
            }
        }
    }

    fn parse_group(&mut self) -> Result<Node> {
        self.pos += 1; // consume '('
        let mut name = None;
        if self.eat(b'?') {
            if self.eat(b':') {
                // non-capturing
            } else if self.eat(b'P') {
                if !self.eat(b'<') {
                    return Err(self.err("expected '<' after '(?P'"));
                }
                let start = self.pos;
                while let Some(c) = self.peek() {
                    if c == b'>' {
                        break;
                    }
                    self.pos += 1;
                }
                if !self.eat(b'>') {
                    return Err(self.err("unterminated group name"));
                }
                let parsed = String::from_utf8_lossy(&self.bytes[start..self.pos - 1]).into_owned();
                self.names.push(parsed.clone());
                name = Some(parsed);
            } else if self.eat(b'<') {
                let start = self.pos;
                while let Some(c) = self.peek() {
                    if c == b'>' {
                        break;
                    }
                    self.pos += 1;
                }
                if !self.eat(b'>') {
                    return Err(self.err("unterminated group name"));
                }
                let parsed = String::from_utf8_lossy(&self.bytes[start..self.pos - 1]).into_owned();
                self.names.push(parsed.clone());
                name = Some(parsed);
            } else {
                return Err(self.err("unsupported group flag"));
            }
        }
        let inner = self.parse_alternation()?;
        if !self.eat(b')') {
            return Err(self.err("unterminated group"));
        }
        Ok(Node::Group {
            name,
            node: Box::new(inner),
        })
    }

    fn parse_class(&mut self) -> Result<Node> {
        self.pos += 1; // consume '['
        let negated = self.eat(b'^');
        let mut ranges = Vec::new();
        let mut first = true;
        loop {
            let Some(c) = self.peek() else {
                return Err(self.err("unterminated character class"));
            };
            if c == b']' && !first {
                self.pos += 1;
                break;
            }
            first = false;
            self.pos += 1;
            let lo = if c == b'\\' {
                let Some(e) = self.peek() else {
                    return Err(self.err("trailing backslash in class"));
                };
                self.pos += 1;
                match escape_node(e) {
                    Node::Literal(b) => b,
                    Node::Class { ranges: r, .. } => {
                        ranges.extend(r);
                        continue;
                    }
                    _ => return Err(self.err("unsupported escape in class")),
                }
            } else {
                c
            };
            if self.peek() == Some(b'-') && self.bytes.get(self.pos + 1) != Some(&b']') {
                self.pos += 1;
                let Some(hi) = self.peek() else {
                    return Err(self.err("unterminated range in class"));
                };
                self.pos += 1;
                ranges.push((lo, hi));
            } else {
                ranges.push((lo, lo));
            }
        }
        Ok(Node::Class { negated, ranges })
    }
}

fn escape_node(e: u8) -> Node {
    match e {
        b'd' => Node::Class {
            negated: false,
            ranges: vec![(b'0', b'9')],
        },
        b'D' => Node::Class {
            negated: true,
            ranges: vec![(b'0', b'9')],
        },
        b'w' => Node::Class {
            negated: false,
            ranges: vec![(b'a', b'z'), (b'A', b'Z'), (b'0', b'9'), (b'_', b'_')],
        },
        b'W' => Node::Class {
            negated: true,
            ranges: vec![(b'a', b'z'), (b'A', b'Z'), (b'0', b'9'), (b'_', b'_')],
        },
        b's' => Node::Class {
            negated: false,
            ranges: vec![(b' ', b' '), (b'\t', b'\r')],
        },
        b'S' => Node::Class {
            negated: true,
            ranges: vec![(b' ', b' '), (b'\t', b'\r')],
        },
        b'n' => Node::Literal(b'\n'),
        b't' => Node::Literal(b'\t'),
        b'r' => Node::Literal(b'\r'),
        other => Node::Literal(other),
    }
}

/// Captured spans accumulated during a match attempt.
type Captures = BTreeMap<String, Span>;

impl Pattern {
    /// Parses `pattern`.
    ///
    /// # Errors
    ///
    /// Returns [`PortError::Format`] for a syntactically invalid pattern, or
    /// for an approximate quantifier applied to anything but a single literal
    /// character. Python's equivalent is `regex.error`, which
    /// `Utilities.validateExtractOptions` converts into
    /// `ValueError: --bc-pattern '<pattern>' is not a valid regex`.
    ///
    /// # Examples
    ///
    /// ```
    /// use umi_tools_core::pattern::Pattern;
    /// let p = Pattern::parse("^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})")?;
    /// let m = p.match_at_start("ACGTTTTGG").expect("matches");
    /// assert_eq!(m.group("umi_1", "ACGTTTTGG"), Some("ACG"));
    /// assert_eq!(m.group("umi_2", "ACGTTTTGG"), Some("GG"));
    /// # Ok::<(), umi_tools_core::PortError>(())
    /// ```
    pub fn parse(pattern: &str) -> Result<Self> {
        let mut parser = Parser::new(pattern);
        let root = parser.parse_alternation()?;
        if parser.pos != parser.bytes.len() {
            return Err(PortError::format(format!(
                "invalid barcode pattern: unexpected ')' at offset {}",
                parser.pos
            )));
        }
        Ok(Self {
            source: pattern.to_string(),
            root,
            group_names: parser.names,
        })
    }

    /// The pattern text this was parsed from.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.source
    }

    /// Every named group in the pattern, in the order they were declared
    /// (Python source: `options.pattern.groupindex`, which
    /// `validateExtractOptions` iterates to decide whether the pattern extracts
    /// a cell barcode, a UMI, or both).
    #[must_use]
    pub fn group_names(&self) -> &[String] {
        &self.group_names
    }

    /// Matches the pattern against the start of `subject`
    /// (Python: `compiled.match(subject)`).
    ///
    /// Anchored at the start, unanchored at the end — so a pattern can match a
    /// prefix and leave the rest of the read alone, which is exactly how the
    /// 5-prime extraction patterns are written.
    ///
    /// Returns `None` if the pattern does not match.
    #[must_use]
    pub fn match_at_start(&self, subject: &str) -> Option<Match> {
        let bytes = subject.as_bytes();
        let mut captures = Captures::new();
        let mut end = None;
        let matched = match_node(&self.root, bytes, 0, &mut captures, &mut |e, _| {
            end = Some(e);
            true
        });
        if !matched {
            return None;
        }
        Some(Match {
            span: (0, end?),
            groups: captures,
        })
    }
}

/// The continuation a match step calls with the position it reached.
///
/// Returning `false` tells the step to try its next alternative — this is what
/// makes backtracking work, and it is why every combinator here threads a
/// continuation rather than returning a position.
type Cont<'a> = &'a mut dyn FnMut(usize, &mut Captures) -> bool;

/// Attempts to match `node` at `pos`, handing each candidate end position to
/// `k` until one is accepted.
///
/// Returns whether some path through `node` was accepted by `k`.
fn match_node(
    node: &Node,
    subject: &[u8],
    pos: usize,
    captures: &mut Captures,
    k: Cont<'_>,
) -> bool {
    match node {
        Node::Empty => k(pos, captures),

        Node::Start => pos == 0 && k(pos, captures),

        Node::End => pos == subject.len() && k(pos, captures),

        Node::Literal(b) => subject.get(pos) == Some(b) && k(pos + 1, captures),

        // Python's `.` excludes newline unless `re.DOTALL` is set, which
        // umi_tools never sets.
        Node::AnyByte => matches!(subject.get(pos), Some(&c) if c != b'\n') && k(pos + 1, captures),

        Node::Class { negated, ranges } => match subject.get(pos) {
            Some(&c) => {
                let inside = ranges.iter().any(|&(lo, hi)| c >= lo && c <= hi);
                inside != *negated && k(pos + 1, captures)
            }
            None => false,
        },

        Node::Fuzzy { byte, kind, n } => {
            for length in fuzzy_lengths(*kind, *n) {
                let end = pos + length;
                if end > subject.len() || !fuzzy_admits(*byte, *kind, *n, &subject[pos..end]) {
                    continue;
                }
                if k(end, captures) {
                    return true;
                }
            }
            false
        }

        Node::Group { name, node } => {
            let start = pos;
            let saved = captures.clone();
            let matched = match_node(node, subject, pos, captures, &mut |end, caps| {
                let previous = name.as_ref().and_then(|n| caps.get(n).copied());
                if let Some(name) = name {
                    caps.insert(name.clone(), (start, end));
                }
                if k(end, caps) {
                    return true;
                }
                // Roll the capture back so a rejected path leaves no trace.
                if let Some(name) = name {
                    match previous {
                        Some(span) => caps.insert(name.clone(), span),
                        None => caps.remove(name),
                    };
                }
                false
            });
            if !matched {
                *captures = saved;
            }
            matched
        }

        Node::Sequence(items) => match_sequence(items, subject, pos, captures, k),

        Node::Alternation(branches) => {
            for branch in branches {
                let saved = captures.clone();
                if match_node(branch, subject, pos, captures, k) {
                    return true;
                }
                *captures = saved;
            }
            false
        }

        Node::Repeat {
            node,
            min,
            max,
            greedy,
        } => match_repeat(node, *min, *max, *greedy, subject, pos, 0, captures, k),
    }
}

/// Matches a concatenation, threading the continuation through each element.
fn match_sequence(
    items: &[Node],
    subject: &[u8],
    pos: usize,
    captures: &mut Captures,
    k: Cont<'_>,
) -> bool {
    let Some((head, tail)) = items.split_first() else {
        return k(pos, captures);
    };
    match_node(head, subject, pos, captures, &mut |next, caps| {
        match_sequence(tail, subject, next, caps, k)
    })
}

/// Matches `node` repeatedly, honouring `min`, `max` and greediness.
///
/// A greedy quantifier tries one more repetition before trying to stop; a lazy
/// one tries to stop first. Both back off to the other order when the
/// continuation rejects, which is what lets `(?P<cell_1>.{8,12})GGG` shorten
/// the range until the literal fits.
#[allow(clippy::too_many_arguments)]
fn match_repeat(
    node: &Node,
    min: usize,
    max: Option<usize>,
    greedy: bool,
    subject: &[u8],
    pos: usize,
    count: usize,
    captures: &mut Captures,
    k: Cont<'_>,
) -> bool {
    let may_continue = max.is_none_or(|m| count < m);
    let may_stop = count >= min;

    let take_one = |caps: &mut Captures, k: Cont<'_>| -> bool {
        if !may_continue {
            return false;
        }
        match_node(node, subject, pos, caps, &mut |next, c| {
            // A zero-width element would otherwise repeat forever.
            next != pos && match_repeat(node, min, max, greedy, subject, next, count + 1, c, k)
        })
    };

    if greedy {
        if take_one(captures, k) {
            return true;
        }
        may_stop && k(pos, captures)
    } else {
        if may_stop && k(pos, captures) {
            return true;
        }
        take_one(captures, k)
    }
}

/// The candidate consumed-lengths for a fuzzy single-character atom, fewest
/// edits first — the order Python's engine tries them in.
fn fuzzy_lengths(kind: FuzzyKind, n: u32) -> &'static [usize] {
    if n == 0 {
        return &[1];
    }
    match kind {
        FuzzyKind::Substitution => &[1],
        FuzzyKind::Deletion => &[1, 0],
        FuzzyKind::Insertion => &[1, 2],
        FuzzyKind::Edit => &[1, 0, 2],
    }
}

/// Whether a fuzzy single-character atom admits the consumed slice.
fn fuzzy_admits(byte: u8, kind: FuzzyKind, n: u32, consumed: &[u8]) -> bool {
    if n == 0 {
        return consumed == [byte];
    }
    match (kind, consumed.len()) {
        // A substitution admits any single character (the exact one
        // included); a deletion admits consuming nothing. `Edit` allows both,
        // so the two share an arm.
        (FuzzyKind::Substitution | FuzzyKind::Edit, 1)
        | (FuzzyKind::Deletion | FuzzyKind::Edit, 0) => true,
        // For the insertion-only kind, one character is admissible only when
        // it is the exact one.
        (_, 1) => consumed[0] == byte,
        // An insertion: the character plus one extra, on either side.
        (FuzzyKind::Insertion | FuzzyKind::Edit, 2) => consumed[0] == byte || consumed[1] == byte,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(pattern: &str, subject: &str) -> Option<Match> {
        Pattern::parse(pattern)
            .expect("pattern parses")
            .match_at_start(subject)
    }

    #[test]
    fn literal_and_wildcard() {
        assert!(m("ACGT", "ACGTTTT").is_some());
        assert!(m("ACGT", "ACGA").is_none());
        assert!(m("A.GT", "ACGT").is_some());
        assert!(m(".{4}", "ACGT").is_some());
        assert!(m(".{5}", "ACGT").is_none());
    }

    #[test]
    fn anchors_behave_like_pythons_match() {
        // `.match` anchors at the start, so `^` is redundant but harmless.
        assert!(m("^ACGT", "ACGTTT").is_some());
        assert!(m("ACGT$", "ACGTTT").is_none());
        assert!(m("ACGT$", "ACGT").is_some());
        assert!(m(".*ACGT$", "TTTACGT").is_some());
    }

    #[test]
    fn named_groups_capture_spans() {
        let p = Pattern::parse("^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})").unwrap();
        let subject = "ACGTTTTGG";
        let hit = p.match_at_start(subject).unwrap();
        assert_eq!(hit.span("umi_1"), Some((0, 3)));
        assert_eq!(hit.span("umi_2"), Some((7, 9)));
        assert_eq!(hit.group("umi_1", subject), Some("ACG"));
        assert_eq!(hit.group("umi_2", subject), Some("GG"));
        assert_eq!(p.group_names(), ["umi_1", "umi_2"]);
    }

    #[test]
    fn three_prime_pattern_with_leading_wildcard() {
        let p = Pattern::parse(".*(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})$").unwrap();
        let subject = "AAAAACGTTTTGG";
        let hit = p.match_at_start(subject).unwrap();
        assert_eq!(hit.group("umi_1", subject), Some("ACG"));
        assert_eq!(hit.group("umi_2", subject), Some("GG"));
    }

    #[test]
    fn bounded_repeat_ranges() {
        assert!(m("(?P<cell_1>.{8,12})", "ACGTACGT").is_some());
        assert!(m("(?P<cell_1>.{8,12})", "ACGTACG").is_none());
        let p = Pattern::parse("(?P<cell_1>.{8,12})GGG").unwrap();
        // The greedy range must back off so the literal can match.
        let subject = "ACGTACGTAGGG";
        let hit = p.match_at_start(subject).unwrap();
        assert_eq!(hit.group("cell_1", subject), Some("ACGTACGTA"));
    }

    #[test]
    fn star_plus_and_optional() {
        assert!(m("A*", "").is_some());
        assert!(m("A+", "").is_none());
        assert!(m("A+", "AAA").is_some());
        assert!(m("AB?C", "AC").is_some());
        assert!(m("AB?C", "ABC").is_some());
    }

    #[test]
    fn character_classes() {
        assert!(m("[ACGT]{4}", "ACGT").is_some());
        assert!(m("[ACGT]{4}", "ACGN").is_none());
        assert!(m("[^N]{4}", "ACGT").is_some());
        assert!(m("[^N]{4}", "ACGN").is_none());
        assert!(m("[A-C]{2}", "AB").is_some());
        assert!(m(r"\d{2}", "42").is_some());
        assert!(m(r"\d{2}", "4A").is_none());
    }

    #[test]
    fn alternation() {
        assert!(m("ACGT|TTTT", "TTTT").is_some());
        assert!(m("ACGT|TTTT", "GGGG").is_none());
        let p = Pattern::parse("(?P<umi_1>ACGT|TTTT)").unwrap();
        assert_eq!(
            p.match_at_start("TTTT").unwrap().group("umi_1", "TTTT"),
            Some("TTTT")
        );
    }

    /// The semantics confirmed against the live `PyPI` `regex` module; see the
    /// module docs for the exact probes.
    #[test]
    fn fuzzy_substitution_matches_any_single_character() {
        assert!(m("(A{s<=1})X", "AX").is_some());
        assert!(m("(A{s<=1})X", "BX").is_some());
        assert!(m("(A{s<=1})X", "X").is_none());
        assert!(m("(A{s<=1})X", "ABX").is_none());
    }

    #[test]
    fn fuzzy_edit_allows_deletion_and_insertion() {
        assert!(m("(A{e<=1})X", "AX").is_some());
        assert!(m("(A{e<=1})X", "BX").is_some());
        assert!(m("(A{e<=1})X", "X").is_some());
        assert!(m("(A{e<=1})X", "ABX").is_some());
        assert!(m("(A{e<=1})X", "BCX").is_none());
    }

    #[test]
    fn fuzzy_insertion_and_deletion() {
        assert!(m("(A{i<=1})X", "AX").is_some());
        assert!(m("(A{i<=1})X", "AAX").is_some());
        assert!(m("(A{i<=1})X", "BX").is_none());
        assert!(m("(A{i<=1})X", "X").is_none());

        assert!(m("(A{d<=1})X", "AX").is_some());
        assert!(m("(A{d<=1})X", "X").is_some());
        assert!(m("(A{d<=1})X", "BX").is_none());
    }

    /// The inDrop pattern from the package's own test suite. Its `{s<=2}`
    /// binds to the final `T` only, so substituting the last base of the
    /// linker still matches while substituting the first does not.
    #[test]
    fn indrop_pattern_fuzzy_scope_matches_python() {
        let linker = "GAGTGATTGCTTGTGACGCCTT";
        let pattern = format!(
            "(?P<cell_1>.{{8,12}})(?P<discard_2>{linker}{{s<=2}})\
             (?P<cell_3>.{{8}})(?P<umi_1>.{{6}})T{{3}}.*"
        );
        let p = Pattern::parse(&pattern).unwrap();

        let build = |link: &str| format!("ACGTACGTAC{link}CCCCCCCCTTTTTTTTTAAAA");

        assert!(p.match_at_start(&build(linker)).is_some(), "exact linker");

        let mut last_sub = linker.to_string();
        last_sub.pop();
        last_sub.push('A');
        assert!(
            p.match_at_start(&build(&last_sub)).is_some(),
            "last base substituted"
        );

        let first_sub = format!("A{}", &linker[1..]);
        assert!(
            p.match_at_start(&build(&first_sub)).is_none(),
            "first base substituted must not match: {{s<=2}} binds to the last T only"
        );
    }

    #[test]
    fn fuzzy_on_a_group_is_rejected() {
        let err = Pattern::parse("(ACGT){s<=2}").unwrap_err();
        assert!(matches!(err, PortError::Format(ref msg)
            if msg.contains("single literal character")));
    }

    #[test]
    fn malformed_patterns_are_format_errors() {
        for bad in ["(ACGT", "ACGT)", "[ACGT", "A{", "A{,}", r"ACGT\"] {
            assert!(
                Pattern::parse(bad).is_err(),
                "expected {bad:?} to be rejected"
            );
        }
    }

    #[test]
    fn group_names_are_reported_in_declaration_order() {
        let p = Pattern::parse("(?P<cell_1>.{6})(?P<umi_1>.{10})(?P<discard_1>.*)").unwrap();
        assert_eq!(p.group_names(), ["cell_1", "umi_1", "discard_1"]);
    }

    #[test]
    fn empty_subject_only_matches_an_empty_pattern() {
        assert!(m("", "").is_some());
        assert!(m("A", "").is_none());
        assert!(m(".*", "").is_some());
    }
}
