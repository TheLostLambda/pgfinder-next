// FIXME: Order functions in the same way as the EBNF

use std::collections::HashMap;
use std::fmt;

use miette::{Diagnostic, LabeledSpan, SourceSpan};
use nom::combinator::cut;
use nom::{
    branch::alt,
    character::complete::{char, one_of, satisfy, u32},
    combinator::{all_consuming, complete, consumed, map, map_res, not, opt, recognize, success},
    error::{ErrorKind, FromExternalError, ParseError},
    multi::many1,
    sequence::{delimited, pair, preceded},
    Finish, IResult, Parser,
};
use thiserror::Error;

use super::{
    chemical_database::ChemicalDatabase, ChemicalComposition, ChemicalLookupError, Count, Element,
    MassNumber, OffsetKind, Particle,
};

// FIXME: Check this over, just jotting down ideas — names are awful...
// FIXME: Move all of this error handling and context wrapping to another crate to be shared
#[derive(Debug, Clone, Eq, PartialEq, Error)]
#[error("{error}")]
pub struct CompositionError {
    full_input: String,
    labels: Vec<LabeledSpan>,
    error: ErrorNode,
}

// FIXME: Bleh, naming!
#[derive(Debug, Clone, Eq, PartialEq, Error)]
pub enum ErrorNode {
    #[error("{kind}")]
    Node {
        kind: CompositionErrorKind,
        #[source]
        source: Option<Box<CompositionError>>,
    },
    #[error("attempted {} parse branches unsuccessfully", .0.len())]
    Branch(Vec<CompositionError>),
}

impl Diagnostic for CompositionError {
    fn source_code(&self) -> Option<&dyn miette::SourceCode> {
        Some(&self.full_input)
    }

    fn help<'a>(&'a self) -> Option<Box<dyn fmt::Display + 'a>> {
        if let ErrorNode::Node { kind, .. } = &self.error {
            kind.help()
        } else {
            None
        }
    }

    fn labels(&self) -> Option<Box<dyn Iterator<Item = miette::LabeledSpan> + '_>> {
        dbg!(&self.labels);
        Some(Box::new(self.labels.iter().cloned()))
    }

    fn related<'a>(&'a self) -> Option<Box<dyn Iterator<Item = &'a dyn Diagnostic> + 'a>> {
        if let ErrorNode::Branch(related) = &self.error {
            Some(Box::new(related.iter().map(|e| e as &dyn Diagnostic)))
        } else {
            None
        }
    }

    fn diagnostic_source(&self) -> Option<&dyn Diagnostic> {
        // FIXME: Better way to do this?
        if let ErrorNode::Node { source, .. } = &self.error {
            source.as_ref().map(|e| &**e as &dyn Diagnostic)
        } else {
            None
        }
    }
}

impl CompositionError {
    // FIXME: Build the whole labelled span here?
    // FIXME: Is mutable best here?
    fn bubble_labels(&mut self) {
        if self.labels.is_empty() {
            match &mut self.error {
                ErrorNode::Node {
                    source: Some(child),
                    ..
                } => {
                    child.bubble_labels();
                    self.labels = child.labels.drain(..).collect();
                }
                ErrorNode::Branch(alternatives) => {
                    let new_labels = alternatives.iter_mut().flat_map(|child| {
                        child.bubble_labels();
                        child.labels.drain(..)
                    });
                    self.labels = Self::merge_labels(new_labels);
                }
                ErrorNode::Node { .. } => (),
            }
        }
    }

    // FIXME: Where in the world does this belong...
    // FIXME: This eventually needs testing — showing that labels with different spans *don't* get merged
    fn merge_labels(labels: impl Iterator<Item = LabeledSpan>) -> Vec<LabeledSpan> {
        let mut span_map: HashMap<SourceSpan, Vec<String>> = HashMap::new();
        for labeled_span in labels {
            let span = labeled_span.inner();
            // FIXME: Gross with the clone() and to_owned() in here...
            let label = labeled_span.label().unwrap().to_owned();
            span_map
                .entry(*span)
                .and_modify(|l| l.push(label.clone()))
                .or_insert_with(|| vec![label]);
        }
        span_map
            .into_iter()
            .map(|(span, labels)| {
                let label = labels.join(" or ");
                LabeledSpan::new_with_span(Some(label), span)
            })
            .collect()
    }
}

// FIXME: Abstract this out into it's own module (or nom-supreme?)
// FIXME: Check that field ordering everywhere matches this!
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct CompositionParseError<'a> {
    input: &'a str,
    length: usize,
    reported: bool,
    fatal: bool,
    kind: CompositionErrorKind,
    alternatives: Vec<CompositionParseError<'a>>,
    source: Option<Box<CompositionParseError<'a>>>,
}

trait LabelledError {
    fn label(&self) -> Option<&'static str> {
        None
    }
}

// FIXME: Oh lord, what a mess...
impl<'a> CompositionParseError<'a> {
    // FIXME: Add a version that doesn't keep the source?
    // FIXME: Instead of all of this rebuilding with .., just take `mut self` like miette's GraphicalReportHandler
    fn kind(self, kind: CompositionErrorKind) -> Self {
        // FIXME: This erases any non-final kinds... Is that okay?
        let source = if self.reported {
            Some(Box::new(self.clone()))
        } else {
            None
        };
        Self {
            kind,
            reported: true,
            source,
            // FIXME: Boy, we really need a default constructor...
            alternatives: Vec::new(),
            ..self
        }
    }

    fn map_kind(self, from: CompositionErrorKind, to: CompositionErrorKind) -> Self {
        if self.kind == from {
            self.kind(to)
        } else {
            self
        }
    }

    fn fatal(self, fatal: bool) -> Self {
        Self { fatal, ..self }
    }

    // FIXME: Either get rid of these, or write setters for all of the properties...
    fn length(self, length: usize) -> Self {
        Self { length, ..self }
    }

    fn input(self, input: &'a str) -> Self {
        Self { input, ..self }
    }

    // FIXME: God help with the naming... Arguments too...
    fn into_final_error(self, full_input: &str) -> CompositionError {
        let source = self
            .source
            .clone()
            .map(|e| Box::new(e.into_final_error(full_input)));
        let alternatives: Vec<_> = self
            .alternatives
            .iter()
            .cloned()
            .map(|e| e.into_final_error(full_input))
            .collect();
        // NOTE: The additional space is added so that Diagnostic labels can point to the end of an input
        let input = format!("{full_input} ");
        let span = span_from_input(full_input, self.input, self.length);
        let label = self.kind.label().map(str::to_string);
        let labels = label
            .into_iter()
            .map(|l| LabeledSpan::new_with_span(Some(l), span))
            .collect();
        // FIXME: WET CODE!!!
        if alternatives.is_empty() {
            CompositionError {
                full_input: input,
                labels,
                error: ErrorNode::Node {
                    kind: self.kind,
                    source,
                },
            }
        } else {
            // FIXME: Yuckie!
            let mut other_self = self;
            other_self.alternatives.clear();
            CompositionError {
                full_input: input,
                labels: Vec::new(),
                error: ErrorNode::Branch(
                    [vec![other_self.into_final_error(full_input)], alternatives].concat(),
                ),
            }
        }
    }
}

// FIXME: OMG refactor... Also, can I do better than stealing here?
fn span_from_input(full_input: &str, input: &str, length: usize) -> SourceSpan {
    let base_addr = full_input.as_ptr() as usize;
    let substr_addr = input.as_ptr() as usize;
    // FIXME: Keep this?
    assert!(
        substr_addr >= base_addr,
        "tried to get the span of a non-substring!"
    );
    let start = substr_addr - base_addr;
    let end = start + length;
    SourceSpan::from(start..end)
}

pub fn final_parser<'a, O, P>(parser: P) -> impl FnMut(&'a str) -> Result<O, CompositionError>
where
    P: Parser<&'a str, O, CompositionParseError<'a>>,
{
    // FIXME: I can't inline this because of some borrow-checker closure witchcraft...
    let mut parser = all_consuming(complete(parser));
    move |input| {
        parser.parse(input).finish().map(|(_, c)| c).map_err(|e| {
            // FIXME: Ew
            let mut error = e.into_final_error(input);
            error.bubble_labels();
            error
        })
    }
}

pub fn map_res_span<'a, O1, O2, E2, F, G>(parser: F, f: G) -> impl FnMut(&'a str) -> ParseResult<O2>
where
    O1: Clone,
    F: Copy + Parser<&'a str, O1, CompositionParseError<'a>>,
    G: Copy + FnMut(O1) -> Result<O2, E2>,
    CompositionParseError<'a>: FromExternalError<&'a str, E2>,
{
    move |input| {
        let i = input;
        let (input, (consumed, o1)) = consumed(parser)(input)?;
        report_err(map_res(success(o1), f), |e| {
            e.input(i).length(consumed.len())
        })(input)
    }
}

// FIXME: Check if I'm being consistent about using `impl` or generics... I think I should avoid any generics I don't
// use, as long as this is just a sort of "internal" library
// FIXME: See if this signature can be simplified (elide lifetimes?)
// FIXME: Do I really need a closure? Or just a kind?
fn report_err<'a, O, P, F>(mut parser: P, f: F) -> impl FnMut(&'a str) -> ParseResult<O>
where
    P: Parser<&'a str, O, CompositionParseError<'a>>,
    F: Copy + FnOnce(CompositionParseError<'a>) -> CompositionParseError<'a>,
{
    move |i| {
        parser.parse(i).map_err(|e| match e.map(f) {
            nom::Err::Error(e) if e.fatal => nom::Err::Failure(e),
            rest => rest,
        })
    }
}

impl<'a> ParseError<&'a str> for CompositionParseError<'a> {
    fn from_error_kind(input: &'a str, kind: nom::error::ErrorKind) -> Self {
        // FIXME: Trying implementing Default and doing ..Self::default()
        Self {
            input,
            length: 0,
            kind: kind.into(),
            reported: false,
            fatal: false,
            alternatives: Vec::new(),
            source: None,
        }
    }

    fn append(_input: &str, _kind: nom::error::ErrorKind, other: Self) -> Self {
        other
    }

    // FIXME: Probably delete this!
    // FIXME: Just kidding... This probably needs to the part of the trait that the sub ErrorKind implements, along
    // with the labels, I suppose. And the from Nom error, and the from external?
    // FIXME: But, can I just replace this function with that ErrorKind -> MyError from impl?
    fn from_char(input: &'a str, c: char) -> Self {
        Self {
            input,
            length: 0,
            kind: CompositionErrorKind::Expected(c),
            reported: false,
            fatal: false,
            alternatives: Vec::new(),
            source: None,
        }
    }

    fn or(self, other: Self) -> Self {
        // FIXME: Oh god...
        let (new_self, alternative) = match (self.reported, other.reported) {
            (true, true) => (self.clone(), vec![other.clone()]),
            (true, false) => (self.clone(), Vec::new()),
            _ => (other.clone(), Vec::new()),
        };
        let alternatives = [self.alternatives, alternative].concat();
        Self {
            alternatives,
            ..new_self
        }
    }
}

// FIXME: Keep this in this file as well, since it's specific to ChemicalLookupError
// FIXME: Maybe merge with the LabeledError trait?
impl<'a> FromExternalError<&'a str, ChemicalLookupError> for CompositionParseError<'a> {
    fn from_external_error(input: &'a str, kind: ErrorKind, e: ChemicalLookupError) -> Self {
        Self::from_error_kind(input, kind)
            .kind(CompositionErrorKind::LookupError(e))
            .fatal(true)
    }
}

// FIXME: API guidelines, check word ordering
// FIXME: Keep this in this file when moving all of the other parser error stuff elsewhere!
// FIXME: MAKE THIS PRIVATE!
// FIXME: Also order these according to the EBNF
#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
pub enum CompositionErrorKind {
    #[error("expected a lowercase ASCII letter")]
    ExpectedLowercase,
    #[error("expected an uppercase ASCII letter")]
    ExpectedUppercase,
    #[diagnostic(help("a 0 value doesn't make sense here, if you've mistakenly included a leading zero, like NH02, try just NH2 instead"))]
    #[error("counts cannot start with 0")]
    ExpectedNoLeadingZero,
    // FIXME: Can I get rid of this?
    // FIXME: Move this out of the user-facing error type and make it something internal
    // FIXME: Actually, using an `expect` combinator, I might be able to actually get rid of this!
    #[diagnostic(help("this is an internal error that you shouldn't ever see! If you have gotten this error, then please report it as a bug!"))]
    #[error("expected '{0}'")]
    Expected(char),
    #[error("expected a chemical formula (optionally followed by a '+' or '-' and a particle offset), or a standalone particle offset")]
    ExpectedChemicalComposition,
    #[error(
        "expected an element (like Au) or an isotope (like [15N]) optionally followed by a number"
    )]
    ExpectedAtomicOffset,
    #[error("expected a particle (like p or e), optionally preceded by a number")]
    ExpectedParticleOffset,
    #[diagnostic(help("you've probably forgotten to close an earlier '[' bracket"))]
    #[error("expected ']' to close isotope brackets")]
    ExpectedIsotopeEnd,
    #[error("expected '[' to open isotope brackets")]
    ExpectedIsotopeStart,
    #[error("expected an isotopic mass number")]
    ExpectedMassNumber,
    #[error("expected an element symbol")]
    ExpectedElementSymbol,
    #[error("expected a particle symbol")]
    ExpectedParticleSymbol,
    #[diagnostic(help("double-check for typos, or add a new entry to the chemical database"))]
    #[error(transparent)]
    LookupError(ChemicalLookupError),
    #[diagnostic(help("this is an internal error that you shouldn't ever see! If you have gotten this error, then please report it as a bug!"))]
    #[error("internal `nom` error: {0:?}")]
    Nom(ErrorKind),
    #[diagnostic(help(
        "check the unparsed region for errors, or remove it from the rest of the composition"
    ))]
    #[error("could not interpret the full input as a valid chemical composition")]
    IncompleteParse,
}

impl From<ErrorKind> for CompositionErrorKind {
    fn from(value: ErrorKind) -> Self {
        match value {
            ErrorKind::Eof => Self::IncompleteParse,
            kind => Self::Nom(kind),
        }
    }
}

impl LabelledError for CompositionErrorKind {
    fn label(&self) -> Option<&'static str> {
        match self {
            Self::LookupError(ChemicalLookupError::Element(_)) => Some("element not found"),
            Self::LookupError(ChemicalLookupError::Isotope(_, _, _, _)) => {
                Some("isotope not found")
            }
            Self::LookupError(ChemicalLookupError::Particle(_)) => Some("particle not found"),
            Self::ExpectedUppercase => Some("expected uppercase"),
            Self::ExpectedLowercase => Some("expected lowercase"),
            Self::ExpectedIsotopeStart => Some("'['"),
            Self::ExpectedIsotopeEnd => Some("expected ']'"),
            Self::ExpectedMassNumber => Some("expected a number"),
            Self::ExpectedNoLeadingZero => Some("expected non-zero"),
            Self::IncompleteParse => Some("input was valid up until this point"),
            Self::Nom(_) => Some("the region that triggered this bug!"),
            _ => None,
        }
    }
}

// In the mod.rs, in the ChemicalComposition::new(), convert this ParserError to something with a pretty span for
// miette! Here all we're concerned about is the offending string and the error it produces.
// That code in mod.rs is also where the &str can be gotten rid of so that there isn't a lifetime in the error type
// NOTE: I can use the `consumed` combinator to get a nice &str source for things like element lookup errors

type ParseResult<'a, O> = IResult<&'a str, O, CompositionParseError<'a>>;

/// Element = uppercase , [ lowercase ] ;
fn element_symbol(i: &str) -> ParseResult<&str> {
    report_err(recognize(pair(uppercase, opt(lowercase))), |e| {
        e.kind(CompositionErrorKind::ExpectedElementSymbol)
    })(i)
}

/// Particle = lowercase ;
fn particle_symbol(i: &str) -> ParseResult<&str> {
    report_err(recognize(lowercase), |e| {
        e.kind(CompositionErrorKind::ExpectedParticleSymbol)
    })(i)
}

/// Isotope = "[" , Count , Element , "]" ;
fn isotope_expr(i: &str) -> ParseResult<(MassNumber, &str)> {
    report_err(
        delimited(char('['), cut(pair(count, element_symbol)), cut(char(']'))),
        |e| match e.kind {
            // FIXME: Replace match with map_kind?
            CompositionErrorKind::Expected('[') => {
                e.kind(CompositionErrorKind::ExpectedIsotopeStart)
            }
            CompositionErrorKind::Expected(']') => e.kind(CompositionErrorKind::ExpectedIsotopeEnd),
            CompositionErrorKind::Nom(ErrorKind::Digit) => {
                e.kind(CompositionErrorKind::ExpectedMassNumber)
            }
            _ => e,
        },
    )(i)
}

/// Element = uppercase , [ lowercase ] ;
fn element<'a>(db: &'a ChemicalDatabase) -> impl FnMut(&'a str) -> ParseResult<Element> {
    map_res_span(element_symbol, |symbol| Element::new(db, symbol))
}

/// Isotope = "[" , Count , Element , "]" ;
fn isotope<'a>(db: &'a ChemicalDatabase) -> impl FnMut(&'a str) -> ParseResult<Element> {
    map_res_span(isotope_expr, |(mass_number, symbol)| {
        Element::new_isotope(db, symbol, mass_number)
    })
}

/// Particle = lowercase ;
fn particle<'a>(db: &'a ChemicalDatabase) -> impl FnMut(&'a str) -> ParseResult<Particle> {
    map_res_span(particle_symbol, |symbol| Particle::new(db, symbol))
}

/// Offset Kind = "+" | "-" ;
fn offset_kind(i: &str) -> ParseResult<OffsetKind> {
    map(one_of("+-"), |c| match c {
        '+' => OffsetKind::Add,
        '-' => OffsetKind::Remove,
        _ => unreachable!(),
    })(i)
}

/// Count = digit - "0" , { digit } ;
fn count(i: &str) -> ParseResult<Count> {
    report_err(preceded(cut(not(char('0'))), u32), |e| {
        e.map_kind(
            CompositionErrorKind::Nom(ErrorKind::Not),
            CompositionErrorKind::ExpectedNoLeadingZero,
        )
    })(i)
}

/// Particle Offset = [ Count ] , Particle ;
fn particle_offset<'a>(
    db: &'a ChemicalDatabase,
) -> impl FnMut(&'a str) -> ParseResult<(Count, Particle)> {
    let optional_count = opt(count).map(|o| o.unwrap_or(1));
    report_err(pair(optional_count, particle(db)), |e| {
        e.kind(CompositionErrorKind::ExpectedParticleOffset)
    })
}

/// Atomic Offset = ( Element | Isotope ) , [ Count ] ;
fn atomic_offset<'a>(
    db: &'a ChemicalDatabase,
) -> impl FnMut(&'a str) -> ParseResult<(Element, Count)> {
    let optional_count = opt(count).map(|o| o.unwrap_or(1));
    report_err(pair(alt((element(db), isotope(db))), optional_count), |e| {
        e.kind(CompositionErrorKind::ExpectedAtomicOffset)
    })
}

/// Chemical Composition
///   = { Atomic Offset }- , [ Offset Kind , Particle Offset ]
///   | Particle Offset
///   ;
pub fn chemical_composition<'a>(
    db: &'a ChemicalDatabase,
) -> impl FnMut(&'a str) -> ParseResult<ChemicalComposition> {
    let atoms_and_particles = map(
        pair(
            many1(atomic_offset(db)),
            opt(pair(offset_kind, cut(particle_offset(db)))),
        ),
        |(chemical_formula, particle_offset)| {
            let particle_offset = particle_offset.map(|(k, (c, p))| (k, c, p));
            ChemicalComposition {
                chemical_formula,
                particle_offset,
            }
        },
    );
    let just_particles = map(particle_offset(db), |(c, p)| {
        let particle_offset = Some((OffsetKind::Add, c, p));
        ChemicalComposition {
            chemical_formula: Vec::new(),
            particle_offset,
        }
    });
    report_err(alt((atoms_and_particles, just_particles)), |e| {
        e.kind(CompositionErrorKind::ExpectedChemicalComposition)
    })
}

/// uppercase
///   = "A" | "B" | "C" | "D" | "E" | "F" | "G"
///   | "H" | "I" | "J" | "K" | "L" | "M" | "N"
///   | "O" | "P" | "Q" | "R" | "S" | "T" | "U"
///   | "V" | "W" | "X" | "Y" | "Z"
///   ;
fn uppercase(i: &str) -> ParseResult<char> {
    report_err(satisfy(|c| c.is_ascii_uppercase()), |e| {
        e.kind(CompositionErrorKind::ExpectedUppercase)
    })(i)
}

/// lowercase
///   = "a" | "b" | "c" | "d" | "e" | "f" | "g"
///   | "h" | "i" | "j" | "k" | "l" | "m" | "n"
///   | "o" | "p" | "q" | "r" | "s" | "t" | "u"
///   | "v" | "w" | "x" | "y" | "z"
///   ;
fn lowercase(i: &str) -> ParseResult<char> {
    // FIXME: Maybe make a let binding for all of the parsers that comes before the `report_err` — that way you can see
    // what the final nom code is, then you can see, seperately, where I start adding error information
    report_err(satisfy(|c| c.is_ascii_lowercase()), |e| {
        e.kind(CompositionErrorKind::ExpectedLowercase)
    })(i)
}

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use once_cell::sync::Lazy;

    use crate::testing_tools::assert_miette_snapshot;

    use super::*;

    static DB: Lazy<ChemicalDatabase> = Lazy::new(|| {
        ChemicalDatabase::from_kdl("chemistry.kdl", include_str!("../chemistry.kdl")).unwrap()
    });

    #[test]
    fn test_uppercase() {
        // Ensure the complete uppercase ASCII alphabet is present
        for c in 'A'..='Z' {
            assert_eq!(uppercase(&c.to_string()), Ok(("", c)));
        }
        // Ensure the complete lowercase ASCII alphabet is absent
        for c in 'a'..='z' {
            assert!(uppercase(&c.to_string()).is_err());
        }
        // Ensure only one character is parsed
        assert_eq!(uppercase("Hg"), Ok(("g", 'H')));
        assert_eq!(uppercase("HG"), Ok(("G", 'H')));
    }

    #[test]
    fn test_lowercase() {
        // Ensure the complete lowercase ASCII alphabet is present
        for c in 'a'..='z' {
            assert_eq!(lowercase(&c.to_string()), Ok(("", c)));
        }
        // Ensure the complete uppercase ASCII alphabet is absent
        for c in 'A'..='Z' {
            assert!(lowercase(&c.to_string()).is_err());
        }
        // Ensure only one character is parsed
        assert_eq!(lowercase("hg"), Ok(("g", 'h')));
        assert_eq!(lowercase("hG"), Ok(("G", 'h')));
    }

    #[test]
    fn test_count() {
        // Valid Counts
        assert_eq!(count("1"), Ok(("", 1)));
        assert_eq!(count("10"), Ok(("", 10)));
        assert_eq!(count("422"), Ok(("", 422)));
        assert_eq!(count("9999"), Ok(("", 9999)));
        // Invalid Counts
        assert!(count("0").is_err());
        assert!(count("01").is_err());
        assert!(count("00145").is_err());
        assert!(count("H").is_err());
        assert!(count("p").is_err());
        assert!(count("+H").is_err());
        assert!(count("[H]").is_err());
        // Multiple Counts
        assert_eq!(count("1OH"), Ok(("OH", 1)));
        assert_eq!(count("42HeH"), Ok(("HeH", 42)));
    }

    #[test]
    fn test_element_symbol() {
        // Valid Element Symbols
        assert_eq!(element_symbol("H"), Ok(("", "H")));
        assert_eq!(element_symbol("He"), Ok(("", "He")));
        // Invalid Element Symbols
        assert!(element_symbol("p").is_err());
        assert!(element_symbol("ep").is_err());
        assert!(element_symbol("1H").is_err());
        assert!(element_symbol("+H").is_err());
        assert!(element_symbol("[H]").is_err());
        // Multiple Element Symbols
        assert_eq!(element_symbol("OH"), Ok(("H", "O")));
        assert_eq!(element_symbol("HeH"), Ok(("H", "He")));
    }

    #[test]
    fn test_element() {
        let mut element = element(&DB);
        macro_rules! assert_element_name {
            ($input:literal, $output:literal, $name:literal) => {
                assert_eq!(
                    element($input).map(|(r, e)| (r, e.name)),
                    Ok(($output, $name.to_string()))
                );
            };
        }
        // Valid Elements
        assert_element_name!("H", "", "Hydrogen");
        assert_element_name!("He", "", "Helium");
        // Invalid Elements
        assert!(element("p").is_err());
        assert!(element("ep").is_err());
        assert!(element("1H").is_err());
        assert!(element("+H").is_err());
        assert!(element("[H]").is_err());
        // Non-Existent Elements
        assert!(element("X").is_err());
        assert!(element("To").is_err());
        // Multiple Elements
        assert_element_name!("OH", "H", "Oxygen");
        assert_element_name!("HeH", "H", "Helium");
    }

    #[test]
    fn test_isotope_expr() {
        // Valid Isotope Expressions
        assert_eq!(isotope_expr("[1H]"), Ok(("", (1, "H"))));
        assert_eq!(isotope_expr("[18O]"), Ok(("", (18, "O"))));
        assert_eq!(isotope_expr("[37Cl]"), Ok(("", (37, "Cl"))));
        // Invalid Isotope Expressions
        assert!(isotope_expr("H").is_err());
        assert!(isotope_expr("[H]").is_err());
        assert!(isotope_expr("[H2]").is_err());
        assert!(isotope_expr("[18OH]").is_err());
        assert!(isotope_expr("[18]").is_err());
        assert!(isotope_expr("[[18O]]").is_err());
        assert!(isotope_expr("[-18O]").is_err());
        assert!(isotope_expr("[+18O]").is_err());
        // Multiple Isotope Expressions
        assert_eq!(isotope_expr("[13C]O2"), Ok(("O2", (13, "C"))));
        assert_eq!(isotope_expr("[3He]H"), Ok(("H", (3, "He"))));
    }

    #[test]
    fn test_isotope() {
        let mut isotope = isotope(&DB);
        macro_rules! assert_isotope_name {
            ($input:literal, $output:literal, $name:literal) => {
                assert_eq!(
                    isotope($input)
                        .map(|(r, e)| (r, format!("{}-{}", e.name, e.mass_number.unwrap()))),
                    Ok(($output, $name.to_string()))
                );
            };
        }
        // Valid Isotopes
        assert_isotope_name!("[1H]", "", "Hydrogen-1");
        assert_isotope_name!("[18O]", "", "Oxygen-18");
        assert_isotope_name!("[37Cl]", "", "Chlorine-37");
        // Invalid Isotopes
        assert!(isotope("H").is_err());
        assert!(isotope("[H]").is_err());
        assert!(isotope("[H2]").is_err());
        assert!(isotope("[18OH]").is_err());
        assert!(isotope("[18]").is_err());
        assert!(isotope("[[18O]]").is_err());
        assert!(isotope("[-18O]").is_err());
        assert!(isotope("[+18O]").is_err());
        // Non-Existent Isotopes
        assert!(isotope("[42X]").is_err());
        assert!(isotope("[99To]").is_err());
        assert!(isotope("[15C]").is_err());
        assert!(isotope("[100Tc]").is_err());
        // Multiple Isotopes
        assert_isotope_name!("[13C]O2", "O2", "Carbon-13");
        assert_isotope_name!("[3He]H", "H", "Helium-3");
    }

    #[test]
    fn test_particle() {
        let mut particle = particle(&DB);
        macro_rules! assert_particle_name {
            ($input:literal, $output:literal, $name:literal) => {
                assert_eq!(
                    particle($input).map(|(r, e)| (r, e.name)),
                    Ok(($output, $name.to_string()))
                );
            };
        }
        // Valid Particles
        assert_particle_name!("p", "", "Proton");
        assert_particle_name!("e", "", "Electron");
        // Invalid Particles
        assert!(particle("P").is_err());
        assert!(particle("Ep").is_err());
        assert!(particle("1p").is_err());
        assert!(particle("+e").is_err());
        assert!(particle("[p]").is_err());
        // Non-Existent Particles
        assert!(particle("m").is_err());
        assert!(particle("g").is_err());
        // Multiple Particles
        assert_particle_name!("ep", "p", "Electron");
        assert_particle_name!("pe", "e", "Proton");
    }

    #[test]
    fn test_offset_kind() {
        // Valid Offset Kinds
        assert_eq!(offset_kind("+"), Ok(("", OffsetKind::Add)));
        assert_eq!(offset_kind("-"), Ok(("", OffsetKind::Remove)));
        // Invalid Offset Kinds
        assert!(offset_kind("p").is_err());
        assert!(offset_kind("H").is_err());
        assert!(offset_kind("1H").is_err());
        assert!(offset_kind("1+H").is_err());
        assert!(offset_kind("[H]").is_err());
        // Multiple Offset Kinds
        assert_eq!(offset_kind("+-"), Ok(("-", OffsetKind::Add)));
        assert_eq!(offset_kind("--"), Ok(("-", OffsetKind::Remove)));
    }

    #[test]
    fn test_particle_offset() {
        let mut particle_offset = particle_offset(&DB);
        macro_rules! assert_particle_offset {
            ($input:literal, $output:literal, $count:literal, $name:literal ) => {
                assert_eq!(
                    particle_offset($input).map(|(r, (c, p))| (r, (c, p.name))),
                    Ok(($output, ($count, $name.to_string())))
                );
            };
        }
        // Valid Particle Offsets
        assert_particle_offset!("p", "", 1, "Proton");
        assert_particle_offset!("e", "", 1, "Electron");

        assert_particle_offset!("5p", "", 5, "Proton");
        assert_particle_offset!("5e", "", 5, "Electron");

        assert_particle_offset!("100p", "", 100, "Proton");
        assert_particle_offset!("100e", "", 100, "Electron");
        // Invalid Particle Offsets
        assert!(particle_offset("-p").is_err());
        assert!(particle_offset("-e").is_err());
        assert!(particle_offset("+p").is_err());
        assert!(particle_offset("+e").is_err());
        assert!(particle_offset("P").is_err());
        assert!(particle_offset("Ep").is_err());
        assert!(particle_offset("1+p").is_err());
        assert!(particle_offset("+0p").is_err());
        assert!(particle_offset("-02e").is_err());
        assert!(particle_offset("+-e").is_err());
        assert!(particle_offset("+42").is_err());
        assert!(particle_offset("-[p]").is_err());
        // Non-Existent Particle Offsets
        assert!(particle_offset("+m").is_err());
        assert!(particle_offset("-3g").is_err());
        // Multiple Particle Offsets
        assert_particle_offset!("3ep", "p", 3, "Electron");
        assert_particle_offset!("pe", "e", 1, "Proton");
    }

    #[test]
    fn test_atomic_offset() {
        let mut atomic_offset = atomic_offset(&DB);
        macro_rules! assert_atomic_offset {
            ($input:literal, $output:literal, $name:literal, $count:literal) => {
                assert_eq!(
                    atomic_offset($input).map(|(r, (e, c))| {
                        let name = if let Some(a) = e.mass_number {
                            format!("{}-{a}", e.name)
                        } else {
                            e.name
                        };
                        (r, (name, c))
                    }),
                    Ok(($output, ($name.to_string(), $count)))
                );
            };
        }
        // Valid Atomic Offsets
        assert_atomic_offset!("H", "", "Hydrogen", 1);
        assert_atomic_offset!("He", "", "Helium", 1);
        assert_atomic_offset!("H2", "", "Hydrogen", 2);
        assert_atomic_offset!("C18", "", "Carbon", 18);
        assert_atomic_offset!("[1H]", "", "Hydrogen-1", 1);
        assert_atomic_offset!("[18O]", "", "Oxygen-18", 1);
        assert_atomic_offset!("[37Cl]", "", "Chlorine-37", 1);
        assert_atomic_offset!("[1H]2", "", "Hydrogen-1", 2);
        assert_atomic_offset!("[18O]3", "", "Oxygen-18", 3);
        assert_atomic_offset!("[37Cl]5", "", "Chlorine-37", 5);
        // Invalid Atomic Offsets
        assert!(atomic_offset("p").is_err());
        assert!(atomic_offset("-2p").is_err());
        assert!(atomic_offset("ep").is_err());
        assert!(atomic_offset("1H").is_err());
        assert!(atomic_offset("+H").is_err());
        assert!(atomic_offset("2[2H]").is_err());
        assert!(atomic_offset("[H]").is_err());
        assert!(atomic_offset("[H2]").is_err());
        assert!(atomic_offset("[18OH]").is_err());
        assert!(atomic_offset("[18]").is_err());
        assert!(atomic_offset("[[18O]]").is_err());
        assert!(atomic_offset("[-18O]").is_err());
        assert!(atomic_offset("[+18O]").is_err());
        // Non-Existent Atomic Offsets
        assert!(atomic_offset("X2").is_err());
        assert!(atomic_offset("To7").is_err());
        assert!(atomic_offset("[42X]").is_err());
        assert!(atomic_offset("[99To]2").is_err());
        assert!(atomic_offset("[15C]4").is_err());
        assert!(atomic_offset("[100Tc]8").is_err());
        // Multiple Atomic Offsets
        assert_atomic_offset!("OH", "H", "Oxygen", 1);
        assert_atomic_offset!("HeH", "H", "Helium", 1);
        assert_atomic_offset!("H2O", "O", "Hydrogen", 2);
        assert_atomic_offset!("CO2", "O2", "Carbon", 1);
        assert_atomic_offset!("[13C]6O2", "O2", "Carbon-13", 6);
        assert_atomic_offset!("[3He]1H", "H", "Helium-3", 1);
    }

    #[test]
    fn test_chemical_composition() {
        let mut chemical_composition = chemical_composition(&DB);
        macro_rules! check_composition_snapshot {
            ($input:literal, $output:literal) => {
                let (
                    rest,
                    ChemicalComposition {
                        chemical_formula,
                        particle_offset,
                    },
                ) = chemical_composition($input).unwrap();
                assert_eq!(rest, $output);
                let chemical_formula: Vec<_> = chemical_formula
                    .iter()
                    .map(|(e, c)| {
                        let name = if let Some(a) = e.mass_number {
                            format!("{}-{a}", e.name)
                        } else {
                            e.name.clone()
                        };
                        (name, c)
                    })
                    .collect();
                let particle_offset = particle_offset.map(|(k, c, p)| (k, c, p.name));
                assert_debug_snapshot!((chemical_formula, particle_offset));
            };
        }
        // Valid Chemical Compositions
        check_composition_snapshot!("H2O", "");
        check_composition_snapshot!("C11H12N2O2", "");
        check_composition_snapshot!("OH+e", "");
        check_composition_snapshot!("Na-e", "");
        check_composition_snapshot!("NH3+p", "");
        check_composition_snapshot!("Cr2O7+2e", "");
        check_composition_snapshot!("NH2+2p", "");
        check_composition_snapshot!("D2O", "");
        check_composition_snapshot!("[2H]2O", "");
        check_composition_snapshot!("[37Cl]5-2p", "");
        check_composition_snapshot!("NH2[99Tc]", "");
        // Invalid Chemical Compositions
        assert!(chemical_composition(" ").is_err());
        assert!(chemical_composition("-2p").is_err());
        assert!(chemical_composition("+H").is_err());
        assert!(chemical_composition("2[2H]").is_err());
        assert!(chemical_composition("[H+p]O").is_err());
        assert!(chemical_composition("NH2[100Tc]").is_err());
        // Multiple Chemical Compositions
        check_composition_snapshot!("[37Cl]5-2p10", "10");
        check_composition_snapshot!("[2H]2O*H2O", "*H2O");
        check_composition_snapshot!("NH2{100Tc", "{100Tc");
        check_composition_snapshot!("C11H12N2O2 H2O", " H2O");
    }

    #[test]
    fn test_composition_errors() {
        let mut chemical_composition = final_parser(chemical_composition(&DB));

        // Looking up non-existant isotopes, elements, and particles
        assert_miette_snapshot!(chemical_composition("NH2[100Tc]O4"));
        assert_miette_snapshot!(chemical_composition("NH2[99Tc]YhO4"));
        assert_miette_snapshot!(chemical_composition("NH2[99Tc]O4-8m+2p"));

        // Starting a composition without an element or isotope
        assert_miette_snapshot!(chemical_composition("+H2O"));
        assert_miette_snapshot!(chemical_composition("-H2O"));
        assert_miette_snapshot!(chemical_composition("]H2O"));

        // Check counts are non-zero (no leading zeroes either!)
        assert_miette_snapshot!(chemical_composition("C3H0N4"));
        assert_miette_snapshot!(chemical_composition("C3H06N4"));

        // Ensure that particles are lowercase
        assert_miette_snapshot!(chemical_composition("H2O+P"));

        // Ensure that isotope expressions are valid
        assert_miette_snapshot!(chemical_composition("[H2O"));
        assert_miette_snapshot!(chemical_composition("[0H2O"));
        assert_miette_snapshot!(chemical_composition("[10H2O"));
        assert_miette_snapshot!(chemical_composition("[10]H2O"));
        assert_miette_snapshot!(chemical_composition("[37Cl"));

        // Check labels at the end of an input
        assert_miette_snapshot!(chemical_composition("[37Cl]5-"));
        assert_miette_snapshot!(chemical_composition("[37Cl]5+10"));

        // Check for partially valid input
        assert_miette_snapshot!(chemical_composition("[37Cl]52p"));
        assert_miette_snapshot!(chemical_composition("NH2[99Tc]O,4-2e+3p"));
        assert_miette_snapshot!(chemical_composition("eH2O"));
    }
}
