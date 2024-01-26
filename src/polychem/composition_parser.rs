// FIXME: Starting with getting the "happy path" working, then adding in nice error reporting! Don't forget to test!
// FIXME: Make sure to update all of the `is_err()` tests to check that the error contains the correct, rich info...
// FIXME: Order functions in the same way as the EBNF

use std::fmt;
use std::{cmp::max, iter};

use miette::{Diagnostic, LabeledSpan, SourceSpan};
use nom::{
    branch::alt,
    character::complete::{char, one_of, satisfy, u32},
    combinator::{all_consuming, complete, consumed, map, map_res, not, opt, recognize, success},
    error::{ErrorKind, FromExternalError, ParseError},
    multi::{many0, many1},
    sequence::{delimited, pair, preceded, tuple},
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
#[error("{kind}")]
pub struct CompositionError {
    input: String,
    // FIXME: Think about making that label dynamic?
    span: SourceSpan,
    // #[transparent] or something?
    // FIXME: Hopefully `help` gets passed through this, otherwise it's a manual impl for me...
    kind: CompositionErrorKind,
    // Append these!
    // Should failed alt paths be collected here?
    #[source]
    source: Option<Box<CompositionError>>,
}

#[derive(Debug, Clone, Eq, PartialEq, Error)]
#[error(transparent)]
struct SubError<E>(E);

impl Diagnostic for CompositionError {
    fn source_code(&self) -> Option<&dyn miette::SourceCode> {
        Some(&self.input)
    }

    fn labels(&self) -> Option<Box<dyn Iterator<Item = miette::LabeledSpan> + '_>> {
        // FIXME: This should probably only print the label of the error the lowest down the stack that has a non-None
        // label — higher labels should be ignored? But float it up to the top?
        // FIXME: Unsafe indexing again, maybe drop the whole Vec<> again...
        let label = if let Some(l) = self.kind.label() {
            l
        } else {
            // FIXME: This needs to work for an arbitrary depth!
            self.source.as_ref()?.kind.label()?
        };
        Some(Box::new(iter::once(LabeledSpan::new_with_span(
            Some(label.to_string()),
            self.span,
        ))))
    }

    fn diagnostic_source(&self) -> Option<&dyn Diagnostic> {
        // Some(&self.kind)
        // TODO: Merge all of the sources here!
        // FIXME: Also be sure to implement for Error manually...
        // FIXME: Need to actually merge errors!
        // FIXME: Shouldn't be indexing, can panic!
        // self.sources.get(0).map(|e| &e.errors[0] as &dyn Diagnostic)
        self.source.as_ref().map(|s| s as &dyn Diagnostic)
    }

    fn help<'a>(&'a self) -> Option<Box<dyn fmt::Display + 'a>> {
        // FIXME: Need to also merge help strings!
        self.kind.help()
    }
}

impl Diagnostic for Box<CompositionError> {
    fn help<'a>(&'a self) -> Option<Box<dyn fmt::Display + 'a>> {
        self.kind.help()
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
    // FIXME: Add a field for a label string!
    kind: CompositionErrorKind,
    source: Option<Box<CompositionParseError<'a>>>,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Ord, PartialOrd)]
enum ErrorStatus {
    Unhandled,
    Reported,
    Fatal,
}

trait LabelledError {
    fn label(&self) -> Option<&'static str> {
        None
    }
}

// FIXME: Maybe make this into a combinator like `map_res`? Or maybe an extension trait...
// FIXME: Oh lord, what a mess...
// FIXME: also write append_error / kind
// FIXME: I might need to create my own version of `all_consuming` that's capable of reporting all of the `alt`
// branches that failed there! Otherwise I just get a pretty unhelpful Eof error
// FIXME: I should do that extension trait, so I can chain methods on nom::Err<CompositionParseError>
// FIXME: Not needed! The combinator and do the wrapping and unwrapping; impl CompositionParseError will do!
impl<'a> CompositionParseError<'a> {
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
            ..self
        }
    }

    fn fatal_if(self, condition: impl Fn(&Self) -> bool) -> Self {
        let fatal = condition(&self);
        Self { fatal, ..self }
    }

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
            .map(|e| Box::new(e.into_final_error(full_input)));
        CompositionError {
            // Extra space allows for labels at the end of an input
            input: format!("{full_input} "),
            span: span_from_input(full_input, self.input, self.length),
            kind: self.kind,
            source,
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
        parser
            .parse(input)
            .finish()
            .map(|(_, c)| c)
            .map_err(|e| e.into_final_error(input))
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
            source: None,
        }
    }

    fn append(_input: &str, _kind: nom::error::ErrorKind, other: Self) -> Self {
        // FIXME: What does this actually do?
        // This seems to just add stack context, but I don't think I actually care about which higher-level errors
        // accumulate, I think I'm only interest in the base cause and the other paths that also failed. If I want to
        // wrap an error at a higher level, I'll do that manually!
        other
    }

    fn from_char(input: &'a str, c: char) -> Self {
        Self {
            input,
            length: 0,
            kind: CompositionErrorKind::Expected(c),
            reported: false,
            fatal: false,
            source: None,
        }
    }

    // FIXME: I still need to filter out the rubbish (keep final errors) but I can just return one of the final ones
    fn or(self, other: Self) -> Self {
        let kind = if self.reported { self.kind } else { other.kind };
        // FIXME: This logic is probably wrong — merging locations? Asserting they are the same?
        let status = max(self.reported, other.reported);
        let length = max(self.length, other.length);
        Self {
            kind,
            length,
            reported: status,
            ..other
        }
    }
}

impl<'a> FromExternalError<&'a str, ChemicalLookupError> for CompositionParseError<'a> {
    fn from_external_error(input: &'a str, _kind: ErrorKind, e: ChemicalLookupError) -> Self {
        Self {
            input,
            length: 0,
            kind: CompositionErrorKind::LookupError(e),
            reported: true,
            fatal: true,
            source: None,
        }
    }
}

// FIXME: API guidelines, check word ordering
// FIXME: Keep this in this file when moving all of the other parser error stuff elsewhere!
// FIXME: Need to make all of these error messages start with lowercase and word to compose better! expected...
#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
enum CompositionErrorKind {
    // ... #[help] [#error] etc
    #[error("Should be lowercase, yo")]
    ExpectedLowercase,
    #[error("Should be uppercase, yo")]
    ExpectedUppercase,
    #[error("Should be a count (can't start with a zero!)")]
    ExpectedCount,
    #[error("Should be {0:?}, yo")]
    Expected(char),
    #[diagnostic(help("Have you tried..."))]
    #[error("Should be an element (Au) or and isotope [15N] optionally followed by a number")]
    ExpectedAtomicOffset,
    #[diagnostic(help("Have you tried being better?"))]
    #[error("Should a ± number and particle")]
    ExpectedParticleOffset,
    // ExpectedOffsetKind
    // etc...
    #[diagnostic(help("Uh, trying using something that actually exists, dumbo..."))]
    #[error(transparent)]
    LookupError(ChemicalLookupError),
    #[error("Nommy mommy wet itself: {0:?}")]
    Nom(ErrorKind),
    #[error("Unexpectedly ran into gibberish!")]
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
            // FIXME: Should I just implement LabelledError for ChemicalLookupError? I don't think so, since that error
            // isn't normally formatted using miette (it's not parsing anything, so what source do you point to?)
            // Here I've just used one string for all of the sub-types...
            Self::LookupError(ChemicalLookupError::Element(_)) => Some("element not found"),
            Self::LookupError(ChemicalLookupError::Isotope(_, _)) => Some("isotope not found"),
            Self::LookupError(ChemicalLookupError::Particle(_)) => Some("particle not found"),
            Self::ExpectedUppercase => Some("expected uppercase"),
            Self::ExpectedLowercase => Some("expected lowercase"),
            Self::ExpectedCount => Some("expected non-zero"),
            Self::IncompleteParse => Some("input was valid up until this point"),
            Self::Nom(_) => Some("here"),
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
    recognize(pair(uppercase, opt(lowercase)))(i)
}

fn particle_symbol(i: &str) -> ParseResult<&str> {
    recognize(lowercase)(i)
}

/// Isotope = "[" , Integer , Element , "]" ;
fn isotope_expr(i: &str) -> ParseResult<(MassNumber, &str)> {
    delimited(char('['), pair(count, element_symbol), char(']'))(i)
}

/// Element = uppercase , [ lowercase ] ;
fn element<'a>(db: &'a ChemicalDatabase) -> impl FnMut(&'a str) -> ParseResult<Element> {
    map_res_span(element_symbol, |symbol| Element::new(db, symbol))
}

/// Isotope = "[" , Integer , Element , "]" ;
fn isotope<'a>(db: &'a ChemicalDatabase) -> impl FnMut(&'a str) -> ParseResult<Element> {
    map_res_span(isotope_expr, |(mass_number, symbol)| {
        Element::new_isotope(db, symbol, mass_number)
    })
}

/// Particle = lowercase ;
fn particle<'a>(db: &'a ChemicalDatabase) -> impl FnMut(&'a str) -> ParseResult<Particle> {
    map_res_span(particle_symbol, |symbol| Particle::new(db, symbol))
}

/// Offset errors = "+" | "-" ;
fn offset_kind(i: &str) -> ParseResult<OffsetKind> {
    map(one_of("+-"), |c| match c {
        '+' => OffsetKind::Add,
        '-' => OffsetKind::Remove,
        _ => unreachable!(),
    })(i)
}

/// Count = digit - "0" , { digit } ;
fn count(i: &str) -> ParseResult<Count> {
    report_err(preceded(not(char('0')), u32), |e| {
        // FIXME: Maybe change `fatal_if` just to take a nom ErrorKind? Then it does the checking itself
        e.fatal_if(|e| e.kind == CompositionErrorKind::Nom(ErrorKind::Not))
            .kind(CompositionErrorKind::ExpectedCount)
    })(i)
}

/// Particle Offset = Offset Kind , [ Integer ] ,
///   Particle ;
fn particle_offset<'a>(
    db: &'a ChemicalDatabase,
) -> impl FnMut(&'a str) -> ParseResult<(OffsetKind, Count, Particle)> {
    let optional_count = opt(count).map(|o| o.unwrap_or(1));
    report_err(tuple((offset_kind, optional_count, particle(db))), |e| {
        e.fatal_if(|e| e.kind == CompositionErrorKind::ExpectedLowercase)
            .kind(CompositionErrorKind::ExpectedParticleOffset)
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

/// Chemical Composition = { Atomic Offset }- ,
///   { Particle Offset } ;
pub fn chemical_composition<'a>(
    db: &'a ChemicalDatabase,
) -> impl FnMut(&'a str) -> ParseResult<ChemicalComposition> {
    let parts = pair(many1(atomic_offset(db)), many0(particle_offset(db)));
    map(parts, |(chemical_formula, charged_particles)| {
        ChemicalComposition {
            chemical_formula,
            charged_particles,
        }
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
        e.kind(CompositionErrorKind::ExpectedUppercase).length(1)
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
        e.kind(CompositionErrorKind::ExpectedLowercase).length(1)
    })(i)
}

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use once_cell::sync::Lazy;

    use crate::polychem::testing_tools::assert_miette_snapshot;

    use super::*;

    static DB: Lazy<ChemicalDatabase> = Lazy::new(|| {
        ChemicalDatabase::from_kdl("chemistry.kdl", include_str!("chemistry.kdl")).unwrap()
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
        // FIXME: Add a check for a richer error message when parsing fails
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
        // FIXME: Add a check for a richer error message when parsing fails
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
            ($input:literal, $output:literal, $kind:expr, $count:literal, $name:literal ) => {
                assert_eq!(
                    particle_offset($input).map(|(r, (k, c, p))| (r, (k, c, p.name))),
                    Ok(($output, ($kind, $count, $name.to_string())))
                );
            };
        }
        // Valid Particle Offsets
        assert_particle_offset!("+p", "", OffsetKind::Add, 1, "Proton");
        assert_particle_offset!("-p", "", OffsetKind::Remove, 1, "Proton");
        assert_particle_offset!("+e", "", OffsetKind::Add, 1, "Electron");
        assert_particle_offset!("-e", "", OffsetKind::Remove, 1, "Electron");

        assert_particle_offset!("+5p", "", OffsetKind::Add, 5, "Proton");
        assert_particle_offset!("-5p", "", OffsetKind::Remove, 5, "Proton");
        assert_particle_offset!("+5e", "", OffsetKind::Add, 5, "Electron");
        assert_particle_offset!("-5e", "", OffsetKind::Remove, 5, "Electron");

        assert_particle_offset!("+100p", "", OffsetKind::Add, 100, "Proton");
        assert_particle_offset!("-100p", "", OffsetKind::Remove, 100, "Proton");
        assert_particle_offset!("+100e", "", OffsetKind::Add, 100, "Electron");
        assert_particle_offset!("-100e", "", OffsetKind::Remove, 100, "Electron");
        // Invalid Particle Offsets
        assert!(particle_offset("P").is_err());
        assert!(particle_offset("Ep").is_err());
        assert!(particle_offset("1p").is_err());
        assert!(particle_offset("1+p").is_err());
        assert!(particle_offset("p").is_err());
        assert!(particle_offset("+0p").is_err());
        assert!(particle_offset("-02e").is_err());
        assert!(particle_offset("+-e").is_err());
        assert!(particle_offset("+42").is_err());
        assert!(particle_offset("-[p]").is_err());
        // Non-Existent Particle Offsets
        assert!(particle_offset("+m").is_err());
        assert!(particle_offset("-3g").is_err());
        // Multiple Particle Offsets
        assert_particle_offset!("+3ep", "p", OffsetKind::Add, 3, "Electron");
        assert_particle_offset!("-pe", "e", OffsetKind::Remove, 1, "Proton");
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
                        charged_particles,
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
                let charged_particles: Vec<_> = charged_particles
                    .iter()
                    .map(|(k, c, p)| (k, c, &p.name))
                    .collect();
                assert_debug_snapshot!((chemical_formula, charged_particles));
            };
        }
        // Valid Chemical Compositions
        check_composition_snapshot!("H2O", "");
        check_composition_snapshot!("C11H12N2O2", "");
        check_composition_snapshot!("OH+e", "");
        check_composition_snapshot!("Na-e", "");
        check_composition_snapshot!("NH3+p", "");
        check_composition_snapshot!("Cr2O7+2e", "");
        check_composition_snapshot!("NH2+2p+e", "");
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
        check_composition_snapshot!("NH2[100Tc", "[100Tc");
        check_composition_snapshot!("C11H12N2O2 H2O", " H2O");
    }

    // FIXME: Dirty! Rename and refactor
    #[test]
    fn test_errors() -> miette::Result<()> {
        let mut chemical_composition = final_parser(chemical_composition(&DB));
        // FIXME: Don't forget to uncomment or delete!
        assert_miette_snapshot!(chemical_composition("NH2[100Tc]O4"));
        assert_miette_snapshot!(chemical_composition("NH2[99Tc]YhO4"));
        assert_miette_snapshot!(chemical_composition("NH2[99Tc]O4-8m+2p"));
        assert_miette_snapshot!(chemical_composition("xH2O"));
        assert_miette_snapshot!(chemical_composition("-H2O"));
        assert_miette_snapshot!(chemical_composition("NH2[99Tc]O,4-2e+3p"));
        // Check counts are non-zero (no leading zeroes either!)
        assert_miette_snapshot!(chemical_composition("C3H0N4"));
        assert_miette_snapshot!(chemical_composition("C3H06N4"));
        assert_miette_snapshot!(chemical_composition("H2O+P"));
        // Check labels at the end of an input
        assert_miette_snapshot!(chemical_composition("[37Cl]5-2p+10"));
        Ok(())
    }
}
