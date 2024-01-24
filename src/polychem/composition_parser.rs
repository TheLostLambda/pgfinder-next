// FIXME: Starting with getting the "happy path" working, then adding in nice error reporting! Don't forget to test!
// FIXME: Make sure to update all of the `is_err()` tests to check that the error contains the correct, rich info...
// FIXME: Order functions in the same way as the EBNF

use miette::{Diagnostic, SourceSpan};
use nom::{
    branch::alt,
    character::complete::{char, one_of, satisfy, u32},
    combinator::{map, map_res, not, opt, recognize},
    error::{ErrorKind, FromExternalError, ParseError},
    multi::{many0, many1},
    sequence::{delimited, pair, preceded, tuple},
    IResult, Parser,
};
use thiserror::Error;

use super::{
    chemical_database::ChemicalDatabase, ChemicalComposition, ChemicalLookupError, Count, Element,
    MassNumber, OffsetKind, Particle,
};

// FIXME: Check this over, just jotting down ideas — names are awful...
// FIXME: Move all of this error handling and context wrapping to another crate to be shared
#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
#[error("{kind}")]
pub struct CompositionError {
    #[source_code]
    input: String,
    // FIXME: Think about making that label dynamic?
    #[label("here")]
    span: SourceSpan,
    // #[transparent] or something?
    // FIXME: Hopefully `help` gets passed through this, otherwise it's a manual impl for me...
    #[source]
    #[diagnostic_source]
    kind: CompositionErrorKind,
    // Append these!
    #[related]
    related: Vec<CompositionErrorKind>,
}

// FIXME: Abstract this out into it's own module (or nom-supreme?)
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct CompositionParseError<'a> {
    input: &'a str,
    length: usize,
    failures: Vec<CompositionErrorKind>,
    // FIXME: Maybe don't even both with this... If I do, I might want to add the location info?
    related: Vec<CompositionErrorKind>,
}

// FIXME: Maybe make this into a combinator like `map_res`? Or maybe an extension trait...
// FIXME: Oh lord, what a mess...
// FIXME: also write append_error / kind
// FIXME: I might need to create my own version of `all_consuming` that's capable of reporting all of the `alt`
// branches that failed there! Otherwise I just get a pretty unhelpful Eof error
// FIXME: I should do that extension trait, so I can chain methods on nom::Err<CompositionParseError>
// FIXME: Not needed! The combinator and do the wrapping and unwrapping; impl CompositionParseError will do!
fn set_kind(
    error: nom::Err<CompositionParseError>,
    kind: CompositionErrorKind,
) -> nom::Err<CompositionParseError> {
    error.map(|e| CompositionParseError {
        failures: vec![kind],
        ..e
    })
}

impl<'a> ParseError<&'a str> for CompositionParseError<'a> {
    fn from_error_kind(input: &'a str, kind: nom::error::ErrorKind) -> Self {
        Self {
            input,
            length: 0,
            failures: vec![CompositionErrorKind::Nom(kind)],
            related: Vec::new(),
        }
    }

    fn append(_input: &str, _kind: nom::error::ErrorKind, other: Self) -> Self {
        // FIXME: What does this actually do?
        // This seems to just add stack context, but I don't think I actually care about which higher-level errors
        // accumulate, I think I'm only interest in the base cause and the other paths that also failed. If I want to
        // wrap an error at a higher level, I'll do that manually!
        dbg!(_input, _kind, &other);
        other
    }

    fn from_char(input: &'a str, c: char) -> Self {
        Self {
            input,
            length: 0,
            failures: vec![CompositionErrorKind::Expected(c)],
            related: Vec::new(),
        }
    }

    fn or(self, other: Self) -> Self {
        dbg!(&self, &other);
        let failures = [self.failures, other.failures].concat();
        // FIXME: This logic is probably wrong
        Self { failures, ..other }
    }
}

impl<'a> FromExternalError<&'a str, ChemicalLookupError> for CompositionParseError<'a> {
    fn from_external_error(input: &'a str, _kind: ErrorKind, e: ChemicalLookupError) -> Self {
        Self {
            input,
            length: 0,
            failures: vec![CompositionErrorKind::LookupError(e)],
            related: Vec::new(),
        }
    }
}

// FIXME: API guidelines, check word ordering
#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
enum CompositionErrorKind {
    // ... #[help] [#error] etc
    #[error("Should be lowercase, yo")]
    ExpectedLowercase,
    #[error("Should be uppercase, yo")]
    ExpectedUppercase,
    #[error("Should be {0:?}, yo")]
    Expected(char),
    // ExpectedOffsetKind
    // etc...
    #[error(transparent)]
    LookupError(ChemicalLookupError),
    #[error("Nommy mommy wet itself: {0:?}")]
    Nom(ErrorKind),
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

/// Isotope = "[" , Integer , Element , "]" ;
fn isotope_expr(i: &str) -> ParseResult<(MassNumber, &str)> {
    delimited(char('['), pair(count, element_symbol), char(']'))(i)
}

/// Element = uppercase , [ lowercase ] ;
fn element<'a>(db: &'a ChemicalDatabase) -> impl FnMut(&'a str) -> ParseResult<Element> {
    map_res(element_symbol, |symbol| Element::new(db, symbol))
}

/// Isotope = "[" , Integer , Element , "]" ;
fn isotope<'a>(db: &'a ChemicalDatabase) -> impl FnMut(&'a str) -> ParseResult<Element> {
    map_res(isotope_expr, |(mass_number, symbol)| {
        Element::new_isotope(db, symbol, mass_number)
    })
}

/// Particle = lowercase ;
fn particle<'a>(db: &'a ChemicalDatabase) -> impl FnMut(&'a str) -> ParseResult<Particle> {
    map_res(recognize(lowercase), |symbol| Particle::new(db, symbol))
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
    preceded(not(char('0')), u32)(i)
}

/// Particle Offset = Offset Kind , [ Integer ] ,
///   Particle ;
fn particle_offset<'a>(
    db: &'a ChemicalDatabase,
) -> impl FnMut(&'a str) -> ParseResult<(OffsetKind, Count, Particle)> {
    let optional_count = opt(count).map(|o| o.unwrap_or(1));
    tuple((offset_kind, optional_count, particle(db)))
}

/// Atomic Offset = ( Element | Isotope ) , [ Count ] ;
fn atomic_offset<'a>(
    db: &'a ChemicalDatabase,
) -> impl FnMut(&'a str) -> ParseResult<(Element, Count)> {
    let optional_count = opt(count).map(|o| o.unwrap_or(1));
    pair(alt((element(db), isotope(db))), optional_count)
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
    satisfy(|c| c.is_ascii_uppercase())(i)
        .map_err(|e| set_kind(e, CompositionErrorKind::ExpectedUppercase))
}

/// lowercase
///   = "a" | "b" | "c" | "d" | "e" | "f" | "g"
///   | "h" | "i" | "j" | "k" | "l" | "m" | "n"
///   | "o" | "p" | "q" | "r" | "s" | "t" | "u"
///   | "v" | "w" | "x" | "y" | "z"
///   ;
fn lowercase(i: &str) -> ParseResult<char> {
    satisfy(|c| c.is_ascii_lowercase())(i)
        .map_err(|e| set_kind(e, CompositionErrorKind::ExpectedLowercase))
}

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use nom::combinator::all_consuming;
    use once_cell::sync::Lazy;

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
        let mut chemical_composition = all_consuming(chemical_composition(&DB));
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
        eprintln!("What is up, my dude?\n\n\n");
        dbg!(chemical_composition("Au-2p"));
        panic!();
        // Multiple Chemical Compositions
        check_composition_snapshot!("[37Cl]5-2p+10", "+10");
        check_composition_snapshot!("[2H]2O+H2O", "+H2O");
        check_composition_snapshot!("NH2[100Tc]", "[100Tc]");
        check_composition_snapshot!("C11H12N2O2 H2O", " H2O");
    }
}
