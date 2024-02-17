// External Crate Imports
use miette::Diagnostic;
use nom::{
    branch::alt,
    character::complete::{char, one_of, satisfy, u32},
    combinator::{cut, map, not, opt, recognize},
    error::ErrorKind,
    multi::many1,
    sequence::{delimited, pair, preceded},
    IResult, Parser,
};
use nom_miette::{
    expect, map_res, wrap_err, FromExternalError, LabeledError, LabeledErrorKind, LabeledParseError,
};
use thiserror::Error;

// Local Module Imports
use super::{
    atomic_database::AtomicDatabase, AtomicLookupError, ChemicalComposition, Count, Element,
    MassNumber, OffsetKind, Particle,
};

// Public API ==========================================================================================================

pub type CompositionError = LabeledError<CompositionErrorKind>;

/// Chemical Composition
///   = { Atomic Offset }- , [ Offset Kind , Particle Offset ]
///   | Particle Offset
///   ;
pub fn chemical_composition<'a>(
    db: &'a AtomicDatabase,
) -> impl FnMut(&'a str) -> ParseResult<ChemicalComposition> {
    let chemical_formula = many1(atomic_offset(db));
    let optional_particle_offset = opt(pair(offset_kind, cut(particle_offset(db))));
    let atoms_and_particles = map(
        pair(chemical_formula, optional_particle_offset),
        |(chemical_formula, particle_offset)| {
            let particle_offset = particle_offset.map(|(k, (c, p))| (k, c, p));
            ChemicalComposition {
                chemical_formula,
                particle_offset,
            }
        },
    );

    let just_particles = map(particle_offset(db), |(count, particle)| {
        let particle_offset = Some((OffsetKind::Add, count, particle));
        ChemicalComposition {
            chemical_formula: Vec::new(),
            particle_offset,
        }
    });

    let parser = alt((atoms_and_particles, just_particles));
    wrap_err(parser, CompositionErrorKind::ExpectedChemicalComposition)
}

// Private Sub-Parsers =================================================================================================

/// Atomic Offset = ( Element | Isotope ) , [ Count ] ;
fn atomic_offset<'a>(
    db: &'a AtomicDatabase,
) -> impl FnMut(&'a str) -> ParseResult<(Element, Count)> {
    let element_or_isotope = alt((element(db), isotope(db)));
    let optional_count = opt(count).map(|o| o.unwrap_or(1));
    let parser = pair(element_or_isotope, optional_count);
    wrap_err(parser, CompositionErrorKind::ExpectedAtomicOffset)
}

/// Offset Kind = "+" | "-" ;
fn offset_kind(i: &str) -> ParseResult<OffsetKind> {
    map(one_of("+-"), |c| match c {
        '+' => OffsetKind::Add,
        '-' => OffsetKind::Remove,
        _ => unreachable!(),
    })(i)
}

/// Particle Offset = [ Count ] , Particle ;
fn particle_offset<'a>(
    db: &'a AtomicDatabase,
) -> impl FnMut(&'a str) -> ParseResult<(Count, Particle)> {
    let optional_count = opt(count).map(|o| o.unwrap_or(1));
    let parser = pair(optional_count, particle(db));
    wrap_err(parser, CompositionErrorKind::ExpectedParticleOffset)
}

// ---------------------------------------------------------------------------------------------------------------------

/// Element = uppercase , [ lowercase ] ;
fn element<'a>(db: &'a AtomicDatabase) -> impl FnMut(&'a str) -> ParseResult<Element> {
    map_res(element_symbol, |symbol| Element::new(db, symbol))
}

/// Isotope = "[" , Count , Element , "]" ;
fn isotope<'a>(db: &'a AtomicDatabase) -> impl FnMut(&'a str) -> ParseResult<Element> {
    map_res(isotope_expr, |(mass_number, symbol)| {
        Element::new_isotope(db, symbol, mass_number)
    })
}

/// Count = digit - "0" , { digit } ;
fn count(i: &str) -> ParseResult<Count> {
    let not_zero = expect(
        cut(not(char('0'))),
        CompositionErrorKind::ExpectedNoLeadingZero,
    );
    let digits = expect(u32, CompositionErrorKind::ExpectedDigit);
    preceded(not_zero, digits)(i)
}

/// Particle = lowercase ;
fn particle<'a>(db: &'a AtomicDatabase) -> impl FnMut(&'a str) -> ParseResult<Particle> {
    map_res(particle_symbol, |symbol| Particle::new(db, symbol))
}

// ---------------------------------------------------------------------------------------------------------------------

/// Element = uppercase , [ lowercase ] ;
fn element_symbol(i: &str) -> ParseResult<&str> {
    let parser = recognize(pair(uppercase, opt(lowercase)));
    wrap_err(parser, CompositionErrorKind::ExpectedElementSymbol)(i)
}

/// Isotope = "[" , Count , Element , "]" ;
fn isotope_expr(i: &str) -> ParseResult<(MassNumber, &str)> {
    let opening_bracket = expect(char('['), CompositionErrorKind::ExpectedIsotopeStart);
    let mass_number = wrap_err(count, CompositionErrorKind::ExpectedMassNumber);
    let closing_bracket = expect(cut(char(']')), CompositionErrorKind::ExpectedIsotopeEnd);
    delimited(
        opening_bracket,
        cut(pair(mass_number, element_symbol)),
        closing_bracket,
    )(i)
}

/// Particle = lowercase ;
fn particle_symbol(i: &str) -> ParseResult<&str> {
    let parser = recognize(lowercase);
    wrap_err(parser, CompositionErrorKind::ExpectedParticleSymbol)(i)
}

// ---------------------------------------------------------------------------------------------------------------------

/// uppercase
///   = "A" | "B" | "C" | "D" | "E" | "F" | "G"
///   | "H" | "I" | "J" | "K" | "L" | "M" | "N"
///   | "O" | "P" | "Q" | "R" | "S" | "T" | "U"
///   | "V" | "W" | "X" | "Y" | "Z"
///   ;
fn uppercase(i: &str) -> ParseResult<char> {
    let parser = satisfy(|c| c.is_ascii_uppercase());
    expect(parser, CompositionErrorKind::ExpectedUppercase)(i)
}

/// lowercase
///   = "a" | "b" | "c" | "d" | "e" | "f" | "g"
///   | "h" | "i" | "j" | "k" | "l" | "m" | "n"
///   | "o" | "p" | "q" | "r" | "s" | "t" | "u"
///   | "v" | "w" | "x" | "y" | "z"
///   ;
fn lowercase(i: &str) -> ParseResult<char> {
    let parser = satisfy(|c| c.is_ascii_lowercase());
    expect(parser, CompositionErrorKind::ExpectedLowercase)(i)
}

// Parse Error Types and Trait Implementations =========================================================================

type ParseResult<'a, O> = IResult<&'a str, O, LabeledParseError<'a, CompositionErrorKind>>;

#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
pub(super) enum CompositionErrorKind {
    #[error(
        "expected a chemical formula (optionally followed by a '+' or '-' and a particle offset), \
        or a standalone particle offset"
    )]
    ExpectedChemicalComposition,

    #[error(
        "expected an element (like Au) or an isotope (like [15N]) optionally followed by a number"
    )]
    ExpectedAtomicOffset,

    #[error("expected a particle (like p or e), optionally preceded by a number")]
    ExpectedParticleOffset,

    #[diagnostic(help(
        "a 0 value doesn't make sense here, if you've mistakenly included a leading zero, like \
        NH02, try just NH2 instead"
    ))]
    #[error("counts cannot start with 0")]
    ExpectedNoLeadingZero,

    #[error("expected a digit 1-9")]
    ExpectedDigit,

    #[error("expected an element symbol")]
    ExpectedElementSymbol,

    #[error("expected '[' to open isotope brackets")]
    ExpectedIsotopeStart,

    #[error("expected an isotopic mass number")]
    ExpectedMassNumber,

    #[diagnostic(help("you've probably forgotten to close an earlier '[' bracket"))]
    #[error("expected ']' to close isotope brackets")]
    ExpectedIsotopeEnd,

    #[error("expected a particle symbol")]
    ExpectedParticleSymbol,

    #[error("expected an uppercase ASCII letter")]
    ExpectedUppercase,

    #[error("expected a lowercase ASCII letter")]
    ExpectedLowercase,

    #[diagnostic(help("double-check for typos, or add a new entry to the atomic database"))]
    #[error(transparent)]
    LookupError(AtomicLookupError),

    #[diagnostic(help(
        "this is an internal error that you shouldn't ever see! If you have gotten this error, \
        then please report it as a bug!"
    ))]
    #[error("internal `nom` error: {0:?}")]
    NomError(ErrorKind),

    #[diagnostic(help(
        "check the unparsed region for errors, or remove it from the rest of the composition"
    ))]
    #[error("could not interpret the full input as a valid chemical composition")]
    Incomplete,
}

impl LabeledErrorKind for CompositionErrorKind {
    fn label(&self) -> Option<&'static str> {
        Some(match self {
            Self::LookupError(AtomicLookupError::Element(_)) => "element not found",
            Self::LookupError(AtomicLookupError::Isotope(_, _, _, _)) => "isotope not found",
            Self::LookupError(AtomicLookupError::Particle(_)) => "particle not found",
            Self::ExpectedUppercase => "expected uppercase",
            Self::ExpectedLowercase => "expected lowercase",
            Self::ExpectedDigit => "expected digit",
            Self::ExpectedIsotopeStart => "'['",
            Self::ExpectedIsotopeEnd => "expected ']'",
            Self::ExpectedMassNumber => "expected a mass number",
            Self::ExpectedNoLeadingZero => "expected non-zero",
            Self::Incomplete => "input was valid up until this point",
            Self::NomError(_) => "the region that triggered this bug!",
            _ => return None,
        })
    }
}

impl<'a> FromExternalError<'a, AtomicLookupError> for LabeledParseError<'a, CompositionErrorKind> {
    const FATAL: bool = true;

    fn from_external_error(input: &'a str, e: AtomicLookupError) -> Self {
        Self::new(input, CompositionErrorKind::LookupError(e))
    }
}

impl From<ErrorKind> for CompositionErrorKind {
    fn from(value: ErrorKind) -> Self {
        match value {
            ErrorKind::Eof => Self::Incomplete,
            kind => Self::NomError(kind),
        }
    }
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use nom_miette::final_parser;
    use once_cell::sync::Lazy;

    use crate::testing_tools::assert_miette_snapshot;

    use super::*;

    static DB: Lazy<AtomicDatabase> = Lazy::new(|| {
        AtomicDatabase::from_kdl(
            "atomic_database.kdl",
            include_str!("../atomic_database.kdl"),
        )
        .unwrap()
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
        // Check that multiple labels are reported for errors with different spans
        assert_miette_snapshot!(chemical_composition("2H2"));
    }
}
