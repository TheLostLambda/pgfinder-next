// FIXME: Not 100% sold on have this as it's own module...

use miette::Diagnostic;
use nom::{character::complete::satisfy, combinator::recognize, error::ErrorKind, IResult};
use nom_miette::{
    expect, map_res, wrap_err, FromExternalError, LabeledErrorKind, LabeledParseError,
};
use polychem::{
    polymerizer::{self, Polymerizer},
    AnyModification, PolychemError, Residue,
};
use thiserror::Error;

use crate::{AminoAcid, Monomer, Monosaccharide};

/// Monomer = Glycan , [ "-" , Peptide ] | Peptide ;
fn monomer<'a, 's>(
    polymerizer: &mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<Monomer<'a>> {
    |_| todo!()
}

// =

/// Glycan = { Monosaccharide }- ;
fn glycan<'a, 's>(
    polymerizer: &mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<Vec<Residue<'a, 'a>>> {
    |_| todo!()
}

/// Peptide = { Amino Acid }- ;
fn peptide<'a, 's>(
    polymerizer: &mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<AminoAcid<'a>> {
    |_| todo!()
}

// =

// FIXME: Don't know if it's a great idea to tie together the lifetimes of the chemical databases and the polymerizer
// instance here? Everything is using 'a...
// FIXME: Add modifications
/// Monosaccharide = lowercase , [ Modifications ] ;
fn monosaccharide<'a, 's>(
    polymerizer: &'a mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<Monosaccharide<'a>> {
    let parser = recognize(lowercase);
    map_res(parser, |abbr| polymerizer.residue(abbr))
}

// FIXME: Add modifications and lateral chains
/// Amino Acid = uppercase , [ Modifications ] , [ Lateral Chain ]
fn amino_acid<'a, 's>(
    polymerizer: &'a mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<Residue<'a, 'a>> {
    let parser = recognize(uppercase);
    map_res(parser, |abbr| polymerizer.residue(abbr))
}

// =

// FIXME: Can this EBNF be made simpler? Or not with comma separators...
/// Modifications = "(" ,
///   ( Predefined Modification
///   | Chemical Offset
///   ) , { { " " } , "," , { " " } ,
///   ( Predefined Modification
///   | Chemical Offset
///   ) } , ")" ;
fn modifications<'a, 's>(
    polymerizer: &'a mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<Vec<AnyModification<'a, 'a>>> {
    |_| todo!()
}

/// uppercase
///   = "A" | "B" | "C" | "D" | "E" | "F" | "G"
///   | "H" | "I" | "J" | "K" | "L" | "M" | "N"
///   | "O" | "P" | "Q" | "R" | "S" | "T" | "U"
///   | "V" | "W" | "X" | "Y" | "Z"
///   ;
fn uppercase(i: &str) -> ParseResult<char> {
    let parser = satisfy(|c| c.is_ascii_uppercase());
    expect(parser, MuropeptideErrorKind::ExpectedUppercase)(i)
}

/// lowercase
///   = "a" | "b" | "c" | "d" | "e" | "f" | "g"
///   | "h" | "i" | "j" | "k" | "l" | "m" | "n"
///   | "o" | "p" | "q" | "r" | "s" | "t" | "u"
///   | "v" | "w" | "x" | "y" | "z"
///   ;
fn lowercase(i: &str) -> ParseResult<char> {
    let parser = satisfy(|c| c.is_ascii_lowercase());
    expect(parser, MuropeptideErrorKind::ExpectedLowercase)(i)
}

type ParseResult<'a, O> = IResult<&'a str, O, LabeledParseError<'a, MuropeptideErrorKind>>;

#[derive(Clone, Eq, PartialEq, Debug, Diagnostic, Error)]
pub enum MuropeptideErrorKind {
    #[error("expected an uppercase ASCII letter")]
    ExpectedUppercase,

    #[error("expected a lowercase ASCII letter")]
    ExpectedLowercase,

    #[diagnostic(help("double-check for typos, or add a new entry to the polymer database"))]
    #[error(transparent)]
    LookupError(Box<PolychemError>),

    #[diagnostic(help(
        "this is an internal error that you shouldn't ever see! If you have gotten this error, \
        then please report it as a bug!"
    ))]
    #[error("internal `nom` error: {0:?}")]
    NomError(ErrorKind),

    #[diagnostic(help(
        "check the unparsed region for errors, or remove it from the rest of the muropeptide"
    ))]
    #[error("could not interpret the full input as a valid muropeptide structure")]
    Incomplete,
}

impl LabeledErrorKind for MuropeptideErrorKind {
    fn label(&self) -> Option<&'static str> {
        Some(match self {
            // FIXME: Replace with `PolychemError`s
            // Self::LookupError(AtomicLookupError::Element(_)) => "element not found",
            // Self::LookupError(AtomicLookupError::Isotope(_, _, _, _)) => "isotope not found",
            // Self::LookupError(AtomicLookupError::Particle(_)) => "particle not found",
            Self::Incomplete => "input was valid up until this point",
            Self::NomError(_) => "the region that triggered this bug!",
            _ => return None,
        })
    }
}

impl<'a> FromExternalError<'a, Box<PolychemError>> for MuropeptideErrorKind {
    const FATAL: bool = true;

    fn from_external_error(input: &'a str, e: Box<PolychemError>) -> LabeledParseError<'_, Self> {
        LabeledParseError::new(input, Self::LookupError(e))
    }
}

impl From<ErrorKind> for MuropeptideErrorKind {
    fn from(value: ErrorKind) -> Self {
        match value {
            ErrorKind::Eof => Self::Incomplete,
            kind => Self::NomError(kind),
        }
    }
}

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;
    use polychem::{AtomicDatabase, PolymerDatabase};

    use super::*;

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(|| {
        AtomicDatabase::from_kdl(
            "atomic_database.kdl",
            include_str!("../../polychem/atomic_database.kdl"),
        )
        .unwrap()
    });

    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::from_kdl(
            &ATOMIC_DB,
            "muropeptide_chemistry.kdl",
            include_str!("../../polychem/muropeptide_chemistry.kdl"),
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
    fn test_monosaccharide() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut monosaccharide = monosaccharide(&mut polymerizer);
        macro_rules! assert_monosaccharide_name {
            ($input:literal, $output:literal, $name:literal) => {
                assert_eq!(
                    monosaccharide($input).map(|(r, e)| (r, e.name())),
                    Ok(($output, $name))
                );
            };
        }
        // Valid Monosaccharides
        assert_monosaccharide_name!("g", "", "N-Acetylglucosamine");
        assert_monosaccharide_name!("m", "", "N-Acetylmuramic Acid");
        // Invalid Monosaccharides
        assert!(monosaccharide("P").is_err());
        assert!(monosaccharide("EP").is_err());
        assert!(monosaccharide("1h").is_err());
        assert!(monosaccharide("+m").is_err());
        assert!(monosaccharide("-g").is_err());
        assert!(monosaccharide("[h]").is_err());
        // Non-Existent Monosaccharides
        assert!(monosaccharide("s").is_err());
        assert!(monosaccharide("f").is_err());
        // Multiple Monosaccharides
        assert_monosaccharide_name!("gm", "m", "N-Acetylglucosamine");
        assert_monosaccharide_name!("m-A", "-A", "N-Acetylmuramic Acid");
    }
}
