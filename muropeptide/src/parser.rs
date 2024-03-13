use miette::Diagnostic;
use nom::{
    branch::alt,
    character::complete::{char, satisfy},
    combinator::recognize,
    error::ErrorKind,
    sequence::terminated,
    IResult,
};
use nom_miette::{
    expect, into, map_res, wrap_err, FromExternalError, LabeledErrorKind, LabeledParseError,
};
use polychem::{
    atoms::chemical_composition::{self, CompositionErrorKind},
    polymerizer::Polymerizer,
    AnyModification, ChemicalComposition, Count, Modification, NamedMod, OffsetKind, OffsetMod,
    PolychemError,
};
use thiserror::Error;

use crate::{AminoAcid, LateralChain, Monomer, Monosaccharide, UnbranchedAminoAcid};

// FIXME: Paste all of these EBNF comments into another file and make sure they are valid!
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
) -> impl FnMut(&'s str) -> ParseResult<Vec<Monosaccharide<'a>>> {
    |_| todo!()
}

/// Peptide = { Amino Acid }- ;
fn peptide<'a, 's>(
    polymerizer: &mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<Vec<AminoAcid<'a>>> {
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

/// Amino Acid = Unbranched Amino Acid , [ Lateral Chain ] ;
fn amino_acid<'a, 's>(
    polymerizer: &'a mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<AminoAcid<'a>> {
    |_| todo!()
}

// =

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

/// Unbranched Amino Acid = uppercase , [ Modifications ] ;
fn unbranched_amino_acid<'a, 's>(
    polymerizer: &'a mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<UnbranchedAminoAcid<'a>> {
    let parser = recognize(uppercase);
    map_res(parser, |abbr| polymerizer.residue(abbr))
}

/// Lateral Chain = "[" , [ "<" (* C-to-N *) | ">" (* N-to-C *) ] ,
///   { Unbranched Amino Acid }- , "]" ;
fn lateral_chain<'a, 's>(
    polymerizer: &'a mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<LateralChain> {
    |_| todo!()
}

// =

/// Predefined Modification = [ Multiplier ] , letter ,
///   { letter | digit | "_" } ;
fn predefined_modification<'a, 's>(
    polymerizer: &'a mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<Modification<NamedMod<'a, 'a>>> {
    |_| todo!()
}

/// Chemical Offset = Offset Kind , [ Multiplier ] ,
///   Chemical Composition ;
fn chemical_offset<'a, 's>(
    polymerizer: &'a mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<Modification<OffsetMod<'a>>> {
    |_| todo!()
}

// =

/// Multiplier = Count , "x" ;
pub fn multiplier(i: &str) -> ParseResult<Count> {
    let parser = terminated(count, char('x'));
    wrap_err(parser, MuropeptideErrorKind::ExpectedMultiplier)(i)
}

// Adapted parsers =

fn chemical_composition<'a, 's>(
    polymerizer: &'a mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<ChemicalComposition<'a>> {
    into(chemical_composition::chemical_composition(
        polymerizer.atomic_db(),
    ))
}

macro_rules! wrap_composition_parsers {
    ($($f:ident -> $t:ty),+ $(,)?) => {
        $(
            fn $f(i: &str) -> ParseResult<$t> {
                into(chemical_composition::$f)(i)
            }
        )+
    };
}

wrap_composition_parsers!(
    count -> Count,
    offset_kind -> OffsetKind,
    uppercase -> char,
    lowercase -> char,
);

type ParseResult<'a, O> = IResult<&'a str, O, LabeledParseError<'a, MuropeptideErrorKind>>;

#[derive(Clone, Eq, PartialEq, Debug, Diagnostic, Error)]
pub enum MuropeptideErrorKind {
    #[error("expected a count followed by 'x'")]
    ExpectedMultiplier,

    #[diagnostic(transparent)]
    #[error(transparent)]
    PolychemError(Box<PolychemError>),

    #[diagnostic(transparent)]
    #[error(transparent)]
    CompositionError(#[from] CompositionErrorKind),

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
            // FIXME: Need to add branches for passing labels through the transparent errors?
            Self::Incomplete => "input was valid up until this point",
            Self::NomError(_) => "the region that triggered this bug!",
            _ => return None,
        })
    }
}

impl<'a> FromExternalError<'a, Box<PolychemError>> for MuropeptideErrorKind {
    const FATAL: bool = true;

    fn from_external_error(input: &'a str, e: Box<PolychemError>) -> LabeledParseError<'_, Self> {
        LabeledParseError::new(input, Self::PolychemError(e))
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
    fn test_multiplier() {
        // Valid Multipliers
        assert_eq!(multiplier("1x"), Ok(("", 1)));
        assert_eq!(multiplier("10x"), Ok(("", 10)));
        assert_eq!(multiplier("422x"), Ok(("", 422)));
        assert_eq!(multiplier("9999x"), Ok(("", 9999)));
        // Invalid Multipliers
        assert!(multiplier("1").is_err());
        assert!(multiplier("10").is_err());
        assert!(multiplier("422").is_err());
        assert!(multiplier("9999").is_err());
        assert!(multiplier("0").is_err());
        assert!(multiplier("01").is_err());
        assert!(multiplier("00145").is_err());
        assert!(multiplier("H").is_err());
        assert!(multiplier("p").is_err());
        assert!(multiplier("+H").is_err());
        assert!(multiplier("[H]").is_err());
        // Multiple Multipliers
        assert_eq!(multiplier("1xOH"), Ok(("OH", 1)));
        assert_eq!(multiplier("42xHeH"), Ok(("HeH", 42)));
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

    // FIXME: Add a test that checks all of the errors using `assert_miette_snapshot`! Maybe make that a crate?
}
