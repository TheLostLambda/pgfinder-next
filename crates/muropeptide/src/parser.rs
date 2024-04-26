use miette::Diagnostic;
use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::{alpha1, alphanumeric1, char, space0},
    combinator::{opt, recognize},
    error::ErrorKind,
    multi::{many0, separated_list1},
    sequence::{delimited, pair},
    IResult,
};
use nom_miette::{wrap_err, FromExternalError, LabeledErrorKind, LabeledParseError};
use polychem::{
    errors::PolychemError,
    parsers::{errors::PolychemErrorKind, modifications, primitives::uppercase},
    ModificationId, Polymer,
};
use thiserror::Error;

use crate::{AminoAcid, LateralChain, Monomer, Monosaccharide, UnbranchedAminoAcid};

// FIXME: Paste all of these EBNF comments into another file and make sure they are valid!
/// Monomer = Glycan , [ "-" , Peptide ] | Peptide ;
fn monomer<'a, 'p, 's>(
    _polymer: &mut Polymer<'a, 'p>,
) -> impl FnMut(&'s str) -> ParseResult<Monomer> {
    |_| todo!()
}

// =

/// Glycan = { Monosaccharide }- ;
fn glycan<'a, 'p, 's>(
    _polymer: &mut Polymer<'a, 'p>,
) -> impl FnMut(&'s str) -> ParseResult<Vec<Monosaccharide>> {
    |_| todo!()
}

/// Peptide = { Amino Acid }- ;
fn peptide<'a, 'p, 's>(
    _polymer: &mut Polymer<'a, 'p>,
) -> impl FnMut(&'s str) -> ParseResult<Vec<AminoAcid>> {
    |_| todo!()
}

// =

// FIXME: Don't know if it's a great idea to tie together the lifetimes of the chemical databases and the polymerizer
// instance here? Everything is using 'a...
// FIXME: Add modifications
/// Monosaccharide = lowercase , [ Modifications ] ;
fn monosaccharide<'a, 'p, 's>(
    _polymer: &'a mut Polymer<'a, 'p>,
) -> impl FnMut(&'s str) -> ParseResult<Monosaccharide> {
    |_| todo!()
}

/// Amino Acid = Unbranched Amino Acid , [ Lateral Chain ] ;
fn amino_acid<'a, 'p, 's>(
    _polymer: &mut Polymer<'a, 'p>,
) -> impl FnMut(&'s str) -> ParseResult<AminoAcid> {
    |_| todo!()
}

// =

/// Modifications = "(" , Any Modification ,
///   { { " " } , "," , { " " } , Any Modification } , ")" ;
fn modifications<'a, 'p, 's>(
    polymer: &mut Polymer<'a, 'p>,
) -> impl FnMut(&'s str) -> ParseResult<Vec<ModificationId>> {
    let separator = delimited(space0, char(','), space0);
    delimited(
        char('('),
        separated_list1(separator, modifications::any(polymer, identifier)),
        char(')'),
    )
}

// FIXME: Make private again!
// FIXME: Switch to a more efficient modification application API
/// Unbranched Amino Acid = uppercase , [ Modifications ] ;
pub fn unbranched_amino_acid<'a, 'p, 's>(
    polymer: &'a mut Polymer<'a, 'p>,
) -> impl FnMut(&'s str) -> ParseResult<UnbranchedAminoAcid> {
    let _parser = pair(recognize(uppercase), opt(modifications(polymer)));
    // TODO: Apply modifications to the residues they follow here!
    |_| todo!()
}

// NOTE: These are not meant to be links, it's just EBNF
#[allow(clippy::doc_link_with_quotes)]
/// Lateral Chain = "[" , [ "<" (* C-to-N *) | ">" (* N-to-C *) ] ,
///   { Unbranched Amino Acid }- , "]" ;
fn lateral_chain<'a, 'p, 's>(
    _polymer: &mut Polymer<'a, 'p>,
) -> impl FnMut(&'s str) -> ParseResult<LateralChain> {
    |_| todo!()
}

// =

/// Identifier = letter , { letter | digit | "_" } ;
fn identifier(i: &str) -> ParseResult<&str> {
    // PERF: Could maybe avoid allocations by using `many0_count` instead, but needs benchmarking
    let parser = recognize(pair(alpha1, many0(alt((alphanumeric1, tag("_"))))));
    wrap_err(parser, MuropeptideErrorKind::ExpectedIdentifier)(i)
}

type ParseResult<'a, O> = IResult<&'a str, O, LabeledParseError<'a, MuropeptideErrorKind>>;

#[derive(Clone, Eq, PartialEq, Debug, Diagnostic, Error)]
pub enum MuropeptideErrorKind {
    #[error("expected an ASCII letter, optionally followed by any number of ASCII letters, digits, and underscores")]
    ExpectedIdentifier,

    // FIXME: Kill this and merge into the error below!
    #[diagnostic(transparent)]
    #[error(transparent)]
    PolychemError(Box<PolychemError>),

    #[diagnostic(transparent)]
    #[error(transparent)]
    CompositionError(#[from] PolychemErrorKind),

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

// FIXME: Can I get rid of this?
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

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);

    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "polymer_database.kdl",
            include_str!("../data/polymer_database.kdl"),
        )
        .unwrap()
    });

    #[ignore]
    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_modifications() {
        // TODO: Restore from git!
        todo!();
    }

    #[test]
    fn test_identifier() {
        // Valid Identifiers
        assert_eq!(identifier("Ac"), Ok(("", "Ac")));
        assert_eq!(identifier("H2O"), Ok(("", "H2O")));
        assert_eq!(identifier("Anh"), Ok(("", "Anh")));
        assert_eq!(identifier("E2E"), Ok(("", "E2E")));
        assert_eq!(identifier("no_way"), Ok(("", "no_way")));
        assert_eq!(identifier("H"), Ok(("", "H")));
        assert_eq!(identifier("p"), Ok(("", "p")));
        // Invalid Identifiers
        assert!(identifier(" H2O").is_err());
        assert!(identifier("1").is_err());
        assert!(identifier("9999").is_err());
        assert!(identifier("0").is_err());
        assert!(identifier("00145").is_err());
        assert!(identifier("+H").is_err());
        assert!(identifier("[H]").is_err());
        assert!(identifier("Øof").is_err());
        assert!(identifier("2xAc").is_err());
        assert!(identifier("-Ac").is_err());
        assert!(identifier("_Ac").is_err());
        // Multiple Identifiers
        assert_eq!(identifier("OH-p"), Ok(("-p", "OH")));
        assert_eq!(identifier("HeH 2slow"), Ok((" 2slow", "HeH")));
        assert_eq!(identifier("Gefählt"), Ok(("ählt", "Gef")));
        // This is a weird unicode 6
        assert!('𝟨'.is_numeric());
        assert!(!'𝟨'.is_ascii_digit());
        assert_eq!(identifier("C2H𝟨O"), Ok(("𝟨O", "C2H")));
    }

    #[ignore]
    #[test]
    fn test_monosaccharide() {
        // TODO: Restore from git!
        todo!();
    }

    // FIXME: Unfininshed! Needs modification support — same with monosaccharide!
    #[ignore]
    #[test]
    fn test_unbranched_amino_acid() {
        // TODO: Restore from git!
        todo!();
    }

    // FIXME: Add a test that checks all of the errors using `assert_miette_snapshot`! Maybe make that a crate?
}
