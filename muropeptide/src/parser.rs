// FIXME: Not 100% sold on have this as it's own module...

use miette::Diagnostic;
use nom::{error::ErrorKind, IResult};
use nom_miette::{FromExternalError, LabeledErrorKind, LabeledParseError};
use polychem::{polymerizer::Polymerizer, PolychemError};
use thiserror::Error;

use crate::Monomer;

// FIXME: Change to match the form of Chemical Composition
/// Monomer
///   = Glycan
///   | Peptide
///   | Glycan , "-" , Peptide
///   ;
fn monomer<'a, 's>(
    polymerizer: &mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<Monomer<'a>> {
    |_| todo!()
}

type ParseResult<'a, O> = IResult<&'a str, O, LabeledParseError<'a, MuropeptideErrorKind>>;

#[derive(Clone, Eq, PartialEq, Debug, Diagnostic, Error)]
pub enum MuropeptideErrorKind {
    #[diagnostic(help("double-check for typos, or add a new entry to the polymer database"))]
    #[error(transparent)]
    LookupError(PolychemError),

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

impl<'a> FromExternalError<'a, PolychemError> for MuropeptideErrorKind {
    const FATAL: bool = true;

    fn from_external_error(input: &'a str, e: PolychemError) -> LabeledParseError<'_, Self> {
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
