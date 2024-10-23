use miette::Diagnostic;
use nom::{error::ErrorKind, IResult};
use nom_miette::{FromExternalError, LabeledError, LabeledErrorKind, LabeledParseError};
use thiserror::Error;

use crate::{atoms::errors::AtomicLookupError, errors::PolychemError};

pub(crate) type CompositionError = LabeledError<PolychemErrorKind>;
pub type ParseResult<'a, O, K = PolychemErrorKind> = IResult<&'a str, O, LabeledParseError<'a, K>>;

pub trait UserErrorKind: From<PolychemErrorKind> + From<ErrorKind> {}
impl<T: From<PolychemErrorKind> + From<ErrorKind>> UserErrorKind for T {}

// NOTE: Public so that other parsers using `chemical_composition` as a building block can inspect errors — I should
// consider finding a way to make this more private in the future, or at least mark it as #[non_exhaustive] so that it
// doesn't end up becoming SemVer nightmare...
#[derive(Clone, Eq, PartialEq, Debug, Diagnostic, Error)]
pub enum PolychemErrorKind {
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

    // FIXME: This needs a label and testing!
    #[error("expected a count followed by 'x'")]
    ExpectedMultiplier,

    #[diagnostic(help(
        "a 0 value doesn't make sense here, if you've mistakenly included a leading zero, like \
        NH02, try just NH2 instead"
    ))]
    #[error("counts cannot start with 0")]
    ExpectedNoLeadingZero,

    #[error("expected an ASCII digit 1-9")]
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

    #[diagnostic(transparent)]
    #[error(transparent)]
    LookupError(Box<AtomicLookupError>),

    // FIXME: This needs a label and testing!
    #[diagnostic(transparent)]
    #[error(transparent)]
    PolychemError(#[from] Box<PolychemError>),

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

impl LabeledErrorKind for PolychemErrorKind {
    fn label(&self) -> Option<&'static str> {
        Some(match self {
            // NOTE: Stuck with this nested match until either `box_patterns` or `deref_patterns` are stabilized.
            // Keep an eye on:
            //   1) https://github.com/rust-lang/rust/issues/29641
            //   2) https://github.com/rust-lang/rust/issues/87121
            Self::LookupError(e) => match **e {
                AtomicLookupError::Element(..) => "element not found",
                AtomicLookupError::Isotope(..) => "isotope not found",
                AtomicLookupError::Particle(..) => "particle not found",
                AtomicLookupError::Abundance(..) => "no natural abundance",
            },
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

impl<'a> FromExternalError<'a, AtomicLookupError> for PolychemErrorKind {
    const FATAL: bool = true;

    fn from_external_error(input: &'a str, e: AtomicLookupError) -> LabeledParseError<'a, Self> {
        LabeledParseError::new(input, Self::LookupError(Box::new(e)))
    }
}

impl<'a> FromExternalError<'a, Box<PolychemError>> for PolychemErrorKind {
    const FATAL: bool = true;

    fn from_external_error(input: &'a str, e: Box<PolychemError>) -> LabeledParseError<'a, Self> {
        LabeledParseError::new(input, Self::PolychemError(e))
    }
}

impl From<ErrorKind> for PolychemErrorKind {
    fn from(value: ErrorKind) -> Self {
        match value {
            ErrorKind::Eof => Self::Incomplete,
            kind => Self::NomError(kind),
        }
    }
}
