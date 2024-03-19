use miette::Diagnostic;
use thiserror::Error;

use crate::MassNumber;

pub(crate) mod atomic_database;
pub mod chemical_composition;
pub mod chemical_composition_parser;
mod element;
mod particle;

// NOTE: Public so that other parsers using `chemical_composition` as a building block can inspect errors — I should
// consider finding a way to make this more private in the future, or at least mark it as #[non_exhaustive] so that it
// doesn't end up becoming SemVer nightmare...
#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
pub enum AtomicLookupError {
    #[diagnostic(help("double-check for typos, or add a new entry to the atomic database"))]
    #[error("the element {0:?} could not be found in the supplied atomic database")]
    Element(String),

    #[diagnostic(help("double-check for typos, or add a new entry to the atomic database"))]
    #[error(
        "the isotope \"{0}-{1}\" could not be found in the supplied atomic database, though the following {2} \
        isotopes were found: {3:?}"
    )]
    Isotope(String, MassNumber, String, Vec<MassNumber>),

    #[diagnostic(help("double-check for typos, or add a new entry to the atomic database"))]
    #[error("the particle {0:?} could not be found in the supplied atomic database")]
    Particle(String),

    #[diagnostic(help(
        "consider explicitly selecting the isotope to be used in mass calculations — e.g. [{}{}]", .2[0], .1
    ))]
    #[error("no natural abundance data could be found for {0} ({1}), though the following isotopes were found: {2:?}")]
    Abundance(String, String, Vec<MassNumber>),
}
