use ahash::HashMap;
use miette::Diagnostic;
use thiserror::Error;

use crate::{Isotope, MassNumber};

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
        isotopes were found: {3}"
    )]
    Isotope(String, MassNumber, String, String),

    #[diagnostic(help("double-check for typos, or add a new entry to the atomic database"))]
    #[error("the particle {0:?} could not be found in the supplied atomic database")]
    Particle(String),

    #[diagnostic(help(
        "consider explicitly selecting the isotope to be used in mass calculations — e.g. [{2}{1}]"
    ))]
    #[error("no natural abundance data could be found for {0} ({1}), though the following isotopes were found: {3}")]
    Abundance(String, String, MassNumber, String),
}

impl AtomicLookupError {
    pub(crate) fn element(symbol: &str) -> Self {
        Self::Element(symbol.to_owned())
    }

    pub(crate) fn isotope(
        symbol: &str,
        mass_number: MassNumber,
        name: &str,
        isotopes: &HashMap<MassNumber, Isotope>,
    ) -> Self {
        Self::Isotope(
            symbol.to_owned(),
            mass_number,
            name.to_owned(),
            Self::display_vec(isotopes.keys()),
        )
    }

    pub(crate) fn particle(symbol: &str) -> Self {
        Self::Particle(symbol.to_owned())
    }

    pub(crate) fn abundance(
        symbol: &str,
        name: &str,
        isotopes: &HashMap<MassNumber, Isotope>,
    ) -> Self {
        Self::Abundance(
            symbol.to_owned(),
            name.to_owned(),
            // SAFETY: Validation of the `AtomicDatabase` ensures there is always at least one isotope per element
            *isotopes.keys().min().unwrap(),
            Self::display_vec(isotopes.keys()),
        )
    }

    // FIXME: Where does this belong?
    pub(crate) fn display_vec<I: ToString>(items: impl IntoIterator<Item = I>) -> String {
        let mut items: Vec<_> = items.into_iter().map(|i| i.to_string()).collect();
        items.sort_unstable();
        format!("[{}]", items.join(", "))
    }
}
