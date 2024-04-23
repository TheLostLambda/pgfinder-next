use miette::Diagnostic;
use thiserror::Error;

use crate::{parsers::errors::CompositionError, FunctionalGroup};

pub type Result<T, E = Box<PolychemError>> = std::result::Result<T, E>;

// FIXME: Maybe there are too many layers of things being wrapped here!
// FIXME: Maybe just rename this to be `Error`?
// FIXME: Check all of the errors returned from public API are wrapped in this!
#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
pub enum PolychemError {
    #[error(transparent)]
    #[diagnostic(transparent)]
    Composition(#[from] CompositionError),

    #[error("the residue {0:?} could not be found in the supplied polymer database")]
    ResidueLookup(String),

    #[error("the modification {0:?} could not be found in the supplied polymer database")]
    ModificationLookup(String),

    #[error("the bond {0:?} could not be found in the supplied polymer database")]
    BondLookup(String),

    #[error("the functional group {0} could not be found on the residue {1} ({2})")]
    GroupLookup(String, String, String),
}

// FIXME: Move this to it's own errors.rs module? Bring the enum along too?
impl PolychemError {
    pub(crate) fn residue_lookup(abbr: &str) -> Self {
        Self::ResidueLookup(abbr.to_owned())
    }

    pub(crate) fn modification_lookup(abbr: &str) -> Self {
        Self::ModificationLookup(abbr.to_owned())
    }

    pub(crate) fn bond_lookup(abbr: &str) -> Self {
        Self::BondLookup(abbr.to_owned())
    }

    pub(crate) fn group_lookup(functional_group: FunctionalGroup, name: &str, abbr: &str) -> Self {
        Self::GroupLookup(
            functional_group.to_string(),
            name.to_owned(),
            abbr.to_owned(),
        )
    }
}
