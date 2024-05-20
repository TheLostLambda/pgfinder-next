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
    Composition {
        #[from]
        error: CompositionError,
    },

    #[error("the residue {abbr:?} could not be found in the supplied polymer database")]
    ResidueLookup { abbr: String },

    #[error("the modification {abbr:?} could not be found in the supplied polymer database")]
    ModificationLookup { abbr: String },

    #[error("the bond {abbr:?} could not be found in the supplied polymer database")]
    BondLookup { abbr: String },

    #[error("the functional group {group_name} could not be found on the residue {name} ({abbr})")]
    GroupLookup {
        group_name: String,
        name: String,
        abbr: String,
    },
}

impl PolychemError {
    pub(crate) fn residue_lookup(abbr: &str) -> Self {
        let abbr = abbr.to_owned();

        Self::ResidueLookup { abbr }
    }

    pub(crate) fn modification_lookup(abbr: &str) -> Self {
        let abbr = abbr.to_owned();

        Self::ModificationLookup { abbr }
    }

    pub(crate) fn bond_lookup(abbr: &str) -> Self {
        let abbr = abbr.to_owned();

        Self::BondLookup { abbr }
    }

    pub(crate) fn group_lookup(group: FunctionalGroup, name: &str, abbr: &str) -> Self {
        let group_name = group.to_string();
        let name = name.to_owned();
        let abbr = abbr.to_owned();

        Self::GroupLookup {
            group_name,
            name,
            abbr,
        }
    }
}
