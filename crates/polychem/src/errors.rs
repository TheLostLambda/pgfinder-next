use miette::Diagnostic;
use thiserror::Error;

use crate::{
    parsers::errors::CompositionError, polymers::errors::FindFreeGroupsError, FunctionalGroup,
    ResidueId,
};

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

    #[error("residue {id} could not be found in the current polymer")]
    #[diagnostic(help("this residue may belong to another polymer, or may have been previously deleted from this one"))]
    ResidueNotInPolymer { id: ResidueId },

    #[error(
        "failed to form {name} bond between residue {donor_id} ({donor_name}) and residue {acceptor_id} \
        ({acceptor_name}) due to an issue with the {donor_or_acceptor}"
    )]
    Bond {
        name: String,
        donor_id: ResidueId,
        donor_name: String,
        acceptor_id: ResidueId,
        acceptor_name: String,
        donor_or_acceptor: String,
        #[source]
        #[diagnostic_source]
        source: FindFreeGroupsError,
    },

    #[error("failed to apply {name} modification to residue {residue_id} ({residue_name})")]
    NamedModification {
        name: String,
        residue_id: ResidueId,
        residue_name: String,
        #[source]
        #[diagnostic_source]
        source: FindFreeGroupsError,
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

    // FIXME: Perhaps all of my error constructors should be marked with #[must_use]... Or the type should be? Clippy
    // has only caught this one, but I need to make this more consistent!
    #[must_use]
    pub(crate) const fn residue_not_in_polymer(id: ResidueId) -> Self {
        Self::ResidueNotInPolymer { id }
    }

    pub(crate) fn bond(
        name: &str,
        donor_id: ResidueId,
        donor_name: &str,
        acceptor_id: ResidueId,
        acceptor_name: &str,
        donor_or_acceptor: &str,
        source: FindFreeGroupsError,
    ) -> Self {
        let name = name.to_owned();
        let donor_name = donor_name.to_owned();
        let acceptor_name = acceptor_name.to_owned();
        let donor_or_acceptor = donor_or_acceptor.to_owned();

        Self::Bond {
            name,
            donor_id,
            donor_name,
            acceptor_id,
            acceptor_name,
            donor_or_acceptor,
            source,
        }
    }

    pub(crate) fn named_modification(
        name: &str,
        residue_id: ResidueId,
        residue_name: &str,
        source: FindFreeGroupsError,
    ) -> Self {
        let name = name.to_owned();
        let residue_name = residue_name.to_owned();

        Self::NamedModification {
            name,
            residue_id,
            residue_name,
            source,
        }
    }
}
