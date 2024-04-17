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

    #[error("the bond kind {0:?} could not be found in the supplied polymer database")]
    BondLookup(String),

    #[error("the functional group {0} could not be found on the residue {1} ({2})")]
    GroupLookup(String, String, String),

    #[error("failed to apply the offset modification {0} to residue {1} ({2})")]
    OffsetModification(
        String,
        Id,
        String,
        #[source]
        #[diagnostic_source]
        // FIXME: This should be hidden behind a private error struct, like `polymerizer::Error`
        OffsetMultiplierError,
    ),

    #[error("failed to apply the named modification {0} ({1}) to residue {2} ({3})")]
    NamedModification(
        String,
        String,
        Id,
        String,
        #[source]
        #[diagnostic_source]
        polymerizer::Error,
    ),

    #[error("failed to form {0:?} bond between residue {1} ({2}) and residue {3} ({4}) due to an issue with the {5}")]
    Bond(
        String,
        Id,
        String,
        Id,
        String,
        String,
        #[source]
        #[diagnostic_source]
        polymerizer::Error,
    ),
}

// FIXME: Move this to it's own errors.rs module? Bring the enum along too?
impl PolychemError {
    fn residue_lookup(abbr: &str) -> Self {
        Self::ResidueLookup(abbr.to_owned())
    }

    fn modification_lookup(abbr: &str) -> Self {
        Self::ModificationLookup(abbr.to_owned())
    }

    fn bond_lookup(kind: &str) -> Self {
        Self::BondLookup(kind.to_owned())
    }

    fn group_lookup(functional_group: FunctionalGroup, name: &str, abbr: &str) -> Self {
        Self::GroupLookup(
            functional_group.to_string(),
            name.to_owned(),
            abbr.to_owned(),
        )
    }

    fn offset_modification(
        count: Count,
        kind: OffsetKind,
        composition: ChemicalComposition,
        residue: &Residue,
        source: OffsetMultiplierError,
    ) -> Self {
        let modification =
            Modification::new(count, Offset::new_with_composition(kind, composition));
        Self::OffsetModification(
            modification.to_string(),
            residue.id(),
            residue.name().to_owned(),
            source,
        )
    }

    fn named_modification(
        name: &str,
        abbr: &str,
        residue: &Residue,
        source: PolymerizerError,
    ) -> Self {
        Self::NamedModification(
            name.to_owned(),
            abbr.to_owned(),
            residue.id(),
            residue.name.to_owned(),
            source.into(),
        )
    }

    fn bond(
        kind: &str,
        donor: &Residue,
        acceptor: &Residue,
        error_with: &str,
        source: PolymerizerError,
    ) -> Self {
        Self::Bond(
            kind.to_owned(),
            donor.id(),
            donor.name.to_owned(),
            acceptor.id(),
            acceptor.name.to_owned(),
            error_with.to_owned(),
            source.into(),
        )
    }
}
