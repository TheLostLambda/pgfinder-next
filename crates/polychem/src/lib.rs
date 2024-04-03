//! An abstraction for building chemically validated polymers

pub mod atoms;
pub mod parsers;
pub mod polymerizer;
pub mod polymers;
#[cfg(test)]
mod testing_tools;

use parsers::errors::CompositionError;
use polymerizer::PolymerizerError;
use polymers::errors::OffsetMultiplierError;
use serde::Serialize;

// External Crate Imports
use ahash::HashMap;
use miette::Diagnostic;
use rust_decimal::Decimal;
use thiserror::Error;

// FIXME: Work on what's publicly exported / part of the API!
pub use atoms::atomic_database::AtomicDatabase;
pub use polymers::polymer_database::PolymerDatabase;

// FIXME: Blocks here need reordering!

// NOTE: For the types in this module, 'a lifetimes indicate references to the AtomicDatabase, whilst 'p lifetimes
// indicate references to the PolymerDatabase
#[derive(Clone, PartialEq, Eq, Debug, Serialize)]
pub struct Residue<'a, 'p> {
    id: Id,
    abbr: &'p str,
    name: &'p str,
    composition: &'p ChemicalComposition<'a>,
    functional_groups: HashMap<FunctionalGroup<'p>, GroupState<'a, 'p>>,
    offset_modifications: HashMap<ChemicalComposition<'a>, OffsetMultiplier>,
}

// ---------------------------------------------------------------------------------------------------------------------

type Id = usize;

#[derive(Clone, PartialEq, Eq, Debug, Default, Serialize)]
pub struct ChemicalComposition<'a> {
    chemical_formula: Vec<(Element<'a>, Count)>,
    particle_offset: Option<(OffsetKind, Count, Particle<'a>)>,
}

// FIXME: Ensure all of the derives in this file derive as much as possible!
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize)]
pub struct FunctionalGroup<'p> {
    name: &'p str,
    location: &'p str,
}

// FIXME: Oh boy, please pick a better name...
// FIXME: Should this really be a tuple struct?...
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize)]
pub struct OffsetMultiplier(OffsetKind, Count);

#[derive(Clone, PartialEq, Eq, Debug, Serialize)]
pub struct Modification<K> {
    multiplier: Count,
    kind: K,
}

// ---------------------------------------------------------------------------------------------------------------------

#[derive(Clone, PartialEq, Eq, Debug, Serialize)]
struct Element<'a> {
    symbol: &'a str,
    name: &'a str,
    mass_number: Option<MassNumber>,
    isotopes: &'a HashMap<MassNumber, Isotope>,
}

// NOTE: Keep an eye on https://github.com/rust-lang/rust/issues/120257 for non-zero types
pub type Count = u32;
pub type SignedCount = i64;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize)]
pub enum OffsetKind {
    Add,
    Remove,
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize)]
struct Particle<'a> {
    symbol: &'a str,
    name: &'a str,
    mass: &'a Decimal,
    charge: &'a Charge,
}

#[derive(Copy, Clone, PartialEq, Eq, Debug, Default, Serialize)]
pub enum GroupState<'a, 'p> {
    #[default]
    Free,
    Modified(NamedMod<'a, 'p>),
    Donor(Bond<'a, 'p>),
    Acceptor,
}

#[derive(Copy, Clone, PartialEq, Eq, Debug, Serialize)]
pub struct NamedMod<'a, 'p> {
    abbr: &'p str,
    name: &'p str,
    lost: &'p ChemicalComposition<'a>,
    gained: &'p ChemicalComposition<'a>,
}

pub type OffsetMod<'a> = Offset<ChemicalComposition<'a>>;
pub type BorrowedOffsetMod<'a> = Offset<&'a ChemicalComposition<'a>>;

#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize)]
pub struct Offset<C> {
    kind: OffsetKind,
    composition: C,
}

// ---------------------------------------------------------------------------------------------------------------------

type MassNumber = u32;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize)]
struct Isotope {
    relative_mass: Decimal,
    abundance: Option<Decimal>,
}

pub type Charge = i64;

#[derive(Copy, Clone, PartialEq, Eq, Debug, Serialize)]
pub struct Bond<'a, 'p> {
    kind: &'p str,
    lost: &'p ChemicalComposition<'a>,
    acceptor: BondTarget<'p>,
}

// ---------------------------------------------------------------------------------------------------------------------

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize)]
pub struct BondTarget<'p> {
    residue: Id,
    // NOTE: Owned because some of these are sometimes constructed from `Target`s returned by the `TargetIndex`
    group: FunctionalGroup<'p>,
}

// FIXME: Better section naming!
// Convenience API? ====================================================================================================

pub type AnyModification<'a, 'p> = Modification<AnyMod<'a, 'p>>;

#[derive(Clone, PartialEq, Eq, Debug, Serialize)]
pub enum AnyMod<'a, 'p> {
    Named(NamedMod<'a, 'p>),
    Offset(OffsetMod<'a>),
}

// =====================================================================================================================

pub trait Massive {
    fn monoisotopic_mass(&self) -> Decimal;
    fn average_mass(&self) -> Decimal;
}

pub trait Charged {
    fn charge(&self) -> Charge;
}

pub trait Mz: Massive + Charged {
    fn monoisotopic_mz(&self) -> Option<Decimal> {
        let charge = Decimal::from(self.charge()).abs();
        (!charge.is_zero()).then(|| self.monoisotopic_mass() / charge)
    }

    fn average_mz(&self) -> Option<Decimal> {
        let charge = Decimal::from(self.charge()).abs();
        (!charge.is_zero()).then(|| self.average_mass() / charge)
    }
}

// Blanket impls

macro_rules! massive_ref_impls {
    ($($ref_type:ty),+ $(,)?) => {
        $(
            impl<T: Massive> Massive for $ref_type {
                fn monoisotopic_mass(&self) -> Decimal {
                    (**self).monoisotopic_mass()
                }

                fn average_mass(&self) -> Decimal {
                    (**self).average_mass()
                }
            }
        )+
    };
}

massive_ref_impls!(&T, &mut T, Box<T>);

macro_rules! charged_ref_impls {
    ($($ref_type:ty),+ $(,)?) => {
        $(
            impl<T: Charged> Charged for $ref_type {
                fn charge(&self) -> Charge {
                    (**self).charge()
                }
            }
        )+
    };
}

charged_ref_impls!(&T, &mut T, Box<T>);

impl<T: Massive + Charged> Mz for T {}

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
