//! An abstraction for building chemically validated polymers

pub mod atoms;
pub mod polymerizer;
pub mod polymers;
#[cfg(test)]
mod testing_tools;

use atoms::chemical_composition_parser::CompositionError;
use polymerizer::PolymerizerError;
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
    offset_modifications: Vec<Modification<OffsetMod<'a>>>,
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

// FIXME: These fields are public because there are no internal invariants to uphold (it's just a straightforward tuple)
// and there are no read-only fields. It makes perfect sense for users to change either field. Be sure to apply this
// reasoning consistently to all of the other structs! Making more public where everything should be read-write! Don't
// forget to think about what might need to hold a mass-cache in the future!
// FIXME: In fact, damn, this could reasonably have a mass-cache... So I should probably re-private those fields...
#[derive(Clone, PartialEq, Eq, Debug, Serialize)]
pub struct Modification<K> {
    pub multiplier: Count,
    pub kind: K,
}

type SignedCount = i64;

// ---------------------------------------------------------------------------------------------------------------------

#[derive(Clone, PartialEq, Eq, Debug, Serialize)]
struct Element<'a> {
    symbol: &'a str,
    name: &'a str,
    mass_number: Option<MassNumber>,
    isotopes: &'a HashMap<MassNumber, Isotope>,
}

pub type Count = u32;

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

// FIXME: Move these "Any" types to their own section, after the main tree of data structures, since they are only used
// *outside* of this crate!
pub type AnyModification<'a, 'p> = Modification<AnyMod<'a, 'p>>;
#[derive(Clone, PartialEq, Eq, Debug, Serialize)]
pub enum AnyMod<'a, 'p> {
    Named(NamedMod<'a, 'p>),
    Offset(OffsetMod<'a>),
}

#[derive(Copy, Clone, PartialEq, Eq, Debug, Serialize)]
pub struct NamedMod<'a, 'p> {
    abbr: &'p str,
    name: &'p str,
    lost: &'p ChemicalComposition<'a>,
    gained: &'p ChemicalComposition<'a>,
}

#[derive(Clone, PartialEq, Eq, Debug, Serialize)]
pub struct OffsetMod<'a> {
    kind: OffsetKind,
    composition: ChemicalComposition<'a>,
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
    group: FunctionalGroup<'p>,
}

// ---------------------------------------------------------------------------------------------------------------------

pub trait Massive {
    fn monoisotopic_mass(&self) -> Decimal;
    fn average_mass(&self) -> Decimal;
}

pub trait Charged {
    fn charge(&self) -> Charge;
}

pub trait Mz: Massive + Charged {
    fn monoisotopic_mz(&self) -> Option<Decimal> {
        let charge = self.charge().abs();
        (charge != 0).then(|| self.monoisotopic_mass() / Decimal::from(charge))
    }

    fn average_mz(&self) -> Option<Decimal> {
        let charge = self.charge().abs();
        (charge != 0).then(|| self.average_mass() / Decimal::from(charge))
    }
}

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

    #[error("failed to apply the modification {0} ({1}) to residue {2} ({3})")]
    Modification(
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

    fn modification(name: &str, abbr: &str, residue: &Residue, source: PolymerizerError) -> Self {
        Self::Modification(
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
