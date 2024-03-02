//! An abstraction for building chemically validated polymers

pub mod atoms;
pub mod polymerizer;
pub mod polymers;
#[cfg(test)]
mod testing_tools;

use atoms::chemical_composition::CompositionError;
use serde::Serialize;

// Standard Library Imports
use std::collections::HashMap;

// External Crate Imports
use miette::Diagnostic;
use rust_decimal::Decimal;
use thiserror::Error;

// FIXME: Blocks here need reordering!

// NOTE: For the types in this module, 'a lifetimes indicate references to the AtomicDatabase, whilst 'p lifetimes
// indicate references to the PolymerDatabase
#[derive(Clone, PartialEq, Eq, Debug, Serialize)]
pub struct Residue<'a, 'p> {
    id: Id,
    abbr: &'p str,
    name: &'p str,
    composition: &'p ChemicalComposition<'a>,
    functional_groups: HashMap<&'p FunctionalGroup, GroupState<'a, 'p>>,
    offset_modifications: Vec<Modification<OffsetMod<'a>>>,
}

// ---------------------------------------------------------------------------------------------------------------------

type Id = usize;

#[derive(Clone, PartialEq, Eq, Debug, Default, Serialize)]
pub struct ChemicalComposition<'a> {
    chemical_formula: Vec<(Element<'a>, Count)>,
    particle_offset: Option<(OffsetKind, Count, Particle<'a>)>,
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize)]
pub struct FunctionalGroup {
    name: String,
    location: String,
}

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

type Count = u32;

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

#[derive(Clone, PartialEq, Eq, Debug, Default, Serialize)]
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

#[derive(Clone, PartialEq, Eq, Debug, Serialize)]
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

type Charge = i64;

#[derive(Clone, PartialEq, Eq, Debug, Serialize)]
pub struct Bond<'a, 'p> {
    kind: &'p str,
    lost: &'p ChemicalComposition<'a>,
    acceptor: BondTarget<'p>,
}

// ---------------------------------------------------------------------------------------------------------------------

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize)]
struct BondTarget<'p> {
    residue: Id,
    group: &'p FunctionalGroup,
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

// FIXME: These OffsetKind impls probably need a better home?
impl From<OffsetKind> for Decimal {
    fn from(value: OffsetKind) -> Self {
        Charge::from(value).into()
    }
}

impl From<OffsetKind> for Charge {
    fn from(value: OffsetKind) -> Self {
        match value {
            OffsetKind::Add => 1,
            OffsetKind::Remove => -1,
        }
    }
}

#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
#[error(transparent)]
#[diagnostic(transparent)]
pub struct Error(PolychemError);

impl<E: Into<PolychemError>> From<E> for Error {
    fn from(value: E) -> Self {
        Self(value.into())
    }
}

pub type Result<T, E = Error> = std::result::Result<T, E>;

// FIXME: Maybe there are too many layers of things being wrapped here!
// FIXME: Maybe just rename this to be `Error`?
// FIXME: Check all of the errors returned from public API are wrapped in this!
#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
enum PolychemError {
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
    GroupLookup(FunctionalGroup, String, String),

    // FIXME: Destroy me
    #[error(
        "the functional group {0} of {1}-{2} was already {3}, but must be free for modification"
    )]
    ModificationGroupOccupied(FunctionalGroup, String, Id, String),

    // FIXME: Destroy me
    #[error(
        "the functional group {0} of {1}-{2} was already {3}, but must be free for bond formation"
    )]
    BondGroupOccupied(FunctionalGroup, String, Id, String),
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

    fn group_lookup(functional_group: &FunctionalGroup, name: &str, abbr: &str) -> Self {
        Self::GroupLookup(functional_group.clone(), name.to_owned(), abbr.to_owned())
    }

    // FIXME: Destroy me
    fn modification_group_occupied(group: &FunctionalGroup, residue: &Residue) -> Self {
        Self::ModificationGroupOccupied(
            group.clone(),
            residue.name.to_owned(),
            residue.id,
            residue.group_state(group).unwrap().to_string(),
        )
    }

    // FIXME: Destroy me
    fn bond_group_occupied(group: &FunctionalGroup, residue: &Residue) -> Self {
        Self::BondGroupOccupied(
            group.clone(),
            residue.name.to_owned(),
            residue.id,
            residue.group_state(group).unwrap().to_string(),
        )
    }
}
