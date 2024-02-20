//! An abstraction for building chemically validated polymers

pub mod atoms;
pub mod polymerizer;
pub mod polymers;
#[cfg(test)]
mod testing_tools;

use atoms::{chemical_composition::CompositionError, AtomicLookupError};
use serde::{Deserialize, Serialize};

// Standard Library Imports
use std::collections::HashMap;

// External Crate Imports
use miette::Diagnostic;
use rust_decimal::Decimal;
use thiserror::Error;

// FIXME: Blocks here need reordering!

// FIXME: A more intense refactor, but things that don't change for residues, like abbr, name, composition, etc, should
// be stored as references to the chemical databases. Otherwise, when creating new residues, I'm doing a *lot* of
// copying that I really shouldn't need to do... I should really go through all of these types and use references for
// anything that's "static" / just comes from a config file
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct Residue {
    id: Id,
    abbr: String,
    name: String,
    composition: ChemicalComposition,
    functional_groups: HashMap<FunctionalGroup, GroupState>,
    offset_modifications: Vec<Modification>,
}

// ---------------------------------------------------------------------------------------------------------------------

type Id = usize;

#[derive(Clone, PartialEq, Eq, Debug, Default, Serialize, Deserialize)]
pub struct ChemicalComposition {
    chemical_formula: Vec<(Element, Count)>,
    particle_offset: Option<(OffsetKind, Count, Particle)>,
}

type Location = String;

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
pub struct FunctionalGroup {
    name: String,
    location: String,
}

#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct Modification {
    multiplier: Count,
    kind: ModificationKind,
}

// ---------------------------------------------------------------------------------------------------------------------

#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
struct Element {
    symbol: String,
    name: String,
    mass_number: Option<MassNumber>,
    isotopes: HashMap<MassNumber, Isotope>,
}

type Count = u32;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
enum OffsetKind {
    Add,
    Remove,
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
struct Particle {
    symbol: String,
    name: String,
    mass: Decimal,
    charge: Charge,
}

#[derive(Clone, PartialEq, Eq, Debug, Default, Serialize, Deserialize)]
enum GroupState {
    #[default]
    Free,
    Modified(Modification),
    Donor(Bond),
    Acceptor,
}

#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
enum ModificationKind {
    Predefined {
        abbr: String,
        name: String,
        lost: ChemicalComposition,
        gained: ChemicalComposition,
    },
    ChemicalOffset {
        kind: OffsetKind,
        composition: ChemicalComposition,
    },
}

// ---------------------------------------------------------------------------------------------------------------------

type MassNumber = u32;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
struct Isotope {
    relative_mass: Decimal,
    abundance: Option<Decimal>,
}

type Charge = i64;

#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
struct Bond {
    kind: String,
    lost: ChemicalComposition,
    acceptor: BondTarget,
}

// ---------------------------------------------------------------------------------------------------------------------

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
struct BondTarget {
    residue: Id,
    group_location: Location,
}

// ---------------------------------------------------------------------------------------------------------------------

// FIXME: These OffsetKind impls probably need a better home?
impl From<&OffsetKind> for Decimal {
    fn from(value: &OffsetKind) -> Self {
        Charge::from(value).into()
    }
}

impl From<&OffsetKind> for Charge {
    fn from(value: &OffsetKind) -> Self {
        match value {
            OffsetKind::Add => 1,
            OffsetKind::Remove => -1,
        }
    }
}

#[derive(Debug, Diagnostic, Clone, Error)]
#[error(transparent)]
#[diagnostic(transparent)]
pub struct Error(PolychemError);

impl<E: Into<PolychemError>> From<E> for Error {
    fn from(value: E) -> Self {
        Self(value.into())
    }
}

pub type Result<T> = std::result::Result<T, Error>;

// FIXME: Maybe there are too many layers of things being wrapped here!
// FIXME: Maybe just rename this to be `Error`?
// FIXME: Check all of the errors returned from public API are wrapped in this!
#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
enum PolychemError {
    #[error(transparent)]
    #[diagnostic(transparent)]
    Composition(#[from] CompositionError),

    // FIXME: Oof, are these even different enough to warrant different errors?
    #[error("failed to fetch isotope abundances for monoisotopic mass calculation")]
    MonoisotopicMass(
        #[source]
        #[diagnostic_source]
        AtomicLookupError,
    ),

    #[error("failed to fetch isotope abundances for average mass calculation")]
    AverageMass(
        #[source]
        #[diagnostic_source]
        AtomicLookupError,
    ),
}
