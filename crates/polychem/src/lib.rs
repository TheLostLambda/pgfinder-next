//! An abstraction for building chemically validated polymers

pub mod atoms;
pub mod errors;
pub mod parsers;
pub mod polymerizer;
pub mod polymers;
#[cfg(test)]
mod testing_tools;

use std::num::NonZeroU32;

use polymers::target::Index;
use serde::Serialize;

// External Crate Imports
use ahash::{HashMap, HashSet};
use rust_decimal::Decimal;

// FIXME: Work on what's publicly exported / part of the API! — maybe create a prelude?
pub use atoms::atomic_database::AtomicDatabase;
pub use polymers::polymer_database::PolymerDatabase;

// FIXME: Blocks here need reordering!

// NOTE: For the types in this module, 'a lifetimes indicate references to the AtomicDatabase, whilst 'p lifetimes
// indicate references to the PolymerDatabase
// FIXME: Add exhaustive derives to everything!
// Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Display, Default, Serialize
// #[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize)]
pub struct Polymer<'a, 'p> {
    atomic_db: &'a AtomicDatabase,
    polymer_db: &'p PolymerDatabase<'a>,
    next_id: Id,
    groups: Index<'p, HashMap<ResidueId, bool>>,
    residues: HashMap<ResidueId, Residue<'a, 'p>>,
    modifications: HashMap<ModificationId, ModificationInfo<'a, 'p>>,
    bonds: HashMap<BondId, BondInfo<'a, 'p>>,
}

enum ModificationInfo<'a, 'p> {
    Named(NamedMod<'a, 'p>, ResidueId, &'p FunctionalGroup<'p>),
    Offset(Modification<OffsetMod<'a>>, ResidueId),
    Unlocalized(AnyModification<'a, 'p>),
}

// FIXME: Perhaps I should consider changing these `*Info` structs to have named fields? Is the donor -> acceptor order
// obvious enough for internal use? Users of `polychem` should never see this...
struct BondInfo<'a, 'p>(ResidueId, Bond<'a, 'p>, ResidueId);

struct Residue<'a, 'p> {
    abbr: &'p str,
    name: &'p str,
    composition: &'p ChemicalComposition<'a>,
    functional_groups: HashMap<FunctionalGroup<'p>, GroupState>,
    offset_modifications: HashSet<ModificationId>,
}

// ---------------------------------------------------------------------------------------------------------------------

// NOTE: This underlying `Id` type is just a synonym since it's private to this crate — there is no way for users to
// directly provide or modify `Id`s, so I don't need to worry about a newtype wrapper here
type Id = usize;

// FIXME: Pass through Display implementations for all Id newtypes!
// MISSING: `*Id` types intentionally don't implement `Default`, since it should not be possible for users to construct
// their own — only `Polymer` should be capable of forging new `Id`s
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct ResidueId(Id);

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct ModificationId(Id);

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct BondId(Id);

pub struct ChemicalComposition<'a> {
    chemical_formula: Vec<(Element<'a>, Count)>,
    particle_offset: Option<(OffsetKind, Count, Particle<'a>)>,
}

// FIXME: Ensure all of the derives in this file derive as much as possible!
pub struct FunctionalGroup<'p> {
    name: &'p str,
    location: &'p str,
}

struct Modification<K> {
    multiplier: Count,
    kind: K,
}

// ---------------------------------------------------------------------------------------------------------------------

struct Element<'a> {
    symbol: &'a str,
    name: &'a str,
    mass_number: Option<MassNumber>,
    isotopes: &'a HashMap<MassNumber, Isotope>,
}

// FIXME: Impl `Default` as Count(1)
pub struct Count(NonZeroU32);

pub enum OffsetKind {
    Add,
    Remove,
}

struct Particle<'a> {
    symbol: &'a str,
    name: &'a str,
    mass: &'a Decimal,
    charge: &'a Charge,
}

pub enum GroupState {
    #[default]
    Free,
    Modified(ModificationId),
    Donor(BondId),
    Acceptor(BondId),
}

pub struct NamedMod<'a, 'p> {
    abbr: &'p str,
    name: &'p str,
    lost: &'p ChemicalComposition<'a>,
    gained: &'p ChemicalComposition<'a>,
}

pub struct OffsetMod<'a> {
    kind: OffsetKind,
    composition: ChemicalComposition<'a>,
}

// ---------------------------------------------------------------------------------------------------------------------

pub struct MassNumber(NonZeroU32);

struct Isotope {
    relative_mass: Decimal,
    abundance: Option<Decimal>,
}

pub struct Charge(i64);

pub struct Bond<'a, 'p> {
    abbr: &'p str,
    name: &'p str,
    lost: &'p ChemicalComposition<'a>,
}

// ---------------------------------------------------------------------------------------------------------------------

// FIXME: Better section naming!
// Convenience API? ====================================================================================================

pub type AnyModification<'a, 'p> = Modification<AnyMod<'a, 'p>>;

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
