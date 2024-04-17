//! An abstraction for building chemically validated polymers

pub mod atoms;
pub mod errors;
pub mod parsers;
pub mod polymerizer;
pub mod polymers;
#[cfg(test)]
mod testing_tools;

use std::num::NonZeroU32;

use polymerizer::Polymerizer;
use serde::Serialize;

// External Crate Imports
use ahash::{HashMap, HashSet};
use rust_decimal::Decimal;

// FIXME: Work on what's publicly exported / part of the API! — maybe create a prelude?
pub use atoms::atomic_database::AtomicDatabase;
pub use polymers::polymer_database::PolymerDatabase;

// FIXME: Blocks here need reordering!
// FIXME: Add exhaustive derives to everything!
// Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Display, Default, Serialize
// #[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize)]

// Core Data Types =====================================================================================================

// NOTE: For the types in this module, 'a lifetimes indicate references to the AtomicDatabase, whilst 'p lifetimes
// indicate references to the PolymerDatabase
#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
pub struct Polymer<'a, 'p> {
    // DESIGN: The `Polymerizer` struct keeps track of fields that are only useful when adding new components to the
    // current polymer. Abstracting those fields out into a `Polymerizer` helps to separate concerns, and makes it
    // clearer that the `residues`, `modifications`, and `bonds` fields are, by themselves, a complete description of
    // the `Polymer`s structure. A discarded alternative was two versions of the `Polymer` struct ("complete" and
    // "non-complete" versions), which would have introduced additional complexity to the user API, and the free group
    // index lost by discarding the `Polymerizer` would have needed to be reconstructed during fragmentation anyways
    // for quick group lookup. Having two structs just for the ID counter and database references didn't seem worth it.
    // This decision accepts the downside of increasing the size of the `Polymer` struct.
    #[serde(skip)]
    polymerizer: Polymerizer<'a, 'p>,
    // NOTE: Whilst the `Polymerizer` struct is defined elsewhere, the following fields are part of the core polymer
    // representation and are therefore defined in this file, with `impl`s added in separate modules as needed.
    residues: HashMap<ResidueId, Residue<'a, 'p>>,
    modifications: HashMap<ModificationId, ModificationInfo<'a, 'p>>,
    bonds: HashMap<BondId, BondInfo<'a, 'p>>,
}

// ---------------------------------------------------------------------------------------------------------------------

// FIXME: Pass through Display implementations for all Id newtypes!
// FIXME: Add tests that fail if `Default` is implemented for `Id`s!
// MISSING: No `Default` — should not be constructable by the user
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct ResidueId(Id);

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
struct Residue<'a, 'p> {
    abbr: &'p str,
    name: &'p str,
    composition: &'p ChemicalComposition<'a>,
    functional_groups: HashMap<FunctionalGroup<'p>, GroupState>,
    offset_modifications: HashSet<ModificationId>,
}

// MISSING: No `Default` — should not be constructable by the user
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct ModificationId(Id);

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
enum ModificationInfo<'a, 'p> {
    Named(NamedMod<'a, 'p>, ResidueId, &'p FunctionalGroup<'p>),
    Offset(Modification<OffsetMod<'a>>, ResidueId),
    Unlocalized(AnyModification<'a, 'p>),
}

// MISSING: No `Default` — should not be constructable by the user
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct BondId(Id);

// FIXME: Perhaps I should consider changing these `*Info` structs to have named fields? Is the donor -> acceptor order
// obvious enough for internal use? Users of `polychem` should never see this...
#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
struct BondInfo<'a, 'p>(ResidueId, Bond<'a, 'p>, ResidueId);

// ---------------------------------------------------------------------------------------------------------------------

// NOTE: This underlying `Id` type is just a synonym since it's private to this crate — there is no way for users to
// directly provide or modify `Id`s, so I don't need to worry about a newtype wrapper here
type Id = usize;

// MISSING: No `Default` — this *should* be user constructable, but there is no sensible default here
#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
pub struct ChemicalComposition<'a> {
    chemical_formula: Vec<(Element<'a>, Count)>,
    particle_offset: Option<(OffsetKind, Count, Particle<'a>)>,
}

// MISSING: No `Default` — should not be constructable by the user
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct FunctionalGroup<'p> {
    name: &'p str,
    location: &'p str,
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize)]
pub enum GroupState {
    #[default]
    Free,
    Modified(ModificationId),
    Donor(BondId),
    Acceptor(BondId),
}

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
struct NamedMod<'a, 'p> {
    abbr: &'p str,
    name: &'p str,
    lost: &'p ChemicalComposition<'a>,
    gained: &'p ChemicalComposition<'a>,
}

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
struct OffsetMod<'a> {
    kind: OffsetKind,
    composition: ChemicalComposition<'a>,
}

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
struct Modification<K> {
    multiplier: Count,
    kind: K,
}

pub type AnyModification<'a, 'p> = Modification<AnyMod<'a, 'p>>;

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
struct Bond<'a, 'p> {
    abbr: &'p str,
    name: &'p str,
    lost: &'p ChemicalComposition<'a>,
}

// ---------------------------------------------------------------------------------------------------------------------

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
struct Element<'a> {
    symbol: &'a str,
    name: &'a str,
    mass_number: Option<MassNumber>,
    isotopes: &'a HashMap<MassNumber, Isotope>,
}

// FIXME: Impl `Default` as Count(1)
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct Count(NonZeroU32);

// MISSING: No `Default` — this *should* be user constructable, but there is no sensible default here
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub enum OffsetKind {
    Add,
    Remove,
}

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
struct Particle<'a> {
    symbol: &'a str,
    name: &'a str,
    mass: &'a Mass,
    charge: &'a Charge,
}

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
enum AnyMod<'a, 'p> {
    Named(NamedMod<'a, 'p>),
    Offset(OffsetMod<'a>),
}

// ---------------------------------------------------------------------------------------------------------------------

// MISSING: No `Default` — this *should* be user constructable, but there is no sensible default here
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct MassNumber(u32);

#[derive(Clone, Eq, PartialEq, Debug, Serialize)]
struct Isotope {
    relative_mass: Mass,
    abundance: Option<Abundance>,
}

// MISSING: No `Default` — should not be constructable by the user
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct Mass(Decimal);

// MISSING: No `Default` — should not be constructable by the user
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct Charge(i64);

// ---------------------------------------------------------------------------------------------------------------------

// MISSING: No `Default` — should not be constructable by the user
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct Abundance(Decimal);

// MISSING: No `Default` — should not be constructable by the user
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct Mz(Decimal);

// =====================================================================================================================

pub trait Massive {
    fn monoisotopic_mass(&self) -> Mass;
    fn average_mass(&self) -> Mass;
}

pub trait Charged {
    fn charge(&self) -> Charge;
}

// FIXME: Not super sold on that trait name...
pub trait ChargedParticle: Massive + Charged {
    fn monoisotopic_mz(&self) -> Option<Mz> {
        let mass = self.monoisotopic_mass();
        let charge = self.charge().abs();
        mass.with_charge(charge)
    }

    fn average_mz(&self) -> Option<Mz> {
        let mass = self.average_mass();
        let charge = self.charge().abs();
        mass.with_charge(charge)
    }
}

// Blanket impls

macro_rules! massive_ref_impls {
    ($($ref_type:ty),+ $(,)?) => {
        $(
            impl<T: Massive> Massive for $ref_type {
                fn monoisotopic_mass(&self) -> Mass {
                    (**self).monoisotopic_mass()
                }

                fn average_mass(&self) -> Mass {
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

impl<T: Massive + Charged> ChargedParticle for T {}
