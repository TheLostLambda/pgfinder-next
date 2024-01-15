//! An abstraction for building chemically validated polymers
pub mod chemical_database;
pub use chemical_database::*;

use rust_decimal::Decimal;
use std::collections::HashMap;

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Residue {
    id: Id,
    abbr: String,
    name: String,
    composition: ChemicalComposition,
    functional_groups: HashMap<Location, FunctionalGroup>,
    offset_modifications: Vec<Modification>,
}

type Id = usize;

#[derive(Clone, PartialEq, Eq, Debug)]
struct ChemicalComposition {
    chemical_formula: Vec<(Element, u32)>,
    charged_particles: Vec<(OffsetKind, u32, Particle)>,
}

type Location = String;

#[derive(Clone, PartialEq, Eq, Debug)]
struct FunctionalGroup {
    name: String,
    state: GroupState,
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Modification {
    multiplier: u32,
    kind: ModificationKind,
}

#[derive(Clone, PartialEq, Eq, Debug)]
struct Element {
    symbol: String,
    name: String,
    mass_number: Option<MassNumber>,
    isotopes: HashMap<MassNumber, Isotope>,
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
enum OffsetKind {
    Add,
    Remove,
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
struct Particle {
    symbol: String,
    name: String,
    mass: Decimal,
    charge: i32,
}

#[derive(Clone, PartialEq, Eq, Debug, Default)]
enum GroupState {
    #[default]
    Free,
    Modified(Modification),
    Donor(Bond),
    Acceptor,
}

#[derive(Clone, PartialEq, Eq, Debug)]
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

type MassNumber = u32;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
struct Isotope {
    relative_mass: Decimal,
    abundance: Option<Decimal>,
}

#[derive(Clone, PartialEq, Eq, Debug)]
struct Bond {
    kind: String,
    lost_mass: ChemicalComposition,
    acceptor: BondTarget,
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
struct BondTarget {
    residue: Id,
    group_location: Location,
}
