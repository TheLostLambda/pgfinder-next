//! An abstraction for building chemically validated polymers

use rust_decimal::Decimal;
use std::collections::HashMap;

pub struct Residue {
    id: Id,
    abbr: String,
    name: String,
    composition: ChemicalComposition,
    functional_groups: HashMap<Location, FunctionalGroup>,
    offset_modifications: Vec<Modification>,
}

type Id = usize;

struct ChemicalComposition {
    chemical_formula: Vec<(Element, u32)>,
    charged_particles: Vec<(OffsetKind, u32, Particle)>,
}

type Location = String;

struct FunctionalGroup {
    name: String,
    state: GroupState,
}

pub struct Modification {
    multiplier: u32,
    kind: ModificationKind,
}

struct Element {
    symbol: String,
    name: String,
    mass_number: Option<MassNumber>,
    isotopes: HashMap<MassNumber, Isotope>,
}

enum OffsetKind {
    Add,
    Remove,
}

struct Particle {
    symbol: String,
    name: String,
    mass: Decimal,
    charge: i32,
}

enum GroupState {
    Free,
    Modified(Modification),
    Donor(Bond),
    Acceptor,
}

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

struct Isotope {
    relative_mass: Decimal,
    abundance: Option<Decimal>,
}

struct Bond {
    kind: String,
    lost_mass: ChemicalComposition,
    acceptor: BondTarget,
}

struct BondTarget {
    residue: Id,
    group_location: Location,
}
