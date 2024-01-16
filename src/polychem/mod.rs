//! An abstraction for building chemically validated polymers

pub mod chemical_database;
pub use chemical_database::*;
use thiserror::Error;

// Standard Library Imports
use std::collections::HashMap;

// External Crate Imports
use miette::{Diagnostic, Result};
use rust_decimal::Decimal;

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

#[derive(Error, Diagnostic, PartialEq, Eq, Debug)]
enum ChemicalLookupError {
    #[error("the element {0:?} could not found in the supplied chemical database")]
    Element(String),
    #[error("the isotope \"{0}-{1}\" could not found in the supplied chemical database")]
    Isotope(String, MassNumber),
    #[error("the particle {0:?} could not found in the supplied chemical database")]
    Particle(String),
}

impl Element {
    pub fn new(db: &ChemicalDatabase, symbol: impl AsRef<str>) -> Result<Self> {
        let symbol = symbol.as_ref();
        db.elements
            .get(symbol)
            .map(|p| p.clone())
            .ok_or_else(|| ChemicalLookupError::Element(symbol.to_owned()).into())
    }

    pub fn new_isotope(
        db: &ChemicalDatabase,
        symbol: impl AsRef<str>,
        mass_number: MassNumber,
    ) -> Result<Self> {
        let symbol = symbol.as_ref();
        let element = Self::new(db, symbol)?;
        if element.isotopes.contains_key(&mass_number) {
            Ok(Self {
                mass_number: Some(mass_number),
                ..element
            })
        } else {
            Err(ChemicalLookupError::Isotope(symbol.to_owned(), mass_number).into())
        }
    }

    pub fn monoisotopic_mass(&self) -> Decimal {
        // TODO: For dealing with explicit isotopes here, just filter the iterator to select only that isotope
        if let Some((_, i)) = self.isotopes.iter().max_by_key(|(_, i)| i.abundance) {
            i.relative_mass
        } else {
            unreachable!("validation of the chemical database ensures that all elements contain at least one isotope")
        }
    }
}

impl Particle {
    pub fn new(db: &ChemicalDatabase, symbol: impl AsRef<str>) -> Result<Self> {
        let symbol = symbol.as_ref();
        db.particles
            .get(symbol)
            .map(|p| p.clone())
            .ok_or_else(|| ChemicalLookupError::Particle(symbol.to_owned()).into())
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use miette::Result;
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::polychem::ChemicalLookupError;

    use super::{ChemicalDatabase, Element, Particle};

    static DB: Lazy<ChemicalDatabase> = Lazy::new(|| {
        ChemicalDatabase::from_kdl("chemistry.kdl", include_str!("chemistry.kdl")).unwrap()
    });

    #[test]
    fn new_particle() -> Result<()> {
        // Sucessfully lookup particles that exist
        assert_debug_snapshot!(Particle::new(&DB, "p")?);
        assert_debug_snapshot!(Particle::new(&DB, "e")?);
        // Fail to lookup particles that don't exist
        let m = Particle::new(&DB, "m");
        assert_eq!(
            m.unwrap_err()
                .downcast_ref::<ChemicalLookupError>()
                .unwrap(),
            &ChemicalLookupError::Particle("m".to_string())
        );
        Ok(())
    }

    #[test]
    fn new_element() -> Result<()> {
        // Sucessfully lookup elements that exist
        let Element {
            symbol,
            name,
            mass_number,
            isotopes,
        } = Element::new(&DB, "C")?;
        assert_eq!(symbol, "C");
        assert_eq!(name, "Carbon");
        assert_eq!(mass_number, None);
        let mut isotopes: Vec<_> = isotopes.into_iter().collect();
        isotopes.sort_unstable_by(|(a, _), (b, _)| a.cmp(b));
        assert_debug_snapshot!(isotopes);
        // Fail to lookup elements that don't exist
        let m = Element::new(&DB, "R");
        assert_eq!(
            m.unwrap_err()
                .downcast_ref::<ChemicalLookupError>()
                .unwrap(),
            &ChemicalLookupError::Element("R".to_string())
        );
        Ok(())
    }

    #[test]
    fn new_isotope() -> Result<()> {
        // Sucessfully lookup isotopes that exist
        let Element {
            symbol,
            name,
            mass_number,
            isotopes,
        } = Element::new_isotope(&DB, "C", 13)?;
        assert_eq!(symbol, "C");
        assert_eq!(name, "Carbon");
        assert_eq!(mass_number, Some(13));
        let mut isotopes: Vec<_> = isotopes.into_iter().collect();
        isotopes.sort_unstable_by(|(a, _), (b, _)| a.cmp(b));
        assert_debug_snapshot!(isotopes);
        // Fail to lookup isotopes for elements that don't exist
        let m = Element::new_isotope(&DB, "R", 42);
        assert_eq!(
            m.unwrap_err()
                .downcast_ref::<ChemicalLookupError>()
                .unwrap(),
            &ChemicalLookupError::Element("R".to_string())
        );
        // Fail to lookup isotopes that don't exist
        let m = Element::new_isotope(&DB, "C", 15);
        assert_eq!(
            m.unwrap_err()
                .downcast_ref::<ChemicalLookupError>()
                .unwrap(),
            &ChemicalLookupError::Isotope("C".to_string(), 15)
        );
        Ok(())
    }

    #[test]
    fn element_monoisotopic_mass() -> Result<()> {
        let c = Element::new(&DB, "C")?;
        assert_eq!(c.monoisotopic_mass(), dec!(12));
        Ok(())
    }
}
