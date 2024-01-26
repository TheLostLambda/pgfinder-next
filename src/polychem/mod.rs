//! An abstraction for building chemically validated polymers

mod chemical_database;
mod composition_parser;
#[cfg(test)]
mod testing_tools;

use chemical_database::ChemicalDatabase;
use rust_decimal_macros::dec;

// Standard Library Imports
use std::collections::HashMap;

// External Crate Imports
use itertools::Itertools;
use miette::{Context, Diagnostic, Result};
use rust_decimal::{prelude::FromPrimitive, Decimal};
use thiserror::Error;

use self::composition_parser::{chemical_composition, final_parser, CompositionError};

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

// FIXME: Keep this public so that people can build mass calculators
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct ChemicalComposition {
    chemical_formula: Vec<(Element, Count)>,
    charged_particles: Vec<(OffsetKind, Count, Particle)>,
}

type Location = String;

#[derive(Clone, PartialEq, Eq, Debug)]
struct FunctionalGroup {
    name: String,
    state: GroupState,
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Modification {
    multiplier: Count,
    kind: ModificationKind,
}

type Count = u32;

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

#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
enum ChemicalLookupError {
    #[error("the element {0:?} could not found in the supplied chemical database")]
    Element(String),
    #[error("the isotope \"{0}-{1}\" could not found in the supplied chemical database")]
    Isotope(String, MassNumber),
    #[error("the particle {0:?} could not found in the supplied chemical database")]
    Particle(String),
    // FIXME: Unforuntately, this error probably doesn't belong here... All of the other errors can be
    // encountered at parse time, but this one is only triggered by a mass calculation...
    #[error("no natural abundance data could be found for {0} ({1}), though the following isotopes were found: {2:?}")]
    Abundance(String, String, Vec<MassNumber>),
}

impl ChemicalComposition {
    // FIXME: If this isn't public API, drop the AsRef — if it is, then add it for `db`
    fn new(db: &ChemicalDatabase, formula: impl AsRef<str>) -> Result<Self, CompositionError> {
        let formula = formula.as_ref();
        final_parser(chemical_composition(db))(formula)
    }

    fn monoisotopic_mass(&self) -> Result<Decimal> {
        self.mass(Element::monoisotopic_mass)
    }

    fn average_mass(&self) -> Result<Decimal> {
        self.mass(Element::average_mass)
    }

    // FIXME: No no no no! Need to properly handle errors here! No `.unwrap()`!
    fn mass(&self, accessor: impl Fn(&Element) -> Result<Decimal>) -> Result<Decimal> {
        let atom_masses: Decimal = self
            .chemical_formula
            .iter()
            .map(|(e, c)| Decimal::from_u32(*c).unwrap() * accessor(e).unwrap())
            .sum();
        let particle_masses: Decimal = self
            .charged_particles
            .iter()
            .map(|(k, c, p)| {
                // FIXME: Probably refactor this out...
                let sign = match k {
                    OffsetKind::Add => dec!(1),
                    OffsetKind::Remove => dec!(-1),
                };
                let c = Decimal::from_u32(*c).unwrap();
                sign * c * p.mass
            })
            .sum();
        Ok(atom_masses + particle_masses)
    }
}

impl Element {
    fn new(db: &ChemicalDatabase, symbol: impl AsRef<str>) -> Result<Self, ChemicalLookupError> {
        let symbol = symbol.as_ref();
        db.elements
            .get(symbol)
            .cloned()
            .ok_or_else(|| ChemicalLookupError::Element(symbol.to_owned()))
    }

    fn new_isotope(
        db: &ChemicalDatabase,
        symbol: impl AsRef<str>,
        mass_number: MassNumber,
    ) -> Result<Self, ChemicalLookupError> {
        let symbol = symbol.as_ref();
        let element = Self::new(db, symbol)?;
        if element.isotopes.contains_key(&mass_number) {
            Ok(Self {
                mass_number: Some(mass_number),
                ..element
            })
        } else {
            Err(ChemicalLookupError::Isotope(symbol.to_owned(), mass_number))
        }
    }

    fn monoisotopic_mass(&self) -> Result<Decimal> {
        if let Some(mass) = self.isotope_mass() {
            Ok(mass)
        } else {
            // SAFETY: The call to `.unwrap()` is safe here since `.isotope_abundances()` is guaranteed to yield at
            // least one isotope
            Ok(self
                .isotope_abundances()
                .wrap_err("failed to fetch isotope abundances for monoisotopic mass calculation")?
                .max_by_key(|i| i.abundance)
                .unwrap()
                .relative_mass)
        }
    }

    fn average_mass(&self) -> Result<Decimal> {
        if let Some(mass) = self.isotope_mass() {
            Ok(mass)
        } else {
            // SAFETY: The call to `.unwrap()` is safe here since `.isotope_abundances()` is guaranteed to yield at
            // only isotopes containing natural abundance data
            Ok(self
                .isotope_abundances()
                .wrap_err("failed to fetch isotope abundances for average mass calculation")?
                .map(|i| i.relative_mass * i.abundance.unwrap())
                .sum())
        }
    }

    fn isotope_mass(&self) -> Option<Decimal> {
        self.mass_number
            .and_then(|a| self.isotopes.get(&a))
            .map(|i| i.relative_mass)
    }

    fn isotope_abundances(&self) -> Result<impl Iterator<Item = &Isotope>> {
        let mut isotopes_with_abundances = self
            .isotopes
            .values()
            .filter(|i| i.abundance.is_some())
            .peekable();
        if isotopes_with_abundances.peek().is_some() {
            Ok(isotopes_with_abundances)
        } else {
            Err(ChemicalLookupError::Abundance(
                self.name.clone(),
                self.symbol.clone(),
                self.isotopes.keys().copied().sorted().collect(),
            )
            .into())
        }
    }
}

impl Particle {
    fn new(db: &ChemicalDatabase, symbol: impl AsRef<str>) -> Result<Self, ChemicalLookupError> {
        let symbol = symbol.as_ref();
        db.particles
            .get(symbol)
            .cloned()
            .ok_or_else(|| ChemicalLookupError::Particle(symbol.to_owned()))
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use miette::Result;
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::polychem::ChemicalLookupError;

    use super::{ChemicalComposition, ChemicalDatabase, Element, Particle};

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
            m.unwrap_err(),
            ChemicalLookupError::Particle("m".to_string())
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
            m.unwrap_err(),
            ChemicalLookupError::Element("R".to_string())
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
            m.unwrap_err(),
            ChemicalLookupError::Element("R".to_string())
        );
        // Fail to lookup isotopes that don't exist
        let m = Element::new_isotope(&DB, "C", 15);
        assert_eq!(
            m.unwrap_err(),
            ChemicalLookupError::Isotope("C".to_string(), 15)
        );
        Ok(())
    }

    #[test]
    fn element_monoisotopic_mass() -> Result<()> {
        // Successfully calculate the monoisotopic mass of elements with natural abundances
        let c = Element::new(&DB, "C")?.monoisotopic_mass();
        assert!(c.is_ok());
        assert_eq!(c.unwrap(), dec!(12));
        let mg = Element::new(&DB, "Mg")?.monoisotopic_mass();
        assert!(mg.is_ok());
        assert_eq!(mg.unwrap(), dec!(23.985041697));
        let mo = Element::new(&DB, "Mo")?.monoisotopic_mass();
        assert!(mo.is_ok());
        assert_eq!(mo.unwrap(), dec!(97.90540482));
        // Fail to calculate the monoisotopic mass of elements without natural abundances
        let tc = Element::new(&DB, "Tc")?.monoisotopic_mass();
        assert!(tc.is_err());
        assert_debug_snapshot!(tc.unwrap_err());
        Ok(())
    }

    #[test]
    fn element_average_mass() -> Result<()> {
        // Successfully calculate the average mass of elements with natural abundances
        let c = Element::new(&DB, "C")?.average_mass();
        assert!(c.is_ok());
        assert_eq!(c.unwrap(), dec!(12.010735896735249));
        let mg = Element::new(&DB, "Mg")?.average_mass();
        assert!(mg.is_ok());
        assert_eq!(mg.unwrap(), dec!(24.3050516198371));
        let mo = Element::new(&DB, "Mo")?.average_mass();
        assert!(mo.is_ok());
        assert_eq!(mo.unwrap(), dec!(95.959788541188));
        // Fail to calculate the monoisotopic mass of elements without natural abundances
        let po = Element::new(&DB, "Po")?.average_mass();
        assert!(po.is_err());
        assert_debug_snapshot!(po.unwrap_err());
        Ok(())
    }

    #[test]
    fn isotope_masses() -> Result<()> {
        // Get masses for an element with natural abundances
        let c13_mono = Element::new_isotope(&DB, "C", 13)?.monoisotopic_mass();
        assert!(c13_mono.is_ok());
        assert_eq!(c13_mono.unwrap(), dec!(13.00335483507));
        let c13_avg = Element::new_isotope(&DB, "C", 13)?.average_mass();
        assert!(c13_avg.is_ok());
        assert_eq!(c13_avg.unwrap(), dec!(13.00335483507));
        // Get masses for an element without natural abundances
        let tc99_mono = Element::new_isotope(&DB, "Tc", 99)?.monoisotopic_mass();
        assert!(tc99_mono.is_ok());
        assert_eq!(tc99_mono.unwrap(), dec!(98.9062508));
        let tc99_avg = Element::new_isotope(&DB, "Tc", 99)?.average_mass();
        assert!(tc99_avg.is_ok());
        assert_eq!(tc99_avg.unwrap(), dec!(98.9062508));
        Ok(())
    }

    // FIXME: Messy! Give this a good second pass! Needs isotopes and particles tested!
    #[test]
    fn compound_masses() -> Result<()> {
        // The masses here have been checked against https://mstools.epfl.ch/info/
        let water = ChemicalComposition::new(&DB, "H2O")?;
        assert_eq!(water.monoisotopic_mass()?, dec!(18.01056468403));
        assert_eq!(water.average_mass()?, dec!(18.01528643242983260));
        let trp_residue = ChemicalComposition::new(&DB, "C11H10ON2")?;
        assert_eq!(trp_residue.monoisotopic_mass()?, dec!(186.07931295073));
        assert_eq!(trp_residue.average_mass()?, dec!(186.21031375185538640));
        // Testing with proton offsets for adducts
        let ca = ChemicalComposition::new(&DB, "Ca-2p")?;
        assert_eq!(ca.monoisotopic_mass()?, dec!(37.948037929758));
        assert_eq!(ca.average_mass()?, dec!(38.06346957777573));
        let k = ChemicalComposition::new(&DB, "K-p")?;
        assert_eq!(k.monoisotopic_mass()?, dec!(37.956430019779));
        assert_eq!(k.average_mass()?, dec!(38.0910244434650062));
        Ok(())
    }
}
