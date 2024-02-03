//! An abstraction for building chemically validated polymers

pub mod chemical_database;
mod composition_parser;
#[cfg(test)]
mod testing_tools;

use chemical_database::ChemicalDatabase;
use composition_parser::CompositionError;
use nom_miette::final_parser;

// Standard Library Imports
use std::collections::HashMap;

// External Crate Imports
use itertools::Itertools;
use miette::{Diagnostic, Result};
use rust_decimal::{prelude::Zero, Decimal};
use thiserror::Error;

use self::composition_parser::chemical_composition;

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Residue {
    id: Id,
    abbr: String,
    name: String,
    composition: ChemicalComposition,
    functional_groups: HashMap<Location, FunctionalGroup>,
    offset_modifications: Vec<Modification>,
}

// =====================================================================================================================

type Id = usize;

// FIXME: Keep this public so that people can build mass calculators
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct ChemicalComposition {
    chemical_formula: Vec<(Element, Count)>,
    particle_offset: Option<(OffsetKind, Count, Particle)>,
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

// =====================================================================================================================

#[derive(Clone, PartialEq, Eq, Debug)]
struct Element {
    symbol: String,
    name: String,
    mass_number: Option<MassNumber>,
    isotopes: HashMap<MassNumber, Isotope>,
}

type Count = u32;

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
    charge: Charge,
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

// =====================================================================================================================

type MassNumber = u32;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
struct Isotope {
    relative_mass: Decimal,
    abundance: Option<Decimal>,
}

type Charge = i64;

#[derive(Clone, PartialEq, Eq, Debug)]
struct Bond {
    kind: String,
    lost_mass: ChemicalComposition,
    acceptor: BondTarget,
}

// =====================================================================================================================

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
struct BondTarget {
    residue: Id,
    group_location: Location,
}

// =====================================================================================================================

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

// FIXME: Should this really be public?
#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
pub enum ChemicalLookupError {
    #[error("the element {0:?} could not be found in the supplied chemical database")]
    Element(String),
    #[error("the isotope \"{0}-{1}\" could not be found in the supplied chemical database, though the following {2} isotopes were found: {3:?}")]
    Isotope(String, MassNumber, String, Vec<MassNumber>),
    #[error("the particle {0:?} could not be found in the supplied chemical database")]
    Particle(String),
    // FIXME: Unforuntately, this error probably doesn't belong here... All of the other errors can be
    // encountered at parse time, but this one is only triggered by a mass calculation...
    #[error("no natural abundance data could be found for {0} ({1}), though the following isotopes were found: {2:?}")]
    Abundance(String, String, Vec<MassNumber>),
}

// FIXME: Maybe there are too many layers of things being wrapped here!
// FIXME: Maybe just rename this to be `Error`?
#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
pub enum PolychemError {
    #[error(transparent)]
    #[diagnostic(transparent)]
    Composition(#[from] CompositionError),
    // FIXME: Oof, are these even different enough to warrant different errors?
    #[error("failed to fetch isotope abundances for monoisotopic mass calculation")]
    MonoisotopicMass(#[diagnostic_source] ChemicalLookupError),
    #[error("failed to fetch isotope abundances for average mass calculation")]
    AverageMass(#[diagnostic_source] ChemicalLookupError),
}

impl ChemicalComposition {
    // FIXME: If this isn't public API, drop the AsRef — if it is, then add it for `db`
    pub fn new(db: &ChemicalDatabase, formula: impl AsRef<str>) -> Result<Self, CompositionError> {
        let formula = formula.as_ref();
        final_parser(chemical_composition(db))(formula)
    }

    pub fn monoisotopic_mass(&self) -> Result<Decimal, PolychemError> {
        self.mass(Element::monoisotopic_mass)
    }

    pub fn average_mass(&self) -> Result<Decimal, PolychemError> {
        self.mass(Element::average_mass)
    }

    pub fn charge(&self) -> Charge {
        self.particle_offset
            .as_ref()
            .map(|(k, c, p)| {
                let sign: Charge = k.into();
                let c = Charge::from(*c);
                sign * c * p.charge
            })
            .unwrap_or_default()
    }

    fn mass(
        &self,
        accessor: impl Fn(&Element) -> Result<Decimal, PolychemError>,
    ) -> Result<Decimal, PolychemError> {
        // NOTE: Not using iterators makes using `?` possible, but might shut me out of `rayon` optimizations
        let mut mass = Decimal::zero();

        for (element, count) in &self.chemical_formula {
            mass += Decimal::from(*count) * accessor(element)?;
        }

        if let Some((offset_kind, count, particle)) = &self.particle_offset {
            mass += Decimal::from(offset_kind) * Decimal::from(*count) * particle.mass;
        }

        Ok(mass)
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
            Err(ChemicalLookupError::Isotope(
                symbol.to_owned(),
                mass_number,
                element.name.clone(),
                element.isotopes.keys().copied().sorted().collect(),
            ))
        }
    }

    fn monoisotopic_mass(&self) -> Result<Decimal, PolychemError> {
        if let Some(mass) = self.isotope_mass() {
            Ok(mass)
        } else {
            // SAFETY: The call to `.unwrap()` is safe here since `.isotope_abundances()` is guaranteed to yield at
            // least one isotope
            Ok(self
                .isotope_abundances()
                .map_err(PolychemError::MonoisotopicMass)?
                .max_by_key(|i| i.abundance)
                .unwrap()
                .relative_mass)
        }
    }

    fn average_mass(&self) -> Result<Decimal, PolychemError> {
        if let Some(mass) = self.isotope_mass() {
            Ok(mass)
        } else {
            // SAFETY: The call to `.unwrap()` is safe here since `.isotope_abundances()` is guaranteed to yield at
            // only isotopes containing natural abundance data
            Ok(self
                .isotope_abundances()
                .map_err(PolychemError::AverageMass)?
                // .wrap_err("failed to fetch isotope abundances for average mass calculation")?
                .map(|i| i.relative_mass * i.abundance.unwrap())
                .sum())
        }
    }

    fn isotope_mass(&self) -> Option<Decimal> {
        self.mass_number
            .and_then(|a| self.isotopes.get(&a))
            .map(|i| i.relative_mass)
    }

    fn isotope_abundances(&self) -> Result<impl Iterator<Item = &Isotope>, ChemicalLookupError> {
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
            ))
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

    use crate::testing_tools::assert_miette_snapshot;

    use super::{ChemicalComposition, ChemicalDatabase, Element, Particle};

    static DB: Lazy<ChemicalDatabase> = Lazy::new(|| {
        ChemicalDatabase::from_kdl("chemistry.kdl", include_str!("../chemistry.kdl")).unwrap()
    });

    #[test]
    fn new_particle() -> Result<()> {
        // Sucessfully lookup particles that exist
        assert_debug_snapshot!(Particle::new(&DB, "p")?);
        assert_debug_snapshot!(Particle::new(&DB, "e")?);
        // Fail to lookup particles that don't exist
        assert_miette_snapshot!(Particle::new(&DB, "m"));
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
        assert_miette_snapshot!(Element::new(&DB, "R"));
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
        assert_miette_snapshot!(Element::new_isotope(&DB, "R", 42));
        // Fail to lookup isotopes that don't exist
        assert_miette_snapshot!(Element::new_isotope(&DB, "C", 15));
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
        assert_miette_snapshot!(Element::new(&DB, "Tc")?.monoisotopic_mass());
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
        assert_miette_snapshot!(Element::new(&DB, "Po")?.average_mass());
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

    #[test]
    fn composition_monoisotopic_mass() -> Result<()> {
        // The masses here have been checked against https://mstools.epfl.ch/info/
        let water = ChemicalComposition::new(&DB, "H2O")?;
        assert_eq!(water.monoisotopic_mass()?, dec!(18.01056468403));
        let trp_residue = ChemicalComposition::new(&DB, "C11H10ON2")?;
        assert_eq!(trp_residue.monoisotopic_mass()?, dec!(186.07931295073));
        let trp_isotopes = ChemicalComposition::new(&DB, "[13C]11H10O[15N]2")?;
        assert_eq!(trp_isotopes.monoisotopic_mass()?, dec!(199.1102859254));
        let gm_aeja = ChemicalComposition::new(&DB, "C37H63N7O21+p")?;
        assert_eq!(gm_aeja.monoisotopic_mass()?, dec!(942.414978539091));

        // Testing with proton offsets for adducts (checked against https://www.unimod.org/modifications_list.php)
        let p2 = ChemicalComposition::new(&DB, "2p")?;
        let ca2 = ChemicalComposition::new(&DB, "Ca-2e")?;
        assert_eq!(
            ca2.monoisotopic_mass()? - p2.monoisotopic_mass()?,
            dec!(37.946940769939870)
        );
        let p1 = ChemicalComposition::new(&DB, "p")?;
        let k1 = ChemicalComposition::new(&DB, "K-e")?;
        assert_eq!(
            k1.monoisotopic_mass()? - p1.monoisotopic_mass()?,
            dec!(37.955881439869935)
        );
        Ok(())
    }

    #[test]
    fn composition_average_mass() -> Result<()> {
        // The masses here have been checked against https://mstools.epfl.ch/info/
        let water = ChemicalComposition::new(&DB, "H2O")?;
        assert_eq!(water.average_mass()?, dec!(18.01528643242983260));
        let trp_residue = ChemicalComposition::new(&DB, "C11H10ON2")?;
        assert_eq!(trp_residue.average_mass()?, dec!(186.21031375185538640));
        let trp_isotopes = ChemicalComposition::new(&DB, "[13C]11H10O[15N]2")?;
        assert_eq!(trp_isotopes.average_mass()?, dec!(199.11593344840605140));
        let gm_aeja = ChemicalComposition::new(&DB, "C37H63N7O21+p")?;
        assert_eq!(gm_aeja.average_mass()?, dec!(942.93919804214360795));

        // Testing with proton offsets for adducts (checked against https://www.unimod.org/modifications_list.php)
        let p2 = ChemicalComposition::new(&DB, "2p")?;
        let ca2 = ChemicalComposition::new(&DB, "Ca-2e")?;
        assert_eq!(
            ca2.average_mass()? - p2.average_mass()?,
            dec!(38.062372417957600)
        );
        let p1 = ChemicalComposition::new(&DB, "p")?;
        let k1 = ChemicalComposition::new(&DB, "K-e")?;
        assert_eq!(
            k1.average_mass()? - p1.average_mass()?,
            dec!(38.0904758635559412)
        );
        Ok(())
    }

    #[test]
    fn composition_charges() -> Result<()> {
        // Return charge 0 for compositions without particle offsets
        assert_eq!(ChemicalComposition::new(&DB, "Ca")?.charge(), 0);
        // Get the charges for chemical formulae with particle offsets
        assert_eq!(ChemicalComposition::new(&DB, "Ca-2e")?.charge(), 2);
        assert_eq!(ChemicalComposition::new(&DB, "Ca+2p")?.charge(), 2);
        assert_eq!(ChemicalComposition::new(&DB, "Ca+p")?.charge(), 1);
        assert_eq!(ChemicalComposition::new(&DB, "Ca-p")?.charge(), -1);
        assert_eq!(ChemicalComposition::new(&DB, "Ca+3e")?.charge(), -3);
        // Get the charges for standalone particle offsets
        assert_eq!(ChemicalComposition::new(&DB, "e")?.charge(), -1);
        assert_eq!(ChemicalComposition::new(&DB, "p")?.charge(), 1);
        assert_eq!(ChemicalComposition::new(&DB, "3e")?.charge(), -3);
        assert_eq!(ChemicalComposition::new(&DB, "5p")?.charge(), 5);
        Ok(())
    }
}
