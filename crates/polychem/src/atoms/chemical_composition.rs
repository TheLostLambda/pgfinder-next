use std::{
    fmt::{self, Display, Formatter},
    hash::{Hash, Hasher},
};

// External Crate Imports
use nom_miette::final_parser;

// Local Crate Imports
use crate::{
    parsers::chemical_composition::chemical_composition, AtomicDatabase, AverageMass, Charge,
    Charged, ChemicalComposition, Element, Mass, Massive, MonoisotopicMass, Result,
};

// Public API ==========================================================================================================

impl<'a> ChemicalComposition<'a> {
    pub fn new(db: &'a AtomicDatabase, formula: impl AsRef<str>) -> Result<Self> {
        let mut parser = final_parser(chemical_composition(db));
        parser(formula.as_ref()).map_err(|e| Box::new(e.into()))
    }
}

// Massive, Charged, and Mz Trait Implementations ======================================================================

impl Massive for ChemicalComposition<'_> {
    fn monoisotopic_mass(&self) -> MonoisotopicMass {
        self.mass(Element::monoisotopic_mass).into()
    }

    fn average_mass(&self) -> AverageMass {
        self.mass(Element::average_mass).into()
    }
}

impl Charged for ChemicalComposition<'_> {
    fn charge(&self) -> Charge {
        self.particle_offset
            .as_ref()
            .map(|&(offset_kind, count, ref particle)| {
                offset_kind.offset(count * particle.charge())
            })
            .unwrap_or(Charge(0))
    }
}

// Display and Hash Trait Implementations ==============================================================================

impl Display for ChemicalComposition<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        for &(ref element, count) in &self.chemical_formula {
            write!(f, "{element}{count}")?;
        }

        if let Some((offset_kind, count, ref particle)) = self.particle_offset {
            if self.chemical_formula.is_empty() {
                write!(f, "{count}{particle}")?;
            } else {
                write!(f, "{offset_kind}{count}{particle}")?;
            }
        }

        Ok(())
    }
}

impl Hash for ChemicalComposition<'_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // PERF: It's possible to do a bit better than this, but we're keeping it straightforward until we need speed
        self.to_string().hash(state);
    }
}

// Private Helper Methods ==============================================================================================

impl<'a> ChemicalComposition<'a> {
    fn mass<T: Into<Mass>>(&self, accessor: impl Fn(&Element<'a>) -> T) -> Mass {
        let element_masses = self
            .chemical_formula
            .iter()
            .map(|&(ref element, count)| count * accessor(element).into());

        let particle_masses = self
            .particle_offset
            .iter()
            .map(|&(offset_kind, count, ref particle)| offset_kind.offset(count * particle.mass()));

        element_masses.chain(particle_masses).sum()
    }
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::testing_tools::assert_miette_snapshot;

    use super::*;

    static DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);

    #[test]
    fn test_composition_errors() {
        let chemical_composition = |formula| ChemicalComposition::new(&DB, formula);
        // Looking up non-existant isotopes, elements, and particles
        assert_miette_snapshot!(chemical_composition("NH2[100Tc]O4"));
        assert_miette_snapshot!(chemical_composition("NH2[99Tc]YhO4"));
        assert_miette_snapshot!(chemical_composition("NH2[99Tc]O4-8m+2p"));
        // Starting a composition without an element or isotope
        assert_miette_snapshot!(chemical_composition("+H2O"));
        assert_miette_snapshot!(chemical_composition("-H2O"));
        assert_miette_snapshot!(chemical_composition("]H2O"));
        // Check counts are non-zero (no leading zeroes either!)
        assert_miette_snapshot!(chemical_composition("C3H0N4"));
        assert_miette_snapshot!(chemical_composition("C3H06N4"));
        // Ensure that particles are lowercase
        assert_miette_snapshot!(chemical_composition("H2O+P"));
        // Ensure that isotope expressions are valid
        assert_miette_snapshot!(chemical_composition("[H2O"));
        assert_miette_snapshot!(chemical_composition("[0H2O"));
        assert_miette_snapshot!(chemical_composition("[10H2O"));
        assert_miette_snapshot!(chemical_composition("[10]H2O"));
        assert_miette_snapshot!(chemical_composition("[37Cl"));
        // Check labels at the end of an input
        assert_miette_snapshot!(chemical_composition("[37Cl]5-"));
        assert_miette_snapshot!(chemical_composition("[37Cl]5+10"));
        // Check for partially valid input
        assert_miette_snapshot!(chemical_composition("[37Cl]52p"));
        assert_miette_snapshot!(chemical_composition("NH2[99Tc]O,4-2e+3p"));
        assert_miette_snapshot!(chemical_composition("eH2O"));
        // Check that multiple labels are reported for errors with different spans
        assert_miette_snapshot!(chemical_composition("2H2"));
        // Report elements lacking any natural abundance data
        assert_miette_snapshot!(chemical_composition("HTcN"));
        assert_miette_snapshot!(chemical_composition("C12H6PoS"));
    }

    #[test]
    fn composition_display() {
        let formulae = [
            "2p",
            "C11H10ON2",
            "C37H63N7O21+p",
            "Ca-2e",
            "Cr2O7+2e",
            "D2O",
            "H2O",
            "K-e",
            "NH2+2p",
            "NH2[99Tc]",
            "NH3+p",
            "Na-e",
            "OH+e",
            "[13C]11H10O[15N]2",
            "[2H]2O",
            "[37Cl]5-2p",
            "p",
        ];
        for formula in formulae {
            let composition = ChemicalComposition::new(&DB, formula).unwrap();
            assert_eq!(composition.to_string(), formula);
        }
    }

    #[test]
    fn composition_monoisotopic_mass() {
        // The masses here have been checked against https://mstools.epfl.ch/info/
        let water = ChemicalComposition::new(&DB, "H2O").unwrap();
        assert_eq!(
            water.monoisotopic_mass(),
            MonoisotopicMass(dec!(18.01056468403))
        );
        let trp_residue = ChemicalComposition::new(&DB, "C11H10ON2").unwrap();
        assert_eq!(
            trp_residue.monoisotopic_mass(),
            MonoisotopicMass(dec!(186.07931295073))
        );
        let trp_isotopes = ChemicalComposition::new(&DB, "[13C]11H10O[15N]2").unwrap();
        assert_eq!(
            trp_isotopes.monoisotopic_mass(),
            MonoisotopicMass(dec!(199.1102859254))
        );
        let gm_aeja = ChemicalComposition::new(&DB, "C37H63N7O21+p").unwrap();
        assert_eq!(
            gm_aeja.monoisotopic_mass(),
            MonoisotopicMass(dec!(942.414978539091))
        );

        // Testing with proton offsets for adducts (checked against https://www.unimod.org/modifications_list.php)
        let p2 = ChemicalComposition::new(&DB, "2p").unwrap();
        let ca2 = ChemicalComposition::new(&DB, "Ca-2e").unwrap();
        assert_eq!(
            ca2.monoisotopic_mass() - p2.monoisotopic_mass(),
            MonoisotopicMass(dec!(37.946940769939870))
        );
        let p1 = ChemicalComposition::new(&DB, "p").unwrap();
        let k1 = ChemicalComposition::new(&DB, "K-e").unwrap();
        assert_eq!(
            k1.monoisotopic_mass() - p1.monoisotopic_mass(),
            MonoisotopicMass(dec!(37.955881439869935))
        );
    }

    #[test]
    fn composition_average_mass() {
        // The masses here have been checked against https://mstools.epfl.ch/info/
        let water = ChemicalComposition::new(&DB, "H2O").unwrap();
        assert_eq!(
            water.average_mass(),
            AverageMass(dec!(18.01528643242983260))
        );
        let trp_residue = ChemicalComposition::new(&DB, "C11H10ON2").unwrap();
        assert_eq!(
            trp_residue.average_mass(),
            AverageMass(dec!(186.21031375185538640))
        );
        let trp_isotopes = ChemicalComposition::new(&DB, "[13C]11H10O[15N]2").unwrap();
        assert_eq!(
            trp_isotopes.average_mass(),
            AverageMass(dec!(199.11593344840605140))
        );
        let gm_aeja = ChemicalComposition::new(&DB, "C37H63N7O21+p").unwrap();
        assert_eq!(
            gm_aeja.average_mass(),
            AverageMass(dec!(942.93919804214360795))
        );

        // Testing with proton offsets for adducts (checked against https://www.unimod.org/modifications_list.php)
        let p2 = ChemicalComposition::new(&DB, "2p").unwrap();
        let ca2 = ChemicalComposition::new(&DB, "Ca-2e").unwrap();
        assert_eq!(
            ca2.average_mass() - p2.average_mass(),
            AverageMass(dec!(38.062372417957600))
        );
        let p1 = ChemicalComposition::new(&DB, "p").unwrap();
        let k1 = ChemicalComposition::new(&DB, "K-e").unwrap();
        assert_eq!(
            k1.average_mass() - p1.average_mass(),
            AverageMass(dec!(38.0904758635559412))
        );
    }

    #[test]
    fn composition_charges() {
        // Return charge 0 for compositions without particle offsets
        assert_eq!(
            ChemicalComposition::new(&DB, "Ca").unwrap().charge(),
            Charge(0)
        );
        // Get the charges for chemical formulae with particle offsets
        assert_eq!(
            ChemicalComposition::new(&DB, "Ca-2e").unwrap().charge(),
            Charge(2)
        );
        assert_eq!(
            ChemicalComposition::new(&DB, "Ca+2p").unwrap().charge(),
            Charge(2)
        );
        assert_eq!(
            ChemicalComposition::new(&DB, "Ca+p").unwrap().charge(),
            Charge(1)
        );
        assert_eq!(
            ChemicalComposition::new(&DB, "Ca-p").unwrap().charge(),
            Charge(-1)
        );
        assert_eq!(
            ChemicalComposition::new(&DB, "Ca+3e").unwrap().charge(),
            Charge(-3)
        );
        // Get the charges for standalone particle offsets
        assert_eq!(
            ChemicalComposition::new(&DB, "e").unwrap().charge(),
            Charge(-1)
        );
        assert_eq!(
            ChemicalComposition::new(&DB, "p").unwrap().charge(),
            Charge(1)
        );
        assert_eq!(
            ChemicalComposition::new(&DB, "3e").unwrap().charge(),
            Charge(-3)
        );
        assert_eq!(
            ChemicalComposition::new(&DB, "5p").unwrap().charge(),
            Charge(5)
        );
    }
}
