// External Crate Imports
use nom_miette::final_parser;
use rust_decimal::Decimal;

// Local Crate Imports
use super::{atomic_database::AtomicDatabase, chemical_composition_parser::chemical_composition};
use crate::{Charge, Charged, ChemicalComposition, Element, Massive, Mz, Result};

// Public API ==========================================================================================================

impl<'a> ChemicalComposition<'a> {
    pub fn new(db: &'a AtomicDatabase, formula: impl AsRef<str>) -> Result<Self> {
        let mut parser = final_parser(chemical_composition(db));
        parser(formula.as_ref()).map_err(|e| Box::new(e.into()))
    }
}

// Massive, Charged, and Mz Trait Implementations ======================================================================

impl Massive for ChemicalComposition<'_> {
    fn monoisotopic_mass(&self) -> Decimal {
        self.mass(Element::monoisotopic_mass)
    }

    fn average_mass(&self) -> Decimal {
        self.mass(Element::average_mass)
    }
}

impl Charged for ChemicalComposition<'_> {
    fn charge(&self) -> Charge {
        self.particle_offset
            .as_ref()
            .map(|(k, c, p)| {
                let sign = Charge::from(*k);
                let c = Charge::from(*c);
                sign * c * p.charge()
            })
            .unwrap_or_default()
    }
}

impl Mz for ChemicalComposition<'_> {}

// Private Helper Methods ==============================================================================================

impl<'a> ChemicalComposition<'a> {
    fn mass(&self, accessor: impl Fn(&Element<'a>) -> Decimal) -> Decimal {
        let element_masses = self
            .chemical_formula
            .iter()
            .map(|(element, count)| Decimal::from(*count) * accessor(element));

        let particle_masses = self
            .particle_offset
            .iter()
            .map(|(offset_kind, count, particle)| {
                Decimal::from(*offset_kind) * Decimal::from(*count) * particle.mass
            });

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
    fn composition_monoisotopic_mass() -> Result<()> {
        // The masses here have been checked against https://mstools.epfl.ch/info/
        let water = ChemicalComposition::new(&DB, "H2O")?;
        assert_eq!(water.monoisotopic_mass(), dec!(18.01056468403));
        let trp_residue = ChemicalComposition::new(&DB, "C11H10ON2")?;
        assert_eq!(trp_residue.monoisotopic_mass(), dec!(186.07931295073));
        let trp_isotopes = ChemicalComposition::new(&DB, "[13C]11H10O[15N]2")?;
        assert_eq!(trp_isotopes.monoisotopic_mass(), dec!(199.1102859254));
        let gm_aeja = ChemicalComposition::new(&DB, "C37H63N7O21+p")?;
        assert_eq!(gm_aeja.monoisotopic_mass(), dec!(942.414978539091));

        // Testing with proton offsets for adducts (checked against https://www.unimod.org/modifications_list.php)
        let p2 = ChemicalComposition::new(&DB, "2p")?;
        let ca2 = ChemicalComposition::new(&DB, "Ca-2e")?;
        assert_eq!(
            ca2.monoisotopic_mass() - p2.monoisotopic_mass(),
            dec!(37.946940769939870)
        );
        let p1 = ChemicalComposition::new(&DB, "p")?;
        let k1 = ChemicalComposition::new(&DB, "K-e")?;
        assert_eq!(
            k1.monoisotopic_mass() - p1.monoisotopic_mass(),
            dec!(37.955881439869935)
        );
        Ok(())
    }

    #[test]
    fn composition_average_mass() -> Result<()> {
        // The masses here have been checked against https://mstools.epfl.ch/info/
        let water = ChemicalComposition::new(&DB, "H2O")?;
        assert_eq!(water.average_mass(), dec!(18.01528643242983260));
        let trp_residue = ChemicalComposition::new(&DB, "C11H10ON2")?;
        assert_eq!(trp_residue.average_mass(), dec!(186.21031375185538640));
        let trp_isotopes = ChemicalComposition::new(&DB, "[13C]11H10O[15N]2")?;
        assert_eq!(trp_isotopes.average_mass(), dec!(199.11593344840605140));
        let gm_aeja = ChemicalComposition::new(&DB, "C37H63N7O21+p")?;
        assert_eq!(gm_aeja.average_mass(), dec!(942.93919804214360795));

        // Testing with proton offsets for adducts (checked against https://www.unimod.org/modifications_list.php)
        let p2 = ChemicalComposition::new(&DB, "2p")?;
        let ca2 = ChemicalComposition::new(&DB, "Ca-2e")?;
        assert_eq!(
            ca2.average_mass() - p2.average_mass(),
            dec!(38.062372417957600)
        );
        let p1 = ChemicalComposition::new(&DB, "p")?;
        let k1 = ChemicalComposition::new(&DB, "K-e")?;
        assert_eq!(
            k1.average_mass() - p1.average_mass(),
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
