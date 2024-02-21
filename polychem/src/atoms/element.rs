use itertools::Itertools;
use rust_decimal::Decimal;

use crate::{Element, Isotope, MassNumber, PolychemError, Result};

use super::{
    atomic_database::{AtomicDatabase, ElementDescription},
    AtomicLookupError,
};

impl<'a> Element<'a> {
    pub(super) fn new(
        db: &'a AtomicDatabase,
        symbol: impl AsRef<str>,
    ) -> std::result::Result<Self, AtomicLookupError> {
        let symbol = symbol.as_ref();
        let (symbol, ElementDescription { name, isotopes }) = db
            .elements
            .get_key_value(symbol)
            .ok_or_else(|| AtomicLookupError::Element(symbol.to_owned()))?;
        Ok(Self {
            symbol,
            name,
            mass_number: None,
            isotopes,
        })
    }

    pub(super) fn new_isotope(
        db: &'a AtomicDatabase,
        symbol: impl AsRef<str>,
        mass_number: MassNumber,
    ) -> std::result::Result<Self, AtomicLookupError> {
        let symbol = symbol.as_ref();
        let element = Self::new(db, symbol)?;
        if element.isotopes.contains_key(&mass_number) {
            Ok(Self {
                mass_number: Some(mass_number),
                ..element
            })
        } else {
            Err(AtomicLookupError::Isotope(
                symbol.to_owned(),
                mass_number,
                element.name.to_owned(),
                element.isotopes.keys().copied().sorted().collect(),
            ))
        }
    }

    pub(super) fn monoisotopic_mass(&self) -> Result<Decimal> {
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

    pub(super) fn average_mass(&self) -> Result<Decimal> {
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

    fn isotope_abundances(
        &self,
        // FIXME: Qualify std::result::Result as just std::Result or something?
    ) -> std::result::Result<impl Iterator<Item = &Isotope>, AtomicLookupError> {
        let mut isotopes_with_abundances = self
            .isotopes
            .values()
            .filter(|i| i.abundance.is_some())
            .peekable();
        if isotopes_with_abundances.peek().is_some() {
            Ok(isotopes_with_abundances)
        } else {
            Err(AtomicLookupError::Abundance(
                self.name.to_owned(),
                self.symbol.to_owned(),
                self.isotopes.keys().copied().sorted().collect(),
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use miette::Result;
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::testing_tools::assert_miette_snapshot;

    use super::{AtomicDatabase, Element};

    static DB: Lazy<AtomicDatabase> = Lazy::new(|| {
        AtomicDatabase::from_kdl(
            "atomic_database.kdl",
            include_str!("../../atomic_database.kdl"),
        )
        .unwrap()
    });

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
        let mut isotopes: Vec<_> = isotopes.iter().collect();
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
        let mut isotopes: Vec<_> = isotopes.iter().collect();
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
}
