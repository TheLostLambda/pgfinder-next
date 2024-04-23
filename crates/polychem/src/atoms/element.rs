use std::{
    convert::identity,
    fmt::{self, Display, Formatter},
};

use crate::{AverageMass, Element, Isotope, Mass, MassNumber, Massive, MonoisotopicMass, Result};

use super::{
    atomic_database::{AtomicDatabase, ElementDescription},
    errors::AtomicLookupError,
};

impl<'a> Element<'a> {
    pub(crate) fn new(
        db: &'a AtomicDatabase,
        symbol: impl AsRef<str>,
    ) -> Result<Self, AtomicLookupError> {
        Self::lookup(db, symbol, None)
    }

    pub(crate) fn new_isotope(
        db: &'a AtomicDatabase,
        symbol: impl AsRef<str>,
        mass_number: MassNumber,
    ) -> Result<Self, AtomicLookupError> {
        Self::lookup(db, symbol, Some(mass_number))
    }

    fn lookup(
        db: &'a AtomicDatabase,
        symbol: impl AsRef<str>,
        mass_number: Option<MassNumber>,
    ) -> Result<Self, AtomicLookupError> {
        let symbol = symbol.as_ref();
        let (symbol, ElementDescription { name, isotopes }) = db
            .elements
            .get_key_value(symbol)
            .ok_or_else(|| AtomicLookupError::element(symbol))?;

        let element = Self {
            symbol,
            name,
            mass_number,
            isotopes,
        };

        Self::validate_isotopes(element)
    }

    fn validate_isotopes(
        element @ Self {
            symbol,
            name,
            mass_number,
            isotopes,
        }: Self,
    ) -> Result<Self, AtomicLookupError> {
        if let Some(mass_number) = mass_number {
            if !isotopes.contains_key(&mass_number) {
                return Err(AtomicLookupError::isotope(
                    symbol,
                    mass_number,
                    name,
                    isotopes,
                ));
            }
        } else if !isotopes.values().any(|i| i.abundance.is_some()) {
            return Err(AtomicLookupError::abundance(name, symbol, isotopes));
        }

        Ok(element)
    }

    fn isotope_mass(&self) -> Option<Mass> {
        self.mass_number
            .and_then(|a| self.isotopes.get(&a))
            .map(|i| i.relative_mass)
    }

    fn isotope_abundances(&self) -> impl Iterator<Item = &Isotope> {
        self.isotopes.values().filter(|i| i.abundance.is_some())
    }
}

impl Display for Element<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let symbol = self.symbol;
        if let Some(mass_number) = self.mass_number {
            write!(f, "[{mass_number}{symbol}]")
        } else {
            write!(f, "{symbol}")
        }
    }
}

impl Massive for Element<'_> {
    fn monoisotopic_mass(&self) -> MonoisotopicMass {
        // SAFETY: The call to `.unwrap()` is safe here since `.isotope_abundances()` is guaranteed to yield at
        // least one isotope
        self.isotope_mass()
            .map_or_else(
                || {
                    self.isotope_abundances()
                        .max_by_key(|i| i.abundance)
                        .unwrap()
                        .relative_mass
                },
                identity,
            )
            .into()
    }

    fn average_mass(&self) -> AverageMass {
        // SAFETY: The call to `.unwrap()` is safe here since `.isotope_abundances()` is guaranteed to yield
        // only isotopes containing natural abundance data
        self.isotope_mass()
            .map_or_else(
                || {
                    self.isotope_abundances()
                        .map(|i| i.relative_mass * i.abundance.unwrap())
                        .sum()
                },
                identity,
            )
            .into()
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::testing_tools::assert_miette_snapshot;

    use super::*;

    static DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);

    fn mn(c: u32) -> MassNumber {
        MassNumber::new(c).unwrap()
    }

    #[test]
    fn new_element() {
        // Sucessfully lookup elements that exist
        let Element {
            symbol,
            name,
            mass_number,
            isotopes,
        } = Element::new(&DB, "C").unwrap();
        assert_eq!(symbol, "C");
        assert_eq!(name, "Carbon");
        assert_eq!(mass_number, None);
        let mut isotopes: Vec<_> = isotopes.iter().collect();
        isotopes.sort_unstable_by(|(a, _), (b, _)| a.cmp(b));
        assert_debug_snapshot!(isotopes);
        // Fail to lookup elements that don't exist
        assert_miette_snapshot!(Element::new(&DB, "R"));
    }

    #[test]
    fn new_isotope() {
        // Sucessfully lookup isotopes that exist
        let Element {
            symbol,
            name,
            mass_number,
            isotopes,
        } = Element::new_isotope(&DB, "C", mn(13)).unwrap();
        assert_eq!(symbol, "C");
        assert_eq!(name, "Carbon");
        assert_eq!(mass_number, MassNumber::new(13));
        let mut isotopes: Vec<_> = isotopes.iter().collect();
        isotopes.sort_unstable_by(|(a, _), (b, _)| a.cmp(b));
        assert_debug_snapshot!(isotopes);
        // Fail to lookup isotopes for elements that don't exist
        assert_miette_snapshot!(Element::new_isotope(&DB, "R", mn(42)));
        // Fail to lookup isotopes that don't exist
        assert_miette_snapshot!(Element::new_isotope(&DB, "C", mn(15)));
    }

    #[test]
    fn element_display() {
        let c = Element::new(&DB, "C").unwrap();
        assert_eq!(c.to_string(), "C");
        let c13 = Element::new_isotope(&DB, "C", mn(13)).unwrap();
        assert_eq!(c13.to_string(), "[13C]");
        let th = Element::new(&DB, "Th").unwrap();
        assert_eq!(th.to_string(), "Th");
        let th230 = Element::new_isotope(&DB, "Th", mn(230)).unwrap();
        assert_eq!(th230.to_string(), "[230Th]");
    }

    #[test]
    fn element_monoisotopic_mass() {
        // Successfully calculate the monoisotopic mass of elements with natural abundances
        let c = Element::new(&DB, "C").unwrap().monoisotopic_mass();
        assert_eq!(c, MonoisotopicMass(dec!(12)));
        let mg = Element::new(&DB, "Mg").unwrap().monoisotopic_mass();
        assert_eq!(mg, MonoisotopicMass(dec!(23.985041697)));
        let mo = Element::new(&DB, "Mo").unwrap().monoisotopic_mass();
        assert_eq!(mo, MonoisotopicMass(dec!(97.90540482)));
        // Fail to construct elements without natural abundances
        assert_miette_snapshot!(Element::new(&DB, "Tc"));
    }

    #[test]
    fn element_average_mass() {
        // Successfully calculate the average mass of elements with natural abundances
        let c = Element::new(&DB, "C").unwrap().average_mass();
        assert_eq!(c, AverageMass(dec!(12.010735896735249)));
        let mg = Element::new(&DB, "Mg").unwrap().average_mass();
        assert_eq!(mg, AverageMass(dec!(24.3050516198371)));
        let mo = Element::new(&DB, "Mo").unwrap().average_mass();
        assert_eq!(mo, AverageMass(dec!(95.959788541188)));
        // Fail to construct elements without natural abundances
        assert_miette_snapshot!(Element::new(&DB, "Po"));
    }

    #[test]
    fn isotope_masses() {
        // Get masses for an element with natural abundances
        let c13_mono = Element::new_isotope(&DB, "C", mn(13))
            .unwrap()
            .monoisotopic_mass();
        assert_eq!(c13_mono, MonoisotopicMass(dec!(13.00335483507)));
        let c13_avg = Element::new_isotope(&DB, "C", mn(13))
            .unwrap()
            .average_mass();
        assert_eq!(c13_avg, AverageMass(dec!(13.00335483507)));
        // Get masses for an element without natural abundances
        let tc99_mono = Element::new_isotope(&DB, "Tc", mn(99))
            .unwrap()
            .monoisotopic_mass();
        assert_eq!(tc99_mono, MonoisotopicMass(dec!(98.9062508)));
        let tc99_avg = Element::new_isotope(&DB, "Tc", mn(99))
            .unwrap()
            .average_mass();
        assert_eq!(tc99_avg, AverageMass(dec!(98.9062508)));
    }
}
