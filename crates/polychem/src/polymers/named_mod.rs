use rust_decimal::Decimal;

use crate::{Charge, Charged, Massive, Modification, NamedMod, PolychemError, Result};

use super::polymer_database::{ModificationDescription, PolymerDatabase};

impl<'a, 'p> NamedMod<'a, 'p> {
    pub fn new(db: &'p PolymerDatabase<'a>, abbr: impl AsRef<str>) -> Result<Self> {
        let (
            abbr,
            ModificationDescription {
                name, lost, gained, ..
            },
        ) = Self::lookup_description(db, abbr)?;
        Ok(Self {
            abbr,
            name,
            lost,
            gained,
        })
    }

    #[must_use]
    pub const fn abbr(&self) -> &'p str {
        self.abbr
    }

    pub(crate) fn lookup_description(
        db: &'p PolymerDatabase<'a>,
        abbr: impl AsRef<str>,
    ) -> Result<(&'p String, &'p ModificationDescription<'a>)> {
        let abbr = abbr.as_ref();
        db.modifications
            .get_key_value(abbr)
            .ok_or_else(|| PolychemError::modification_lookup(abbr).into())
    }
}

impl<'a, 'p> From<NamedMod<'a, 'p>> for Modification<NamedMod<'a, 'p>> {
    fn from(value: NamedMod<'a, 'p>) -> Self {
        Modification::new(1, value)
    }
}

impl Massive for NamedMod<'_, '_> {
    fn monoisotopic_mass(&self) -> Decimal {
        self.gained.monoisotopic_mass() - self.lost.monoisotopic_mass()
    }

    fn average_mass(&self) -> Decimal {
        self.gained.average_mass() - self.lost.average_mass()
    }
}

impl Charged for NamedMod<'_, '_> {
    fn charge(&self) -> Charge {
        self.gained.charge() - self.lost.charge()
    }
}

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::{testing_tools::assert_miette_snapshot, AtomicDatabase, Mz};

    use super::*;

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);

    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "polymer_database.kdl",
            include_str!("../../tests/data/polymer_database.kdl"),
        )
        .unwrap()
    });

    #[test]
    fn errors() {
        let magnesium = NamedMod::new(&POLYMER_DB, "Mg");
        assert_miette_snapshot!(magnesium);
        let potassium = NamedMod::new(&POLYMER_DB, "K");
        assert_miette_snapshot!(potassium);
    }

    #[test]
    fn from_impls() {
        let named_mod = NamedMod::new(&POLYMER_DB, "Am").unwrap();
        let named_modification: Modification<NamedMod> = named_mod.into();
        assert_eq!(
            named_mod.monoisotopic_mass(),
            named_modification.monoisotopic_mass()
        );
    }

    #[test]
    fn monoisotopic_mass() {
        // Masses checked against https://www.unimod.org/modifications_list.php
        let amidation = NamedMod::new(&POLYMER_DB, "Am").unwrap();
        assert_eq!(amidation.monoisotopic_mass(), dec!(-0.98401558291));
        let acetylation = NamedMod::new(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(acetylation.monoisotopic_mass(), dec!(42.01056468403));
        let deacetylation = NamedMod::new(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(deacetylation.monoisotopic_mass(), dec!(-42.01056468403));
        let calcium = NamedMod::new(&POLYMER_DB, "Ca").unwrap();
        assert_eq!(calcium.monoisotopic_mass(), dec!(38.954217236560870));
    }

    #[test]
    fn average_mass() {
        // Masses checked against https://www.unimod.org/modifications_list.php
        let amidation = NamedMod::new(&POLYMER_DB, "Am").unwrap();
        assert_eq!(amidation.average_mass(), dec!(-0.98476095881670255));
        let acetylation = NamedMod::new(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(acetylation.average_mass(), dec!(42.03675822590033060));
        let deacetylation = NamedMod::new(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(deacetylation.average_mass(), dec!(-42.03675822590033060));
        let calcium = NamedMod::new(&POLYMER_DB, "Ca").unwrap();
        assert_eq!(calcium.average_mass(), dec!(39.069648884578600));
    }

    #[test]
    fn charge() {
        let amidation = NamedMod::new(&POLYMER_DB, "Am").unwrap();
        assert_eq!(amidation.charge(), 0);
        let acetylation = NamedMod::new(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(acetylation.charge(), 0);
        let deacetylation = NamedMod::new(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(deacetylation.charge(), 0);
        let calcium = NamedMod::new(&POLYMER_DB, "Ca").unwrap();
        assert_eq!(calcium.charge(), 1);
    }

    #[test]
    fn monoisotopic_mz() {
        let amidation = NamedMod::new(&POLYMER_DB, "Am").unwrap();
        assert_eq!(amidation.monoisotopic_mz(), None);
        let acetylation = NamedMod::new(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(acetylation.monoisotopic_mz(), None);
        let deacetylation = NamedMod::new(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(deacetylation.monoisotopic_mz(), None);
        let calcium = NamedMod::new(&POLYMER_DB, "Ca").unwrap();
        assert_eq!(calcium.monoisotopic_mz(), Some(dec!(38.954217236560870)));
    }

    #[test]
    fn average_mz() {
        let amidation = NamedMod::new(&POLYMER_DB, "Am").unwrap();
        assert_eq!(amidation.average_mz(), None);
        let acetylation = NamedMod::new(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(acetylation.average_mz(), None);
        let deacetylation = NamedMod::new(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(deacetylation.average_mz(), None);
        let calcium = NamedMod::new(&POLYMER_DB, "Ca").unwrap();
        assert_eq!(calcium.average_mz(), Some(dec!(39.069648884578600)));
    }
}
