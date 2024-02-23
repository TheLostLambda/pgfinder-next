use rust_decimal::Decimal;

use crate::{Charge, Charged, Massive, Mz, NamedMod, PolychemError, Result};

use super::polymer_database::{ModificationDescription, PolymerDatabase};

impl<'a, 'p> NamedMod<'a, 'p> {
    // TODO: Write new_with_targets that returns a (Self, &Vec<Target>) for Polymerizer to use
    pub fn new(db: &'p PolymerDatabase<'a>, abbr: impl AsRef<str>) -> Result<Self> {
        let abbr = abbr.as_ref();
        let (
            abbr,
            ModificationDescription {
                name, lost, gained, ..
            },
        ) = db
            .modifications
            .get_key_value(abbr)
            .ok_or_else(|| PolychemError::ModificationLookup(abbr.to_owned()))?;
        Ok(Self {
            abbr,
            name,
            lost,
            gained,
        })
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

impl Mz for NamedMod<'_, '_> {}

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::{
        atoms::atomic_database::AtomicDatabase, polymers::polymer_database::PolymerDatabase,
        testing_tools::assert_miette_snapshot, Charged, Massive, Mz, NamedMod,
    };

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(|| {
        AtomicDatabase::from_kdl(
            "atomic_database.kdl",
            include_str!("../../atomic_database.kdl"),
        )
        .unwrap()
    });

    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::from_kdl(
            &ATOMIC_DB,
            "muropeptide_chemistry.kdl",
            include_str!("../../muropeptide_chemistry.kdl"),
        )
        .unwrap()
    });

    #[test]
    fn errors() {
        let calcium = NamedMod::new(&POLYMER_DB, "Ca");
        assert_miette_snapshot!(calcium);
        let potassium = NamedMod::new(&POLYMER_DB, "K");
        assert_miette_snapshot!(potassium);
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
    }

    #[test]
    fn charge() {
        // TODO: It's probably worth creating a database file for these tests that contains some charged modifications!
        let amidation = NamedMod::new(&POLYMER_DB, "Am").unwrap();
        assert_eq!(amidation.charge(), 0);
        let acetylation = NamedMod::new(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(acetylation.charge(), 0);
        let deacetylation = NamedMod::new(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(deacetylation.charge(), 0);
    }

    #[test]
    fn monoisotopic_mz() {
        // TODO: It's probably worth creating a database file for these tests that contains some charged modifications!
        let amidation = NamedMod::new(&POLYMER_DB, "Am").unwrap();
        assert_eq!(amidation.monoisotopic_mz(), None);
        let acetylation = NamedMod::new(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(acetylation.monoisotopic_mz(), None);
        let deacetylation = NamedMod::new(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(deacetylation.monoisotopic_mz(), None);
    }

    #[test]
    fn average_mz() {
        // TODO: It's probably worth creating a database file for these tests that contains some charged modifications!
        let amidation = NamedMod::new(&POLYMER_DB, "Am").unwrap();
        assert_eq!(amidation.average_mz(), None);
        let acetylation = NamedMod::new(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(acetylation.average_mz(), None);
        let deacetylation = NamedMod::new(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(deacetylation.average_mz(), None);
    }
}
