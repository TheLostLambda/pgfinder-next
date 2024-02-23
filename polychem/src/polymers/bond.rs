use rust_decimal::Decimal;

use crate::{Bond, BondTarget, Charge, Charged, Massive, Mz, PolychemError, Result};

use super::polymer_database::{BondDescription, PolymerDatabase};

impl<'a, 'p> Bond<'a, 'p> {
    // TODO: Write new_with_targets that returns a (Self, &Target, &Target) for Polymerizer to use
    pub fn new(
        db: &'p PolymerDatabase<'a>,
        kind: impl AsRef<str>,
        acceptor: BondTarget<'p>,
    ) -> Result<Self> {
        let kind = kind.as_ref();
        let (kind, BondDescription { lost, .. }) = db
            .bonds
            .get_key_value(kind)
            .ok_or_else(|| PolychemError::BondLookup(kind.to_owned()))?;
        Ok(Self {
            kind,
            lost,
            acceptor,
        })
    }
}

impl Massive for Bond<'_, '_> {
    fn monoisotopic_mass(&self) -> Decimal {
        -self.lost.monoisotopic_mass()
    }

    fn average_mass(&self) -> Decimal {
        -self.lost.average_mass()
    }
}

impl Charged for Bond<'_, '_> {
    fn charge(&self) -> Charge {
        -self.lost.charge()
    }
}

impl Mz for Bond<'_, '_> {}

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::{
        atoms::atomic_database::AtomicDatabase, polymers::polymer_database::PolymerDatabase,
        testing_tools::assert_miette_snapshot, Bond, BondTarget, Charged, Massive, Mz,
    };

    const EMPTY_TARGET: BondTarget = BondTarget {
        residue: 0,
        group_location: "",
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
        let disulfide = Bond::new(&POLYMER_DB, "Disulfide", EMPTY_TARGET);
        assert_miette_snapshot!(disulfide);
        let ionic = Bond::new(&POLYMER_DB, "Ionic", EMPTY_TARGET);
        assert_miette_snapshot!(ionic);
    }

    #[test]
    fn monoisotopic_mass() {
        // TODO: It's worth creating a database file for these tests that contains some bonds with different losses!
        let glycosidic = Bond::new(&POLYMER_DB, "Glycosidic", EMPTY_TARGET).unwrap();
        assert_eq!(glycosidic.monoisotopic_mass(), dec!(-18.01056468403));
        let stem = Bond::new(&POLYMER_DB, "Stem", EMPTY_TARGET).unwrap();
        assert_eq!(stem.monoisotopic_mass(), dec!(-18.01056468403));
    }

    #[test]
    fn average_mass() {
        // TODO: It's worth creating a database file for these tests that contains some bonds with different losses!
        let glycosidic = Bond::new(&POLYMER_DB, "Glycosidic", EMPTY_TARGET).unwrap();
        assert_eq!(glycosidic.average_mass(), dec!(-18.01528643242983260));
        let stem = Bond::new(&POLYMER_DB, "Stem", EMPTY_TARGET).unwrap();
        assert_eq!(stem.average_mass(), dec!(-18.01528643242983260));
    }

    #[test]
    fn charge() {
        // TODO: It's probably worth creating a database file for these tests that contains some charged bond losses!
        let glycosidic = Bond::new(&POLYMER_DB, "Glycosidic", EMPTY_TARGET).unwrap();
        assert_eq!(glycosidic.charge(), 0);
        let stem = Bond::new(&POLYMER_DB, "Stem", EMPTY_TARGET).unwrap();
        assert_eq!(stem.charge(), 0);
    }

    #[test]
    fn monoisotopic_mz() {
        // TODO: It's probably worth creating a database file for these tests that contains some charged bond losses!
        let glycosidic = Bond::new(&POLYMER_DB, "Glycosidic", EMPTY_TARGET).unwrap();
        assert_eq!(glycosidic.monoisotopic_mz(), None);
        let stem = Bond::new(&POLYMER_DB, "Stem", EMPTY_TARGET).unwrap();
        assert_eq!(stem.monoisotopic_mz(), None);
    }

    #[test]
    fn average_mz() {
        // TODO: It's probably worth creating a database file for these tests that contains some charged bond losses!
        let glycosidic = Bond::new(&POLYMER_DB, "Glycosidic", EMPTY_TARGET).unwrap();
        assert_eq!(glycosidic.average_mz(), None);
        let stem = Bond::new(&POLYMER_DB, "Stem", EMPTY_TARGET).unwrap();
        assert_eq!(stem.average_mz(), None);
    }
}
