use crate::{
    AverageMass, Bond, Charge, Charged, Massive, MonoisotopicMass, Result, errors::PolychemError,
};

use super::polymer_database::{BondDescription, PolymerDatabase};

impl<'a, 'p> Bond<'a, 'p> {
    // FIXME: Should this be public / exist at all? I'm using `lookup_description()` everywhere anyways?
    pub fn new(db: &'p PolymerDatabase<'a>, abbr: impl AsRef<str>) -> Result<Self> {
        let (
            abbr,
            BondDescription {
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

    // FIXME: Missing testing!
    #[must_use]
    pub const fn abbr(&self) -> &'p str {
        self.abbr
    }

    // FIXME: Missing testing!
    #[must_use]
    pub const fn name(&self) -> &'p str {
        self.name
    }

    pub(crate) fn lookup_description(
        db: &'p PolymerDatabase<'a>,
        abbr: impl AsRef<str>,
    ) -> Result<(&'p String, &'p BondDescription<'a>)> {
        let abbr = abbr.as_ref();
        db.bonds
            .get_key_value(abbr)
            .ok_or_else(|| PolychemError::bond_lookup(abbr).into())
    }
}

impl Massive for Bond<'_, '_> {
    fn monoisotopic_mass(&self) -> MonoisotopicMass {
        self.gained.monoisotopic_mass() - self.lost.monoisotopic_mass()
    }

    fn average_mass(&self) -> AverageMass {
        self.gained.average_mass() - self.lost.average_mass()
    }
}

impl Charged for Bond<'_, '_> {
    fn charge(&self) -> Charge {
        self.gained.charge() - self.lost.charge()
    }
}

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::{
        AtomicDatabase, AverageMz, ChargedParticle, MonoisotopicMz,
        testing_tools::assert_miette_snapshot,
    };

    use super::*;

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);

    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "test_polymer_database.kdl",
            include_str!("../../tests/data/polymer_database.kdl"),
        )
        .unwrap()
    });

    #[test]
    fn errors() {
        let disulfide = Bond::new(&POLYMER_DB, "Disulfide");
        assert_miette_snapshot!(disulfide);
        let ionic = Bond::new(&POLYMER_DB, "Ionic");
        assert_miette_snapshot!(ionic);
    }

    #[test]
    fn monoisotopic_mass() {
        let glycosidic = Bond::new(&POLYMER_DB, "Gly").unwrap();
        assert_eq!(
            glycosidic.monoisotopic_mass(),
            MonoisotopicMass(dec!(-18.01056468403))
        );
        let stem = Bond::new(&POLYMER_DB, "Stem").unwrap();
        assert_eq!(
            stem.monoisotopic_mass(),
            MonoisotopicMass(dec!(-18.01056468403))
        );
        let charged = Bond::new(&POLYMER_DB, "Chr").unwrap();
        assert_eq!(
            charged.monoisotopic_mass(),
            MonoisotopicMass(dec!(-2.014552933242))
        );
        let persulfide = Bond::new(&POLYMER_DB, "Sulf").unwrap();
        assert_eq!(
            persulfide.monoisotopic_mass(),
            MonoisotopicMass(dec!(29.95642110994))
        );
    }

    #[test]
    fn average_mass() {
        let glycosidic = Bond::new(&POLYMER_DB, "Gly").unwrap();
        assert_eq!(
            glycosidic.average_mass(),
            AverageMass(dec!(-18.01528643242983260))
        );
        let stem = Bond::new(&POLYMER_DB, "Stem").unwrap();
        assert_eq!(
            stem.average_mass(),
            AverageMass(dec!(-18.01528643242983260))
        );
        let charged = Bond::new(&POLYMER_DB, "Chr").unwrap();
        assert_eq!(charged.average_mass(), AverageMass(dec!(-2.014552933242)));
        let persulfide = Bond::new(&POLYMER_DB, "Sulf").unwrap();
        assert_eq!(
            persulfide.average_mass(),
            AverageMass(dec!(30.04890589801550530))
        );
    }

    #[test]
    fn charge() {
        let glycosidic = Bond::new(&POLYMER_DB, "Gly").unwrap();
        assert_eq!(glycosidic.charge(), Charge(0));
        let stem = Bond::new(&POLYMER_DB, "Stem").unwrap();
        assert_eq!(stem.charge(), Charge(0));
        let charged = Bond::new(&POLYMER_DB, "Chr").unwrap();
        assert_eq!(charged.charge(), Charge(-2));
        let persulfide = Bond::new(&POLYMER_DB, "Sulf").unwrap();
        assert_eq!(persulfide.charge(), Charge(0));
    }

    #[test]
    fn monoisotopic_mz() {
        let glycosidic = Bond::new(&POLYMER_DB, "Gly").unwrap();
        assert_eq!(glycosidic.monoisotopic_mz(), None);
        let stem = Bond::new(&POLYMER_DB, "Stem").unwrap();
        assert_eq!(stem.monoisotopic_mz(), None);
        let charged = Bond::new(&POLYMER_DB, "Chr").unwrap();
        assert_eq!(
            charged.monoisotopic_mz(),
            Some(MonoisotopicMz(dec!(-1.007276466621)))
        );
        let persulfide = Bond::new(&POLYMER_DB, "Sulf").unwrap();
        assert_eq!(persulfide.monoisotopic_mz(), None);
    }

    #[test]
    fn average_mz() {
        let glycosidic = Bond::new(&POLYMER_DB, "Gly").unwrap();
        assert_eq!(glycosidic.average_mz(), None);
        let stem = Bond::new(&POLYMER_DB, "Stem").unwrap();
        assert_eq!(stem.average_mz(), None);
        let charged = Bond::new(&POLYMER_DB, "Chr").unwrap();
        assert_eq!(charged.average_mz(), Some(AverageMz(dec!(-1.007276466621))));
        let persulfide = Bond::new(&POLYMER_DB, "Sulf").unwrap();
        assert_eq!(persulfide.average_mz(), None);
    }
}
