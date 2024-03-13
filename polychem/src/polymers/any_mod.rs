use rust_decimal::Decimal;

use crate::{
    atoms::atomic_database::AtomicDatabase, AnyMod, Charge, Charged, Massive, Mz, NamedMod,
    OffsetKind, OffsetMod, Result,
};

use super::polymer_database::PolymerDatabase;

impl<'a, 'p> AnyMod<'a, 'p> {
    pub fn named(db: &'p PolymerDatabase<'a>, abbr: impl AsRef<str>) -> Result<Self> {
        Ok(Self::Named(NamedMod::new(db, abbr)?))
    }

    pub fn offset(
        db: &'a AtomicDatabase,
        kind: OffsetKind,
        formula: impl AsRef<str>,
    ) -> Result<Self> {
        Ok(Self::Offset(OffsetMod::new(db, kind, formula)?))
    }
}

// NOTE: There are crates for automating more of this code generation, but odds are that I'll only need to do this sort
// of enum dispatch for AnyMod — it doesn't seem worth a dependency and cluttering lib.rs with attributes
macro_rules! dispatch {
    ($self:expr, $method:ident) => {
        match $self {
            AnyMod::Named(m) => m.$method(),
            AnyMod::Offset(m) => m.$method(),
        }
    };
}

impl Massive for AnyMod<'_, '_> {
    fn monoisotopic_mass(&self) -> Decimal {
        dispatch!(self, monoisotopic_mass)
    }

    fn average_mass(&self) -> Decimal {
        dispatch!(self, average_mass)
    }
}

impl Charged for AnyMod<'_, '_> {
    fn charge(&self) -> Charge {
        dispatch!(self, charge)
    }
}

impl Mz for AnyMod<'_, '_> {}

impl<'a, 'p> From<NamedMod<'a, 'p>> for AnyMod<'a, 'p> {
    fn from(value: NamedMod<'a, 'p>) -> Self {
        Self::Named(value)
    }
}

impl<'a, 'p> From<OffsetMod<'a>> for AnyMod<'a, 'p> {
    fn from(value: OffsetMod<'a>) -> Self {
        Self::Offset(value)
    }
}

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::{
        atoms::atomic_database::AtomicDatabase, polymers::polymer_database::PolymerDatabase,
        testing_tools::assert_miette_snapshot, AnyMod, Charged, Massive, Mz, OffsetKind,
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
        let calcium = AnyMod::named(&POLYMER_DB, "Ca");
        assert_miette_snapshot!(calcium);
        let potassium = AnyMod::named(&POLYMER_DB, "K");
        assert_miette_snapshot!(potassium);
        let water_gained = AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "H[2O]");
        assert_miette_snapshot!(water_gained);
        let water_lost = AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "H[2O]");
        assert_miette_snapshot!(water_lost);
    }

    #[test]
    fn monoisotopic_mass() {
        // Masses checked against https://www.unimod.org/modifications_list.php
        let amidation = AnyMod::named(&POLYMER_DB, "Am").unwrap();
        assert_eq!(amidation.monoisotopic_mass(), dec!(-0.98401558291));
        let acetylation = AnyMod::named(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(acetylation.monoisotopic_mass(), dec!(42.01056468403));
        let deacetylation = AnyMod::named(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(deacetylation.monoisotopic_mass(), dec!(-42.01056468403));

        let water_gained = AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap();
        assert_eq!(water_gained.monoisotopic_mass(), dec!(18.01056468403));
        let water_lost = AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap();
        assert_eq!(water_lost.monoisotopic_mass(), dec!(-18.01056468403));
        // Masses checked against https://bioportal.bioontology.org/ontologies/UBERON
        let ca_gained = AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap();
        assert_eq!(ca_gained.monoisotopic_mass(), dec!(39.961493703181870));
        let ca_lost = AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap();
        assert_eq!(ca_lost.monoisotopic_mass(), dec!(-39.961493703181870));
    }

    #[test]
    fn average_mass() {
        // Masses checked against https://www.unimod.org/modifications_list.php
        let amidation = AnyMod::named(&POLYMER_DB, "Am").unwrap();
        assert_eq!(amidation.average_mass(), dec!(-0.98476095881670255));
        let acetylation = AnyMod::named(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(acetylation.average_mass(), dec!(42.03675822590033060));
        let deacetylation = AnyMod::named(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(deacetylation.average_mass(), dec!(-42.03675822590033060));

        let water_gained = AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap();
        assert_eq!(water_gained.average_mass(), dec!(18.01528643242983260));
        let water_lost = AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap();
        assert_eq!(water_lost.average_mass(), dec!(-18.01528643242983260));
        // Masses checked against https://bioportal.bioontology.org/ontologies/UBERON
        let ca_gained = AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap();
        assert_eq!(ca_gained.average_mass(), dec!(40.076925351199600));
        let ca_lost = AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap();
        assert_eq!(ca_lost.average_mass(), dec!(-40.076925351199600));
    }

    #[test]
    fn charge() {
        // TODO: It's probably worth creating a database file for these tests that contains some charged modifications!
        let amidation = AnyMod::named(&POLYMER_DB, "Am").unwrap();
        assert_eq!(amidation.charge(), 0);
        let acetylation = AnyMod::named(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(acetylation.charge(), 0);
        let deacetylation = AnyMod::named(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(deacetylation.charge(), 0);

        let water_gained = AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap();
        assert_eq!(water_gained.charge(), 0);
        let water_lost = AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap();
        assert_eq!(water_lost.charge(), 0);
        let ca_gained = AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap();
        assert_eq!(ca_gained.charge(), 2);
        let ca_lost = AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap();
        assert_eq!(ca_lost.charge(), -2);
    }

    #[test]
    fn monoisotopic_mz() {
        // TODO: It's probably worth creating a database file for these tests that contains some charged modifications!
        let amidation = AnyMod::named(&POLYMER_DB, "Am").unwrap();
        assert_eq!(amidation.monoisotopic_mz(), None);
        let acetylation = AnyMod::named(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(acetylation.monoisotopic_mz(), None);
        let deacetylation = AnyMod::named(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(deacetylation.monoisotopic_mz(), None);

        let water_gained = AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap();
        assert_eq!(water_gained.monoisotopic_mz(), None);
        let water_lost = AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap();
        assert_eq!(water_lost.monoisotopic_mz(), None);
        let ca_gained = AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap();
        assert_eq!(ca_gained.monoisotopic_mz(), Some(dec!(19.980746851590935)));
        let ca_lost = AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap();
        assert_eq!(ca_lost.monoisotopic_mz(), Some(dec!(-19.980746851590935)));
    }

    #[test]
    fn average_mz() {
        // TODO: It's probably worth creating a database file for these tests that contains some charged modifications!
        let amidation = AnyMod::named(&POLYMER_DB, "Am").unwrap();
        assert_eq!(amidation.average_mz(), None);
        let acetylation = AnyMod::named(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(acetylation.average_mz(), None);
        let deacetylation = AnyMod::named(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(deacetylation.average_mz(), None);

        let water_gained = AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap();
        assert_eq!(water_gained.average_mz(), None);
        let water_lost = AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap();
        assert_eq!(water_lost.average_mz(), None);
        let ca_gained = AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap();
        assert_eq!(ca_gained.average_mz(), Some(dec!(20.0384626755998)));
        let ca_lost = AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap();
        assert_eq!(ca_lost.average_mz(), Some(dec!(-20.0384626755998)));
    }
}
