use rust_decimal::prelude::Decimal;

use crate::{Charge, Charged, Count, Massive, Modification, Mz};

impl<K> Modification<K> {
    // FIXME: Maybe make this `new_with_multiplier`, and make `new(k)` = `new_with_multiplier(1, k)` — depends on if
    // this code ends up as public API (usable outside of the crate)!
    pub const fn new(multiplier: Count, kind: K) -> Self {
        Self { multiplier, kind }
    }
}

impl<K: Massive> Massive for Modification<K> {
    fn monoisotopic_mass(&self) -> Decimal {
        Decimal::from(self.multiplier) * self.kind.monoisotopic_mass()
    }

    fn average_mass(&self) -> Decimal {
        Decimal::from(self.multiplier) * self.kind.average_mass()
    }
}

impl<K: Charged> Charged for Modification<K> {
    fn charge(&self) -> Charge {
        Charge::from(self.multiplier) * self.kind.charge()
    }
}

impl<K: Mz> Mz for Modification<K> {}

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::{
        atoms::atomic_database::AtomicDatabase, polymers::polymer_database::PolymerDatabase,
        AnyMod, Charged, Massive, Modification, Mz, NamedMod, OffsetKind, OffsetMod,
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
    fn monoisotopic_mass() {
        // Masses checked against https://www.unimod.org/modifications_list.php
        let amidation = Modification::new(1, NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.monoisotopic_mass(), dec!(-0.98401558291));
        let acetylation = Modification::new(2, NamedMod::new(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.monoisotopic_mass(), dec!(84.02112936806));
        let deacetylation = Modification::new(3, NamedMod::new(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.monoisotopic_mass(), dec!(-126.03169405209));

        let water_gained = Modification::new(
            4,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap(),
        );
        assert_eq!(water_gained.monoisotopic_mass(), dec!(72.04225873612));
        let water_lost = Modification::new(
            3,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap(),
        );
        assert_eq!(water_lost.monoisotopic_mass(), dec!(-54.03169405209));
        // Masses checked against https://bioportal.bioontology.org/ontologies/UBERON
        let ca_gained = Modification::new(
            2,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_gained.monoisotopic_mass(), dec!(79.92298740636374));
        let ca_lost = Modification::new(
            1,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_lost.monoisotopic_mass(), dec!(-39.961493703181870));

        // Masses checked against https://www.unimod.org/modifications_list.php
        let amidation = Modification::new(3, AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.monoisotopic_mass(), dec!(-2.95204674873));
        let acetylation = Modification::new(2, AnyMod::named(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.monoisotopic_mass(), dec!(84.02112936806));
        let deacetylation = Modification::new(1, AnyMod::named(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.monoisotopic_mass(), dec!(-42.01056468403));

        let water_gained = Modification::new(
            1,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap(),
        );
        assert_eq!(water_gained.monoisotopic_mass(), dec!(18.01056468403));
        let water_lost = Modification::new(
            2,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap(),
        );
        assert_eq!(water_lost.monoisotopic_mass(), dec!(-36.02112936806));
        // Masses checked against https://bioportal.bioontology.org/ontologies/UBERON
        let ca_gained = Modification::new(
            3,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_gained.monoisotopic_mass(), dec!(119.88448110954561));
        let ca_lost = Modification::new(
            4,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_lost.monoisotopic_mass(), dec!(-159.84597481272748));
    }

    #[test]
    fn average_mass() {
        // Masses checked against https://www.unimod.org/modifications_list.php
        let amidation = Modification::new(1, NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.average_mass(), dec!(-0.98476095881670255));
        let acetylation = Modification::new(2, NamedMod::new(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.average_mass(), dec!(84.07351645180066120));
        let deacetylation = Modification::new(3, NamedMod::new(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.average_mass(), dec!(-126.11027467770099180));

        // Masses checked against https://www.unimod.org/modifications_list.php
        let water_gained = Modification::new(
            4,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap(),
        );
        assert_eq!(water_gained.average_mass(), dec!(72.06114572971933040));
        let water_lost = Modification::new(
            3,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap(),
        );
        assert_eq!(water_lost.average_mass(), dec!(-54.04585929728949780));
        // Masses checked against https://bioportal.bioontology.org/ontologies/UBERON
        let ca_gained = Modification::new(
            2,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_gained.average_mass(), dec!(80.153850702399200));
        let ca_lost = Modification::new(
            1,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_lost.average_mass(), dec!(-40.076925351199600));

        // Masses checked against https://www.unimod.org/modifications_list.php
        let amidation = Modification::new(3, AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.average_mass(), dec!(-2.95428287645010765));
        let acetylation = Modification::new(2, AnyMod::named(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.average_mass(), dec!(84.07351645180066120));
        let deacetylation = Modification::new(1, AnyMod::named(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.average_mass(), dec!(-42.03675822590033060));

        let water_gained = Modification::new(
            1,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap(),
        );
        assert_eq!(water_gained.average_mass(), dec!(18.01528643242983260));
        let water_lost = Modification::new(
            2,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap(),
        );
        assert_eq!(water_lost.average_mass(), dec!(-36.03057286485966520));
        // Masses checked against https://bioportal.bioontology.org/ontologies/UBERON
        let ca_gained = Modification::new(
            3,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_gained.average_mass(), dec!(120.230776053598800));
        let ca_lost = Modification::new(
            4,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_lost.average_mass(), dec!(-160.307701404798400));
    }

    #[test]
    fn charge() {
        // TODO: It's probably worth creating a database file for these tests that contains some charged modifications!
        let amidation = Modification::new(1, NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.charge(), 0);
        let acetylation = Modification::new(2, NamedMod::new(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.charge(), 0);
        let deacetylation = Modification::new(3, NamedMod::new(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.charge(), 0);

        let water_gained = Modification::new(
            4,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap(),
        );
        assert_eq!(water_gained.charge(), 0);
        let water_lost = Modification::new(
            3,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap(),
        );
        assert_eq!(water_lost.charge(), 0);
        let ca_gained = Modification::new(
            2,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_gained.charge(), 4);
        let ca_lost = Modification::new(
            1,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_lost.charge(), -2);

        // TODO: It's probably worth creating a database file for these tests that contains some charged modifications!
        let amidation = Modification::new(3, AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.charge(), 0);
        let acetylation = Modification::new(2, AnyMod::named(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.charge(), 0);
        let deacetylation = Modification::new(1, AnyMod::named(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.charge(), 0);

        let water_gained = Modification::new(
            1,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap(),
        );
        assert_eq!(water_gained.charge(), 0);
        let water_lost = Modification::new(
            2,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap(),
        );
        assert_eq!(water_lost.charge(), 0);
        let ca_gained = Modification::new(
            3,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_gained.charge(), 6);
        let ca_lost = Modification::new(
            4,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_lost.charge(), -8);
    }

    #[test]
    fn monoisotopic_mz() {
        // TODO: It's probably worth creating a database file for these tests that contains some charged modifications!
        let amidation = Modification::new(1, NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.monoisotopic_mz(), None);
        let acetylation = Modification::new(2, NamedMod::new(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.monoisotopic_mz(), None);
        let deacetylation = Modification::new(3, NamedMod::new(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.monoisotopic_mz(), None);

        let water_gained = Modification::new(
            4,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap(),
        );
        assert_eq!(water_gained.monoisotopic_mz(), None);
        let water_lost = Modification::new(
            3,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap(),
        );
        assert_eq!(water_lost.monoisotopic_mz(), None);
        let ca_gained = Modification::new(
            2,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_gained.monoisotopic_mz(), Some(dec!(19.980746851590935)));
        let ca_lost = Modification::new(
            1,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_lost.monoisotopic_mz(), Some(dec!(-19.980746851590935)));

        // TODO: It's probably worth creating a database file for these tests that contains some charged modifications!
        let amidation = Modification::new(3, AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.monoisotopic_mz(), None);
        let acetylation = Modification::new(2, AnyMod::named(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.monoisotopic_mz(), None);
        let deacetylation = Modification::new(1, AnyMod::named(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.monoisotopic_mz(), None);

        let water_gained = Modification::new(
            1,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap(),
        );
        assert_eq!(water_gained.monoisotopic_mz(), None);
        let water_lost = Modification::new(
            2,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap(),
        );
        assert_eq!(water_lost.monoisotopic_mz(), None);
        let ca_gained = Modification::new(
            3,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_gained.monoisotopic_mz(), Some(dec!(19.980746851590935)));
        let ca_lost = Modification::new(
            4,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_lost.monoisotopic_mz(), Some(dec!(-19.980746851590935)));
    }

    #[test]
    fn average_mz() {
        // TODO: It's probably worth creating a database file for these tests that contains some charged modifications!
        let amidation = Modification::new(1, NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.average_mz(), None);
        let acetylation = Modification::new(2, NamedMod::new(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.average_mz(), None);
        let deacetylation = Modification::new(3, NamedMod::new(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.average_mz(), None);

        let water_gained = Modification::new(
            4,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap(),
        );
        assert_eq!(water_gained.average_mz(), None);
        let water_lost = Modification::new(
            3,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap(),
        );
        assert_eq!(water_lost.average_mz(), None);
        let ca_gained = Modification::new(
            2,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_gained.average_mz(), Some(dec!(20.0384626755998)));
        let ca_lost = Modification::new(
            1,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_lost.average_mz(), Some(dec!(-20.0384626755998)));

        // TODO: It's probably worth creating a database file for these tests that contains some charged modifications!
        let amidation = Modification::new(3, AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.average_mz(), None);
        let acetylation = Modification::new(2, AnyMod::named(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.average_mz(), None);
        let deacetylation = Modification::new(1, AnyMod::named(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.average_mz(), None);

        let water_gained = Modification::new(
            1,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap(),
        );
        assert_eq!(water_gained.average_mz(), None);
        let water_lost = Modification::new(
            2,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap(),
        );
        assert_eq!(water_lost.average_mz(), None);
        let ca_gained = Modification::new(
            3,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_gained.average_mz(), Some(dec!(20.0384626755998)));
        let ca_lost = Modification::new(
            4,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap(),
        );
        assert_eq!(ca_lost.average_mz(), Some(dec!(-20.0384626755998)));
    }
}
