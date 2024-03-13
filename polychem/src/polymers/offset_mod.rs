use rust_decimal::Decimal;

use crate::{
    atoms::atomic_database::AtomicDatabase, Charge, Charged, ChemicalComposition, Massive, Mz,
    OffsetKind, OffsetMod, Result,
};

impl<'a> OffsetMod<'a> {
    // FIXME: Should this error be wrapped? Probably should have a variant in PolychemError!
    pub fn new(db: &'a AtomicDatabase, kind: OffsetKind, formula: impl AsRef<str>) -> Result<Self> {
        let composition = ChemicalComposition::new(db, formula)?;
        Ok(Self { kind, composition })
    }

    pub fn new_with_composition(kind: OffsetKind, composition: ChemicalComposition<'a>) -> Self {
        Self { kind, composition }
    }
}

impl Massive for OffsetMod<'_> {
    fn monoisotopic_mass(&self) -> Decimal {
        Decimal::from(self.kind) * self.composition.monoisotopic_mass()
    }

    fn average_mass(&self) -> Decimal {
        Decimal::from(self.kind) * self.composition.average_mass()
    }
}

impl Charged for OffsetMod<'_> {
    fn charge(&self) -> Charge {
        Charge::from(self.kind) * self.composition.charge()
    }
}

impl Mz for OffsetMod<'_> {}

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::{
        atoms::atomic_database::AtomicDatabase, testing_tools::assert_miette_snapshot, Charged,
        Massive, Mz, OffsetKind, OffsetMod,
    };

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(|| {
        AtomicDatabase::from_kdl(
            "atomic_database.kdl",
            include_str!("../../atomic_database.kdl"),
        )
        .unwrap()
    });

    #[test]
    fn errors() {
        let water_gained = OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "H[2O]");
        assert_miette_snapshot!(water_gained);
        let water_lost = OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "H[2O]");
        assert_miette_snapshot!(water_lost);
    }

    #[test]
    fn monoisotopic_mass() {
        // Masses checked against https://www.unimod.org/modifications_list.php
        let water_gained = OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap();
        assert_eq!(water_gained.monoisotopic_mass(), dec!(18.01056468403));
        let water_lost = OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap();
        assert_eq!(water_lost.monoisotopic_mass(), dec!(-18.01056468403));
        // Masses checked against https://bioportal.bioontology.org/ontologies/UBERON
        let ca_gained = OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap();
        assert_eq!(ca_gained.monoisotopic_mass(), dec!(39.961493703181870));
        let ca_lost = OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap();
        assert_eq!(ca_lost.monoisotopic_mass(), dec!(-39.961493703181870));
    }

    #[test]
    fn average_mass() {
        // Masses checked against https://www.unimod.org/modifications_list.php
        let water_gained = OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap();
        assert_eq!(water_gained.average_mass(), dec!(18.01528643242983260));
        let water_lost = OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap();
        assert_eq!(water_lost.average_mass(), dec!(-18.01528643242983260));
        // Masses checked against https://bioportal.bioontology.org/ontologies/UBERON
        let ca_gained = OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap();
        assert_eq!(ca_gained.average_mass(), dec!(40.076925351199600));
        let ca_lost = OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap();
        assert_eq!(ca_lost.average_mass(), dec!(-40.076925351199600));
    }

    #[test]
    fn charge() {
        let water_gained = OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap();
        assert_eq!(water_gained.charge(), 0);
        let water_lost = OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap();
        assert_eq!(water_lost.charge(), 0);
        let ca_gained = OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap();
        assert_eq!(ca_gained.charge(), 2);
        let ca_lost = OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap();
        assert_eq!(ca_lost.charge(), -2);
    }

    #[test]
    fn monoisotopic_mz() {
        let water_gained = OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap();
        assert_eq!(water_gained.monoisotopic_mz(), None);
        let water_lost = OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap();
        assert_eq!(water_lost.monoisotopic_mz(), None);
        let ca_gained = OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap();
        assert_eq!(ca_gained.monoisotopic_mz(), Some(dec!(19.980746851590935)));
        let ca_lost = OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap();
        assert_eq!(ca_lost.monoisotopic_mz(), Some(dec!(-19.980746851590935)));
    }

    #[test]
    fn average_mz() {
        let water_gained = OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap();
        assert_eq!(water_gained.average_mz(), None);
        let water_lost = OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap();
        assert_eq!(water_lost.average_mz(), None);
        let ca_gained = OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap();
        assert_eq!(ca_gained.average_mz(), Some(dec!(20.0384626755998)));
        let ca_lost = OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "Ca-2e").unwrap();
        assert_eq!(ca_lost.average_mz(), Some(dec!(-20.0384626755998)));
    }
}
