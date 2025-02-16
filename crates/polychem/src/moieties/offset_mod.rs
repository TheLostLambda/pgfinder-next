use crate::{
    AverageMass, Charge, Charged, ChemicalComposition, Count, Massive, Modification,
    MonoisotopicMass, OffsetKind, OffsetMod,
};

impl<'a> OffsetMod<'a> {
    #[must_use]
    pub(crate) const fn new(kind: OffsetKind, composition: ChemicalComposition<'a>) -> Self {
        Self { kind, composition }
    }

    #[must_use]
    pub const fn kind(&self) -> OffsetKind {
        self.kind
    }
}

impl<'a> From<OffsetMod<'a>> for Modification<OffsetMod<'a>> {
    fn from(value: OffsetMod<'a>) -> Self {
        Self::new(Count::default(), value)
    }
}

impl Massive for OffsetMod<'_> {
    fn monoisotopic_mass(&self) -> MonoisotopicMass {
        self.kind.offset(self.composition.monoisotopic_mass())
    }

    fn average_mass(&self) -> AverageMass {
        self.kind.offset(self.composition.average_mass())
    }
}

impl Charged for OffsetMod<'_> {
    fn charge(&self) -> Charge {
        self.kind.offset(self.composition.charge())
    }
}

#[cfg(test)]
mod tests {
    use rust_decimal_macros::dec;
    use std::sync::LazyLock;

    use crate::{AtomicDatabase, AverageMz, ChargedParticle, MonoisotopicMz};

    use super::*;

    static DB: LazyLock<AtomicDatabase> = LazyLock::new(AtomicDatabase::default);

    static H2O: LazyLock<ChemicalComposition> =
        LazyLock::new(|| ChemicalComposition::new(&DB, "H2O").unwrap());
    static CA: LazyLock<ChemicalComposition> =
        LazyLock::new(|| ChemicalComposition::new(&DB, "Ca-2e").unwrap());

    #[test]
    fn from_impls() {
        let offset_mod = OffsetMod::new(OffsetKind::Add, H2O.clone());
        let offset_modification: Modification<OffsetMod> = offset_mod.clone().into();
        assert_eq!(
            offset_mod.monoisotopic_mass(),
            offset_modification.monoisotopic_mass()
        );

        let offset_mod = OffsetMod::new(OffsetKind::Remove, H2O.clone());
        let offset_modification: Modification<OffsetMod> = offset_mod.clone().into();
        assert_eq!(
            offset_mod.monoisotopic_mass(),
            offset_modification.monoisotopic_mass()
        );
    }

    #[test]
    fn monoisotopic_mass() {
        // Masses checked against https://www.unimod.org/modifications_list.php
        let water_gained = OffsetMod::new(OffsetKind::Add, H2O.clone());
        assert_eq!(
            water_gained.monoisotopic_mass(),
            MonoisotopicMass(dec!(18.01056468403))
        );
        let water_lost = OffsetMod::new(OffsetKind::Remove, H2O.clone());
        assert_eq!(
            water_lost.monoisotopic_mass(),
            MonoisotopicMass(dec!(-18.01056468403))
        );
        // Masses checked against https://bioportal.bioontology.org/ontologies/UBERON
        let ca_gained = OffsetMod::new(OffsetKind::Add, CA.clone());
        assert_eq!(
            ca_gained.monoisotopic_mass(),
            MonoisotopicMass(dec!(39.961493703181870))
        );
        let ca_lost = OffsetMod::new(OffsetKind::Remove, CA.clone());
        assert_eq!(
            ca_lost.monoisotopic_mass(),
            MonoisotopicMass(dec!(-39.961493703181870))
        );
    }

    #[test]
    fn average_mass() {
        // Masses checked against https://www.unimod.org/modifications_list.php
        let water_gained = OffsetMod::new(OffsetKind::Add, H2O.clone());
        assert_eq!(
            water_gained.average_mass(),
            AverageMass(dec!(18.01528643242983260))
        );
        let water_lost = OffsetMod::new(OffsetKind::Remove, H2O.clone());
        assert_eq!(
            water_lost.average_mass(),
            AverageMass(dec!(-18.01528643242983260))
        );
        // Masses checked against https://bioportal.bioontology.org/ontologies/UBERON
        let ca_gained = OffsetMod::new(OffsetKind::Add, CA.clone());
        assert_eq!(
            ca_gained.average_mass(),
            AverageMass(dec!(40.076925351199600))
        );
        let ca_lost = OffsetMod::new(OffsetKind::Remove, CA.clone());
        assert_eq!(
            ca_lost.average_mass(),
            AverageMass(dec!(-40.076925351199600))
        );
    }

    #[test]
    fn charge() {
        let water_gained = OffsetMod::new(OffsetKind::Add, H2O.clone());
        assert_eq!(water_gained.charge(), Charge(0));
        let water_lost = OffsetMod::new(OffsetKind::Remove, H2O.clone());
        assert_eq!(water_lost.charge(), Charge(0));
        let ca_gained = OffsetMod::new(OffsetKind::Add, CA.clone());
        assert_eq!(ca_gained.charge(), Charge(2));
        let ca_lost = OffsetMod::new(OffsetKind::Remove, CA.clone());
        assert_eq!(ca_lost.charge(), Charge(-2));
    }

    #[test]
    fn monoisotopic_mz() {
        let water_gained = OffsetMod::new(OffsetKind::Add, H2O.clone());
        assert_eq!(water_gained.monoisotopic_mz(), None);
        let water_lost = OffsetMod::new(OffsetKind::Remove, H2O.clone());
        assert_eq!(water_lost.monoisotopic_mz(), None);
        let ca_gained = OffsetMod::new(OffsetKind::Add, CA.clone());
        assert_eq!(
            ca_gained.monoisotopic_mz(),
            Some(MonoisotopicMz(dec!(19.980746851590935)))
        );
        let ca_lost = OffsetMod::new(OffsetKind::Remove, CA.clone());
        assert_eq!(
            ca_lost.monoisotopic_mz(),
            Some(MonoisotopicMz(dec!(-19.980746851590935)))
        );
    }

    #[test]
    fn average_mz() {
        let water_gained = OffsetMod::new(OffsetKind::Add, H2O.clone());
        assert_eq!(water_gained.average_mz(), None);
        let water_lost = OffsetMod::new(OffsetKind::Remove, H2O.clone());
        assert_eq!(water_lost.average_mz(), None);
        let ca_gained = OffsetMod::new(OffsetKind::Add, CA.clone());
        assert_eq!(
            ca_gained.average_mz(),
            Some(AverageMz(dec!(20.0384626755998)))
        );
        let ca_lost = OffsetMod::new(OffsetKind::Remove, CA.clone());
        assert_eq!(
            ca_lost.average_mz(),
            Some(AverageMz(dec!(-20.0384626755998)))
        );
    }
}
