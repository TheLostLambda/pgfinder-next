use std::fmt::{self, Display, Formatter};

use crate::{
    AnyMod, AverageMass, Charge, Charged, Count, Massive, Modification, MonoisotopicMass, NamedMod,
    OffsetMod,
};

impl<K> Modification<K> {
    // FIXME: Maybe make this `new_with_multiplier`, and make `new(k)` = `new_with_multiplier(1, k)` — depends on if
    // this code ends up as public API (usable outside of the crate)!
    pub(crate) const fn new(multiplier: Count, kind: K) -> Self {
        Self { multiplier, kind }
    }

    // NOTE: This can't exist as a From impl since it overlaps with `impl From<T> for T` in the standard library. This
    // is awfully annoying, so keep an eye out for specialization, which should make that sort of overlap possible:
    // https://github.com/rust-lang/rust/issues/31844
    pub fn convert<K2: From<K>>(self) -> Modification<K2> {
        let Self { multiplier, kind } = self;
        Modification::new(multiplier, kind.into())
    }
}

impl<K: Massive> Massive for Modification<K> {
    fn monoisotopic_mass(&self) -> MonoisotopicMass {
        self.multiplier * self.kind.monoisotopic_mass()
    }

    fn average_mass(&self) -> AverageMass {
        self.multiplier * self.kind.average_mass()
    }
}

impl<K: Charged> Charged for Modification<K> {
    fn charge(&self) -> Charge {
        self.multiplier * self.kind.charge()
    }
}

fn display_offset_modification(
    f: &mut Formatter<'_>,
    multiplier: Count,
    offset: &OffsetMod,
) -> fmt::Result {
    let OffsetMod { kind, composition } = offset;
    if multiplier == Count::default() {
        write!(f, "{kind}{composition}")
    } else {
        write!(f, "{kind}{multiplier}x{composition}")
    }
}

fn display_named_modification(
    f: &mut Formatter<'_>,
    multiplier: Count,
    named: &NamedMod,
) -> fmt::Result {
    let NamedMod { abbr, .. } = named;
    if multiplier == Count::default() {
        write!(f, "{abbr}")
    } else {
        write!(f, "{multiplier}x{abbr}")
    }
}

impl Display for Modification<OffsetMod<'_>> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let &Self {
            multiplier,
            ref kind,
        } = self;
        display_offset_modification(f, multiplier, kind)
    }
}

impl Display for Modification<NamedMod<'_, '_>> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let &Self {
            multiplier,
            ref kind,
        } = self;
        display_named_modification(f, multiplier, kind)
    }
}

impl Display for Modification<AnyMod<'_, '_>> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let &Self {
            multiplier,
            ref kind,
        } = self;
        match kind {
            AnyMod::Named(kind) => display_named_modification(f, multiplier, kind),
            AnyMod::Offset(kind) => display_offset_modification(f, multiplier, kind),
        }
    }
}

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::{
        AtomicDatabase, AverageMz, ChargedParticle, ChemicalComposition, MonoisotopicMz,
        OffsetKind, OffsetMod, PolymerDatabase,
    };

    use super::*;

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);

    static H2O: Lazy<ChemicalComposition> =
        Lazy::new(|| ChemicalComposition::new(&ATOMIC_DB, "H2O").unwrap());
    static CA: Lazy<ChemicalComposition> =
        Lazy::new(|| ChemicalComposition::new(&ATOMIC_DB, "Ca-2e").unwrap());

    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "test_polymer_database.kdl",
            include_str!("../../tests/data/polymer_database.kdl"),
        )
        .unwrap()
    });

    fn c(c: u32) -> Count {
        Count::new(c).unwrap()
    }

    #[test]
    fn modifcation_display() {
        let water_gain = Modification::new(c(1), OffsetMod::new(OffsetKind::Add, H2O.clone()));
        assert_eq!(water_gain.to_string(), "+H2O");
        let water_loss = Modification::new(c(1), OffsetMod::new(OffsetKind::Remove, H2O.clone()));
        assert_eq!(water_loss.to_string(), "-H2O");
        let double_water_gain =
            Modification::new(c(2), OffsetMod::new(OffsetKind::Add, H2O.clone()));
        assert_eq!(double_water_gain.to_string(), "+2xH2O");
        let double_water_loss =
            Modification::new(c(2), OffsetMod::new(OffsetKind::Remove, H2O.clone()));
        assert_eq!(double_water_loss.to_string(), "-2xH2O");

        let amidation = Modification::new(c(1), NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.to_string(), "Am");
        let double_amidation = Modification::new(c(2), NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(double_amidation.to_string(), "2xAm");

        let water_gain = Modification::new(c(1), AnyMod::offset(OffsetKind::Add, H2O.clone()));
        assert_eq!(water_gain.to_string(), "+H2O");
        let water_loss = Modification::new(c(1), AnyMod::offset(OffsetKind::Remove, H2O.clone()));
        assert_eq!(water_loss.to_string(), "-H2O");
        let double_water_gain =
            Modification::new(c(2), AnyMod::offset(OffsetKind::Add, H2O.clone()));
        assert_eq!(double_water_gain.to_string(), "+2xH2O");
        let double_water_loss =
            Modification::new(c(2), AnyMod::offset(OffsetKind::Remove, H2O.clone()));
        assert_eq!(double_water_loss.to_string(), "-2xH2O");

        let amidation = Modification::new(c(1), AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.to_string(), "Am");
        let double_amidation = Modification::new(c(2), AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(double_amidation.to_string(), "2xAm");
    }

    #[test]
    // FIXME: Maybe split this up some?
    // NOTE: This is just a test, so whilst this isn't great, it's acceptable to be a bit long...
    #[allow(clippy::too_many_lines)]
    fn monoisotopic_mass() {
        // Masses checked against https://www.unimod.org/modifications_list.php
        let amidation = Modification::new(c(1), NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(
            amidation.monoisotopic_mass(),
            MonoisotopicMass(dec!(-0.98401558291))
        );
        let acetylation = Modification::new(c(2), NamedMod::new(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(
            acetylation.monoisotopic_mass(),
            MonoisotopicMass(dec!(84.02112936806))
        );
        let deacetylation = Modification::new(c(3), NamedMod::new(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(
            deacetylation.monoisotopic_mass(),
            MonoisotopicMass(dec!(-126.03169405209))
        );
        let calcium = Modification::new(c(4), NamedMod::new(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(
            calcium.monoisotopic_mass(),
            MonoisotopicMass(dec!(155.81686894624348))
        );

        let water_gained = Modification::new(c(4), OffsetMod::new(OffsetKind::Add, H2O.clone()));
        assert_eq!(
            water_gained.monoisotopic_mass(),
            MonoisotopicMass(dec!(72.04225873612))
        );
        let water_lost = Modification::new(c(3), OffsetMod::new(OffsetKind::Remove, H2O.clone()));
        assert_eq!(
            water_lost.monoisotopic_mass(),
            MonoisotopicMass(dec!(-54.03169405209))
        );
        // Masses checked against https://bioportal.bioontology.org/ontologies/UBERON
        let ca_gained = Modification::new(c(2), OffsetMod::new(OffsetKind::Add, CA.clone()));
        assert_eq!(
            ca_gained.monoisotopic_mass(),
            MonoisotopicMass(dec!(79.92298740636374))
        );
        let ca_lost = Modification::new(c(1), OffsetMod::new(OffsetKind::Remove, CA.clone()));
        assert_eq!(
            ca_lost.monoisotopic_mass(),
            MonoisotopicMass(dec!(-39.961493703181870))
        );

        // Masses checked against https://www.unimod.org/modifications_list.php
        let amidation = Modification::new(c(4), AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(
            amidation.monoisotopic_mass(),
            MonoisotopicMass(dec!(-3.93606233164))
        );
        let acetylation = Modification::new(c(3), AnyMod::named(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(
            acetylation.monoisotopic_mass(),
            MonoisotopicMass(dec!(126.03169405209))
        );
        let deacetylation = Modification::new(c(2), AnyMod::named(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(
            deacetylation.monoisotopic_mass(),
            MonoisotopicMass(dec!(-84.02112936806))
        );
        let calcium = Modification::new(c(1), AnyMod::named(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(
            calcium.monoisotopic_mass(),
            MonoisotopicMass(dec!(38.954217236560870))
        );

        let water_gained = Modification::new(c(1), AnyMod::offset(OffsetKind::Add, H2O.clone()));
        assert_eq!(
            water_gained.monoisotopic_mass(),
            MonoisotopicMass(dec!(18.01056468403))
        );
        let water_lost = Modification::new(c(2), AnyMod::offset(OffsetKind::Remove, H2O.clone()));
        assert_eq!(
            water_lost.monoisotopic_mass(),
            MonoisotopicMass(dec!(-36.02112936806))
        );
        // Masses checked against https://bioportal.bioontology.org/ontologies/UBERON
        let ca_gained = Modification::new(c(3), AnyMod::offset(OffsetKind::Add, CA.clone()));
        assert_eq!(
            ca_gained.monoisotopic_mass(),
            MonoisotopicMass(dec!(119.88448110954561))
        );
        let ca_lost = Modification::new(c(4), AnyMod::offset(OffsetKind::Remove, CA.clone()));
        assert_eq!(
            ca_lost.monoisotopic_mass(),
            MonoisotopicMass(dec!(-159.84597481272748))
        );
    }

    #[test]
    // FIXME: Maybe split this up some?
    // NOTE: This is just a test, so whilst this isn't great, it's acceptable to be a bit long...
    #[allow(clippy::too_many_lines)]
    fn average_mass() {
        // Masses checked against https://www.unimod.org/modifications_list.php
        let amidation = Modification::new(c(1), NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(
            amidation.average_mass(),
            AverageMass(dec!(-0.98476095881670255))
        );
        let acetylation = Modification::new(c(2), NamedMod::new(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(
            acetylation.average_mass(),
            AverageMass(dec!(84.07351645180066120))
        );
        let deacetylation = Modification::new(c(3), NamedMod::new(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(
            deacetylation.average_mass(),
            AverageMass(dec!(-126.11027467770099180))
        );
        let calcium = Modification::new(c(4), NamedMod::new(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(
            calcium.average_mass(),
            AverageMass(dec!(156.278595538314400))
        );

        // Masses checked against https://www.unimod.org/modifications_list.php
        let water_gained = Modification::new(c(4), OffsetMod::new(OffsetKind::Add, H2O.clone()));
        assert_eq!(
            water_gained.average_mass(),
            AverageMass(dec!(72.06114572971933040))
        );
        let water_lost = Modification::new(c(3), OffsetMod::new(OffsetKind::Remove, H2O.clone()));
        assert_eq!(
            water_lost.average_mass(),
            AverageMass(dec!(-54.04585929728949780))
        );
        // Masses checked against https://bioportal.bioontology.org/ontologies/UBERON
        let ca_gained = Modification::new(c(2), OffsetMod::new(OffsetKind::Add, CA.clone()));
        assert_eq!(
            ca_gained.average_mass(),
            AverageMass(dec!(80.153850702399200))
        );
        let ca_lost = Modification::new(c(1), OffsetMod::new(OffsetKind::Remove, CA.clone()));
        assert_eq!(
            ca_lost.average_mass(),
            AverageMass(dec!(-40.076925351199600))
        );

        // Masses checked against https://www.unimod.org/modifications_list.php
        let amidation = Modification::new(c(4), AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(
            amidation.average_mass(),
            AverageMass(dec!(-3.93904383526681020))
        );
        let acetylation = Modification::new(c(3), AnyMod::named(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(
            acetylation.average_mass(),
            AverageMass(dec!(126.11027467770099180))
        );
        let deacetylation = Modification::new(c(2), AnyMod::named(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(
            deacetylation.average_mass(),
            AverageMass(dec!(-84.07351645180066120))
        );
        let calcium = Modification::new(c(1), AnyMod::named(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(
            calcium.average_mass(),
            AverageMass(dec!(39.069648884578600))
        );

        let water_gained = Modification::new(c(1), AnyMod::offset(OffsetKind::Add, H2O.clone()));
        assert_eq!(
            water_gained.average_mass(),
            AverageMass(dec!(18.01528643242983260))
        );
        let water_lost = Modification::new(c(2), AnyMod::offset(OffsetKind::Remove, H2O.clone()));
        assert_eq!(
            water_lost.average_mass(),
            AverageMass(dec!(-36.03057286485966520))
        );
        // Masses checked against https://bioportal.bioontology.org/ontologies/UBERON
        let ca_gained = Modification::new(c(3), AnyMod::offset(OffsetKind::Add, CA.clone()));
        assert_eq!(
            ca_gained.average_mass(),
            AverageMass(dec!(120.230776053598800))
        );
        let ca_lost = Modification::new(c(4), AnyMod::offset(OffsetKind::Remove, CA.clone()));
        assert_eq!(
            ca_lost.average_mass(),
            AverageMass(dec!(-160.307701404798400))
        );
    }

    #[test]
    fn charge() {
        let amidation = Modification::new(c(1), NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.charge(), Charge(0));
        let acetylation = Modification::new(c(2), NamedMod::new(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.charge(), Charge(0));
        let deacetylation = Modification::new(c(3), NamedMod::new(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.charge(), Charge(0));
        let calcium = Modification::new(c(4), NamedMod::new(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(calcium.charge(), Charge(4));

        let water_gained = Modification::new(c(4), OffsetMod::new(OffsetKind::Add, H2O.clone()));
        assert_eq!(water_gained.charge(), Charge(0));
        let water_lost = Modification::new(c(3), OffsetMod::new(OffsetKind::Remove, H2O.clone()));
        assert_eq!(water_lost.charge(), Charge(0));
        let ca_gained = Modification::new(c(2), OffsetMod::new(OffsetKind::Add, CA.clone()));
        assert_eq!(ca_gained.charge(), Charge(4));
        let ca_lost = Modification::new(c(1), OffsetMod::new(OffsetKind::Remove, CA.clone()));
        assert_eq!(ca_lost.charge(), Charge(-2));

        let amidation = Modification::new(c(4), AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.charge(), Charge(0));
        let acetylation = Modification::new(c(3), AnyMod::named(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.charge(), Charge(0));
        let deacetylation = Modification::new(c(2), AnyMod::named(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.charge(), Charge(0));
        let calcium = Modification::new(c(1), AnyMod::named(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(calcium.charge(), Charge(1));

        let water_gained = Modification::new(c(1), AnyMod::offset(OffsetKind::Add, H2O.clone()));
        assert_eq!(water_gained.charge(), Charge(0));
        let water_lost = Modification::new(c(2), AnyMod::offset(OffsetKind::Remove, H2O.clone()));
        assert_eq!(water_lost.charge(), Charge(0));
        let ca_gained = Modification::new(c(3), AnyMod::offset(OffsetKind::Add, CA.clone()));
        assert_eq!(ca_gained.charge(), Charge(6));
        let ca_lost = Modification::new(c(4), AnyMod::offset(OffsetKind::Remove, CA.clone()));
        assert_eq!(ca_lost.charge(), Charge(-8));
    }

    #[test]
    fn monoisotopic_mz() {
        let amidation = Modification::new(c(1), NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.monoisotopic_mz(), None);
        let acetylation = Modification::new(c(2), NamedMod::new(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.monoisotopic_mz(), None);
        let deacetylation = Modification::new(c(3), NamedMod::new(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.monoisotopic_mz(), None);
        let calcium = Modification::new(c(4), NamedMod::new(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(
            calcium.monoisotopic_mz(),
            Some(MonoisotopicMz(dec!(38.954217236560870)))
        );

        let water_gained = Modification::new(c(4), OffsetMod::new(OffsetKind::Add, H2O.clone()));
        assert_eq!(water_gained.monoisotopic_mz(), None);
        let water_lost = Modification::new(c(3), OffsetMod::new(OffsetKind::Remove, H2O.clone()));
        assert_eq!(water_lost.monoisotopic_mz(), None);
        let ca_gained = Modification::new(c(2), OffsetMod::new(OffsetKind::Add, CA.clone()));
        assert_eq!(
            ca_gained.monoisotopic_mz(),
            Some(MonoisotopicMz(dec!(19.980746851590935)))
        );
        let ca_lost = Modification::new(c(1), OffsetMod::new(OffsetKind::Remove, CA.clone()));
        assert_eq!(
            ca_lost.monoisotopic_mz(),
            Some(MonoisotopicMz(dec!(-19.980746851590935)))
        );

        let amidation = Modification::new(c(4), AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.monoisotopic_mz(), None);
        let acetylation = Modification::new(c(3), AnyMod::named(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.monoisotopic_mz(), None);
        let deacetylation = Modification::new(c(2), AnyMod::named(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.monoisotopic_mz(), None);
        let calcium = Modification::new(c(1), AnyMod::named(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(
            calcium.monoisotopic_mz(),
            Some(MonoisotopicMz(dec!(38.954217236560870)))
        );

        let water_gained = Modification::new(c(1), AnyMod::offset(OffsetKind::Add, H2O.clone()));
        assert_eq!(water_gained.monoisotopic_mz(), None);
        let water_lost = Modification::new(c(2), AnyMod::offset(OffsetKind::Remove, H2O.clone()));
        assert_eq!(water_lost.monoisotopic_mz(), None);
        let ca_gained = Modification::new(c(3), AnyMod::offset(OffsetKind::Add, CA.clone()));
        assert_eq!(
            ca_gained.monoisotopic_mz(),
            Some(MonoisotopicMz(dec!(19.980746851590935)))
        );
        let ca_lost = Modification::new(c(4), AnyMod::offset(OffsetKind::Remove, CA.clone()));
        assert_eq!(
            ca_lost.monoisotopic_mz(),
            Some(MonoisotopicMz(dec!(-19.980746851590935)))
        );
    }

    #[test]
    fn average_mz() {
        let amidation = Modification::new(c(1), NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.average_mz(), None);
        let acetylation = Modification::new(c(2), NamedMod::new(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.average_mz(), None);
        let deacetylation = Modification::new(c(3), NamedMod::new(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.average_mz(), None);
        let calcium = Modification::new(c(4), NamedMod::new(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(
            calcium.average_mz(),
            Some(AverageMz(dec!(39.069648884578600)))
        );

        let water_gained = Modification::new(c(4), OffsetMod::new(OffsetKind::Add, H2O.clone()));
        assert_eq!(water_gained.average_mz(), None);
        let water_lost = Modification::new(c(3), OffsetMod::new(OffsetKind::Remove, H2O.clone()));
        assert_eq!(water_lost.average_mz(), None);
        let ca_gained = Modification::new(c(2), OffsetMod::new(OffsetKind::Add, CA.clone()));
        assert_eq!(
            ca_gained.average_mz(),
            Some(AverageMz(dec!(20.0384626755998)))
        );
        let ca_lost = Modification::new(c(1), OffsetMod::new(OffsetKind::Remove, CA.clone()));
        assert_eq!(
            ca_lost.average_mz(),
            Some(AverageMz(dec!(-20.0384626755998)))
        );

        let amidation = Modification::new(c(4), AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.average_mz(), None);
        let acetylation = Modification::new(c(3), AnyMod::named(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.average_mz(), None);
        let deacetylation = Modification::new(c(2), AnyMod::named(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.average_mz(), None);
        let calcium = Modification::new(c(1), AnyMod::named(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(
            calcium.average_mz(),
            Some(AverageMz(dec!(39.069648884578600)))
        );

        let water_gained = Modification::new(c(1), AnyMod::offset(OffsetKind::Add, H2O.clone()));
        assert_eq!(water_gained.average_mz(), None);
        let water_lost = Modification::new(c(2), AnyMod::offset(OffsetKind::Remove, H2O.clone()));
        assert_eq!(water_lost.average_mz(), None);
        let ca_gained = Modification::new(c(3), AnyMod::offset(OffsetKind::Add, CA.clone()));
        assert_eq!(
            ca_gained.average_mz(),
            Some(AverageMz(dec!(20.0384626755998)))
        );
        let ca_lost = Modification::new(c(4), AnyMod::offset(OffsetKind::Remove, CA.clone()));
        assert_eq!(
            ca_lost.average_mz(),
            Some(AverageMz(dec!(-20.0384626755998)))
        );
    }
}
