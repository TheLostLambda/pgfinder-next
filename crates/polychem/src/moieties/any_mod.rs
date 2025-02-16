use crate::{
    AnyMod, AnyModification, AverageMass, Charge, Charged, ChemicalComposition, Count, Massive,
    Modification, MonoisotopicMass, NamedMod, OffsetKind, OffsetMod, PolymerDatabase, Result,
};

impl<'a, 'p> AnyMod<'a, 'p> {
    pub(crate) fn named(db: &'p PolymerDatabase<'a>, abbr: impl AsRef<str>) -> Result<Self> {
        Ok(Self::Named(NamedMod::new(db, abbr)?))
    }

    #[must_use]
    pub const fn offset(kind: OffsetKind, composition: ChemicalComposition<'a>) -> Self {
        Self::Offset(OffsetMod::new(kind, composition))
    }
}

impl<'a, 'p> From<NamedMod<'a, 'p>> for AnyMod<'a, 'p> {
    fn from(value: NamedMod<'a, 'p>) -> Self {
        Self::Named(value)
    }
}

impl<'a> From<OffsetMod<'a>> for AnyMod<'a, '_> {
    fn from(value: OffsetMod<'a>) -> Self {
        Self::Offset(value)
    }
}

impl<'a, 'p, K: Into<AnyMod<'a, 'p>>> From<K> for AnyModification<'a, 'p> {
    fn from(value: K) -> Self {
        Self::new(Count::default(), value.into())
    }
}

// NOTE: I can't merge the following `From` impls since it would overlap with `impl From<T> for T` from the standard
// library. This is awfully annoying, so keep an eye out for specialization, which should make that sort of overlap
// possible: https://github.com/rust-lang/rust/issues/31844. In the meantime, this macro keeps things DRY...
macro_rules! convert_impls {
    ($($kind:ty),+ $(,)?) => {
        $(
            impl<'a, 'p> From<Modification<$kind>> for AnyModification<'a, 'p> {
                fn from(value: Modification<$kind>) -> Self {
                    value.convert()
                }
            }
        )+
    };
}

// NOTE: This should contain all modification kinds *except* for `AnyMod`, since that's the one that triggers overlap
// with the core `impl From<T> for T` blanket impl
convert_impls!(NamedMod<'a, 'p>, OffsetMod<'a>);

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
    fn monoisotopic_mass(&self) -> MonoisotopicMass {
        dispatch!(self, monoisotopic_mass)
    }

    fn average_mass(&self) -> AverageMass {
        dispatch!(self, average_mass)
    }
}

impl Charged for AnyMod<'_, '_> {
    fn charge(&self) -> Charge {
        dispatch!(self, charge)
    }
}

#[cfg(test)]
mod tests {
    use std::panic;

    use rust_decimal_macros::dec;
    use std::sync::LazyLock;

    use crate::{
        AtomicDatabase, AverageMz, ChargedParticle, MonoisotopicMz,
        testing_tools::assert_miette_snapshot,
    };

    use super::*;

    static ATOMIC_DB: LazyLock<AtomicDatabase> = LazyLock::new(AtomicDatabase::default);

    static H2O: LazyLock<ChemicalComposition> =
        LazyLock::new(|| ChemicalComposition::new(&ATOMIC_DB, "H2O").unwrap());
    static CA: LazyLock<ChemicalComposition> =
        LazyLock::new(|| ChemicalComposition::new(&ATOMIC_DB, "Ca-2e").unwrap());

    static POLYMER_DB: LazyLock<PolymerDatabase> = LazyLock::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "test_polymer_database.kdl",
            include_str!("../../tests/data/polymer_database.kdl"),
        )
        .unwrap()
    });

    #[test]
    fn errors() {
        let magnesium = AnyMod::named(&POLYMER_DB, "Mg");
        assert_miette_snapshot!(magnesium);
        let potassium = AnyMod::named(&POLYMER_DB, "K");
        assert_miette_snapshot!(potassium);
    }

    #[test]
    fn from_impls() {
        let named_mod = NamedMod::new(&POLYMER_DB, "Am").unwrap();
        let named_any_mod: AnyMod = named_mod.clone().into();
        let named_any_modification: AnyModification = named_mod.clone().into();
        let named_any_any_modification: AnyModification = named_any_mod.clone().into();
        assert_eq!(
            named_mod.monoisotopic_mass(),
            named_any_mod.monoisotopic_mass()
        );
        assert_eq!(
            named_mod.monoisotopic_mass(),
            named_any_modification.monoisotopic_mass()
        );
        assert_eq!(
            named_mod.monoisotopic_mass(),
            named_any_any_modification.monoisotopic_mass()
        );

        let offset_mod = OffsetMod::new(OffsetKind::Add, H2O.clone());
        let offset_any_mod: AnyMod = offset_mod.clone().into();
        let offset_any_modification: AnyModification = offset_mod.clone().into();
        let offset_any_any_modification: AnyModification = offset_any_mod.clone().into();
        assert_eq!(
            offset_mod.monoisotopic_mass(),
            offset_any_mod.monoisotopic_mass()
        );
        assert_eq!(
            offset_mod.monoisotopic_mass(),
            offset_any_modification.monoisotopic_mass()
        );
        assert_eq!(
            offset_mod.monoisotopic_mass(),
            offset_any_any_modification.monoisotopic_mass()
        );

        let named_modification = Modification::new(
            Count::new(2).unwrap(),
            NamedMod::new(&POLYMER_DB, "Am").unwrap(),
        );
        let named_any_modification: AnyModification = named_modification.clone().into();
        assert_eq!(
            named_modification.monoisotopic_mass(),
            named_any_modification.monoisotopic_mass()
        );

        let offset_modification = Modification::new(
            Count::new(3).unwrap(),
            OffsetMod::new(OffsetKind::Remove, H2O.clone()),
        );
        let offset_any_modification: AnyModification = offset_modification.clone().into();
        assert_eq!(
            offset_modification.monoisotopic_mass(),
            offset_any_modification.monoisotopic_mass()
        );
    }

    #[test]
    fn derive_more_is_variant() {
        let amidation = AnyMod::named(&POLYMER_DB, "Am").unwrap();
        let water_gained = AnyMod::offset(OffsetKind::Add, H2O.clone());

        assert!(amidation.is_named());
        assert!(!amidation.is_offset());

        assert!(!water_gained.is_named());
        assert!(water_gained.is_offset());
    }

    #[test]
    fn derive_more_unwrap() {
        macro_rules! assert_panics {
            ($expr:expr) => {
                assert!(panic::catch_unwind(|| $expr).is_err());
            };
        }

        let amidation = AnyMod::named(&POLYMER_DB, "Am").unwrap();
        let water_gained = AnyMod::offset(OffsetKind::Add, H2O.clone());

        amidation.clone().unwrap_named();
        assert_panics!(amidation.clone().unwrap_offset());

        assert_panics!(water_gained.clone().unwrap_named());
        water_gained.clone().unwrap_offset();
    }

    #[test]
    fn monoisotopic_mass() {
        // Masses checked against https://www.unimod.org/modifications_list.php
        let amidation = AnyMod::named(&POLYMER_DB, "Am").unwrap();
        assert_eq!(
            amidation.monoisotopic_mass(),
            MonoisotopicMass(dec!(-0.98401558291))
        );
        let acetylation = AnyMod::named(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(
            acetylation.monoisotopic_mass(),
            MonoisotopicMass(dec!(42.01056468403))
        );
        let deacetylation = AnyMod::named(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(
            deacetylation.monoisotopic_mass(),
            MonoisotopicMass(dec!(-42.01056468403))
        );
        let calcium = AnyMod::named(&POLYMER_DB, "Ca").unwrap();
        assert_eq!(
            calcium.monoisotopic_mass(),
            MonoisotopicMass(dec!(38.954217236560870))
        );

        let water_gained = AnyMod::offset(OffsetKind::Add, H2O.clone());
        assert_eq!(
            water_gained.monoisotopic_mass(),
            MonoisotopicMass(dec!(18.01056468403))
        );
        let water_lost = AnyMod::offset(OffsetKind::Remove, H2O.clone());
        assert_eq!(
            water_lost.monoisotopic_mass(),
            MonoisotopicMass(dec!(-18.01056468403))
        );
        // Masses checked against https://bioportal.bioontology.org/ontologies/UBERON
        let ca_gained = AnyMod::offset(OffsetKind::Add, CA.clone());
        assert_eq!(
            ca_gained.monoisotopic_mass(),
            MonoisotopicMass(dec!(39.961493703181870))
        );
        let ca_lost = AnyMod::offset(OffsetKind::Remove, CA.clone());
        assert_eq!(
            ca_lost.monoisotopic_mass(),
            MonoisotopicMass(dec!(-39.961493703181870))
        );
    }

    #[test]
    fn average_mass() {
        // Masses checked against https://www.unimod.org/modifications_list.php
        let amidation = AnyMod::named(&POLYMER_DB, "Am").unwrap();
        assert_eq!(
            amidation.average_mass(),
            AverageMass(dec!(-0.98476095881670255))
        );
        let acetylation = AnyMod::named(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(
            acetylation.average_mass(),
            AverageMass(dec!(42.03675822590033060))
        );
        let deacetylation = AnyMod::named(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(
            deacetylation.average_mass(),
            AverageMass(dec!(-42.03675822590033060))
        );
        let calcium = AnyMod::named(&POLYMER_DB, "Ca").unwrap();
        assert_eq!(
            calcium.average_mass(),
            AverageMass(dec!(39.069648884578600))
        );

        let water_gained = AnyMod::offset(OffsetKind::Add, H2O.clone());
        assert_eq!(
            water_gained.average_mass(),
            AverageMass(dec!(18.01528643242983260))
        );
        let water_lost = AnyMod::offset(OffsetKind::Remove, H2O.clone());
        assert_eq!(
            water_lost.average_mass(),
            AverageMass(dec!(-18.01528643242983260))
        );
        // Masses checked against https://bioportal.bioontology.org/ontologies/UBERON
        let ca_gained = AnyMod::offset(OffsetKind::Add, CA.clone());
        assert_eq!(
            ca_gained.average_mass(),
            AverageMass(dec!(40.076925351199600))
        );
        let ca_lost = AnyMod::offset(OffsetKind::Remove, CA.clone());
        assert_eq!(
            ca_lost.average_mass(),
            AverageMass(dec!(-40.076925351199600))
        );
    }

    #[test]
    fn charge() {
        let amidation = AnyMod::named(&POLYMER_DB, "Am").unwrap();
        assert_eq!(amidation.charge(), Charge(0));
        let acetylation = AnyMod::named(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(acetylation.charge(), Charge(0));
        let deacetylation = AnyMod::named(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(deacetylation.charge(), Charge(0));
        let calcium = AnyMod::named(&POLYMER_DB, "Ca").unwrap();
        assert_eq!(calcium.charge(), Charge(1));

        let water_gained = AnyMod::offset(OffsetKind::Add, H2O.clone());
        assert_eq!(water_gained.charge(), Charge(0));
        let water_lost = AnyMod::offset(OffsetKind::Remove, H2O.clone());
        assert_eq!(water_lost.charge(), Charge(0));
        let ca_gained = AnyMod::offset(OffsetKind::Add, CA.clone());
        assert_eq!(ca_gained.charge(), Charge(2));
        let ca_lost = AnyMod::offset(OffsetKind::Remove, CA.clone());
        assert_eq!(ca_lost.charge(), Charge(-2));
    }

    #[test]
    fn monoisotopic_mz() {
        let amidation = AnyMod::named(&POLYMER_DB, "Am").unwrap();
        assert_eq!(amidation.monoisotopic_mz(), None);
        let acetylation = AnyMod::named(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(acetylation.monoisotopic_mz(), None);
        let deacetylation = AnyMod::named(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(deacetylation.monoisotopic_mz(), None);
        let calcium = AnyMod::named(&POLYMER_DB, "Ca").unwrap();
        assert_eq!(
            calcium.monoisotopic_mz(),
            Some(MonoisotopicMz(dec!(38.954217236560870)))
        );

        let water_gained = AnyMod::offset(OffsetKind::Add, H2O.clone());
        assert_eq!(water_gained.monoisotopic_mz(), None);
        let water_lost = AnyMod::offset(OffsetKind::Remove, H2O.clone());
        assert_eq!(water_lost.monoisotopic_mz(), None);
        let ca_gained = AnyMod::offset(OffsetKind::Add, CA.clone());
        assert_eq!(
            ca_gained.monoisotopic_mz(),
            Some(MonoisotopicMz(dec!(19.980746851590935)))
        );
        let ca_lost = AnyMod::offset(OffsetKind::Remove, CA.clone());
        assert_eq!(
            ca_lost.monoisotopic_mz(),
            Some(MonoisotopicMz(dec!(-19.980746851590935)))
        );
    }

    #[test]
    fn average_mz() {
        let amidation = AnyMod::named(&POLYMER_DB, "Am").unwrap();
        assert_eq!(amidation.average_mz(), None);
        let acetylation = AnyMod::named(&POLYMER_DB, "Ac").unwrap();
        assert_eq!(acetylation.average_mz(), None);
        let deacetylation = AnyMod::named(&POLYMER_DB, "DeAc").unwrap();
        assert_eq!(deacetylation.average_mz(), None);
        let calcium = AnyMod::named(&POLYMER_DB, "Ca").unwrap();
        assert_eq!(
            calcium.average_mz(),
            Some(AverageMz(dec!(39.069648884578600)))
        );

        let water_gained = AnyMod::offset(OffsetKind::Add, H2O.clone());
        assert_eq!(water_gained.average_mz(), None);
        let water_lost = AnyMod::offset(OffsetKind::Remove, H2O.clone());
        assert_eq!(water_lost.average_mz(), None);
        let ca_gained = AnyMod::offset(OffsetKind::Add, CA.clone());
        assert_eq!(
            ca_gained.average_mz(),
            Some(AverageMz(dec!(20.0384626755998)))
        );
        let ca_lost = AnyMod::offset(OffsetKind::Remove, CA.clone());
        assert_eq!(
            ca_lost.average_mz(),
            Some(AverageMz(dec!(-20.0384626755998)))
        );
    }
}
