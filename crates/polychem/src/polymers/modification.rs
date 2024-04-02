use std::fmt::{self, Display, Formatter};

use rust_decimal::prelude::Decimal;

use crate::{AnyMod, Charge, Charged, Count, Massive, Modification, NamedMod, Offset, SignedCount};

impl<K> Modification<K> {
    // FIXME: Maybe make this `new_with_multiplier`, and make `new(k)` = `new_with_multiplier(1, k)` — depends on if
    // this code ends up as public API (usable outside of the crate)!
    pub const fn new(multiplier: Count, kind: K) -> Self {
        Self { multiplier, kind }
    }

    // NOTE: This can't exist as a From impl since it overlaps with `impl From<T> for T` in the standard library. This
    // is awfully annoying, so keep an eye out for specialization, which should make that sort of overlap possible:
    // https://github.com/rust-lang/rust/issues/31844
    pub fn into<K2: From<K>>(self) -> Modification<K2> {
        let Self { multiplier, kind } = self;
        Modification::new(multiplier, kind.into())
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
        SignedCount::from(self.multiplier) * self.kind.charge()
    }
}

fn display_offset_modification(
    f: &mut Formatter<'_>,
    multiplier: Count,
    offset: &Offset<impl Display>,
) -> fmt::Result {
    let Offset { kind, composition } = offset;
    if multiplier > 1 {
        write!(f, "{kind}{multiplier}x{composition}")
    } else {
        write!(f, "{kind}{composition}")
    }
}

fn display_named_modification(
    f: &mut Formatter<'_>,
    multiplier: Count,
    named: &NamedMod,
) -> fmt::Result {
    let NamedMod { abbr, .. } = named;
    if multiplier > 1 {
        write!(f, "{multiplier}x{abbr}")
    } else {
        write!(f, "{abbr}")
    }
}

impl<C: Display> Display for Modification<Offset<C>> {
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
        atoms::atomic_database::AtomicDatabase, polymers::polymer_database::PolymerDatabase, Mz,
        OffsetKind, OffsetMod,
    };

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
    fn modifcation_display() {
        let water_gain = Modification::new(
            1,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap(),
        );
        assert_eq!(water_gain.to_string(), "+H2O");
        let water_loss = Modification::new(
            1,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap(),
        );
        assert_eq!(water_loss.to_string(), "-H2O");
        let double_water_gain = Modification::new(
            2,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap(),
        );
        assert_eq!(double_water_gain.to_string(), "+2xH2O");
        let double_water_loss = Modification::new(
            2,
            OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap(),
        );
        assert_eq!(double_water_loss.to_string(), "-2xH2O");

        let amidation = Modification::new(1, NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.to_string(), "Am");
        let double_amidation = Modification::new(2, NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(double_amidation.to_string(), "2xAm");

        let water_gain = Modification::new(
            1,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap(),
        );
        assert_eq!(water_gain.to_string(), "+H2O");
        let water_loss = Modification::new(
            1,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap(),
        );
        assert_eq!(water_loss.to_string(), "-H2O");
        let double_water_gain = Modification::new(
            2,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Add, "H2O").unwrap(),
        );
        assert_eq!(double_water_gain.to_string(), "+2xH2O");
        let double_water_loss = Modification::new(
            2,
            AnyMod::offset(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap(),
        );
        assert_eq!(double_water_loss.to_string(), "-2xH2O");

        let amidation = Modification::new(1, AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.to_string(), "Am");
        let double_amidation = Modification::new(2, AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(double_amidation.to_string(), "2xAm");
    }

    #[test]
    fn monoisotopic_mass() {
        // Masses checked against https://www.unimod.org/modifications_list.php
        let amidation = Modification::new(1, NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.monoisotopic_mass(), dec!(-0.98401558291));
        let acetylation = Modification::new(2, NamedMod::new(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.monoisotopic_mass(), dec!(84.02112936806));
        let deacetylation = Modification::new(3, NamedMod::new(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.monoisotopic_mass(), dec!(-126.03169405209));
        let calcium = Modification::new(4, NamedMod::new(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(calcium.monoisotopic_mass(), dec!(155.81686894624348));

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
        let amidation = Modification::new(4, AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.monoisotopic_mass(), dec!(-3.93606233164));
        let acetylation = Modification::new(3, AnyMod::named(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.monoisotopic_mass(), dec!(126.03169405209));
        let deacetylation = Modification::new(2, AnyMod::named(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.monoisotopic_mass(), dec!(-84.02112936806));
        let calcium = Modification::new(1, AnyMod::named(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(calcium.monoisotopic_mass(), dec!(38.954217236560870));

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
        let calcium = Modification::new(4, NamedMod::new(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(calcium.average_mass(), dec!(156.278595538314400));

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
        let amidation = Modification::new(4, AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.average_mass(), dec!(-3.93904383526681020));
        let acetylation = Modification::new(3, AnyMod::named(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.average_mass(), dec!(126.11027467770099180));
        let deacetylation = Modification::new(2, AnyMod::named(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.average_mass(), dec!(-84.07351645180066120));
        let calcium = Modification::new(1, AnyMod::named(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(calcium.average_mass(), dec!(39.069648884578600));

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
        let amidation = Modification::new(1, NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.charge(), 0);
        let acetylation = Modification::new(2, NamedMod::new(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.charge(), 0);
        let deacetylation = Modification::new(3, NamedMod::new(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.charge(), 0);
        let calcium = Modification::new(4, NamedMod::new(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(calcium.charge(), 4);

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

        let amidation = Modification::new(4, AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.charge(), 0);
        let acetylation = Modification::new(3, AnyMod::named(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.charge(), 0);
        let deacetylation = Modification::new(2, AnyMod::named(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.charge(), 0);
        let calcium = Modification::new(1, AnyMod::named(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(calcium.charge(), 1);

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
        let amidation = Modification::new(1, NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.monoisotopic_mz(), None);
        let acetylation = Modification::new(2, NamedMod::new(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.monoisotopic_mz(), None);
        let deacetylation = Modification::new(3, NamedMod::new(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.monoisotopic_mz(), None);
        let calcium = Modification::new(4, NamedMod::new(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(calcium.monoisotopic_mz(), Some(dec!(38.954217236560870)));

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

        let amidation = Modification::new(4, AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.monoisotopic_mz(), None);
        let acetylation = Modification::new(3, AnyMod::named(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.monoisotopic_mz(), None);
        let deacetylation = Modification::new(2, AnyMod::named(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.monoisotopic_mz(), None);
        let calcium = Modification::new(1, AnyMod::named(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(calcium.monoisotopic_mz(), Some(dec!(38.954217236560870)));

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
        let amidation = Modification::new(1, NamedMod::new(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.average_mz(), None);
        let acetylation = Modification::new(2, NamedMod::new(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.average_mz(), None);
        let deacetylation = Modification::new(3, NamedMod::new(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.average_mz(), None);
        let calcium = Modification::new(4, NamedMod::new(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(calcium.average_mz(), Some(dec!(39.069648884578600)));

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

        let amidation = Modification::new(4, AnyMod::named(&POLYMER_DB, "Am").unwrap());
        assert_eq!(amidation.average_mz(), None);
        let acetylation = Modification::new(3, AnyMod::named(&POLYMER_DB, "Ac").unwrap());
        assert_eq!(acetylation.average_mz(), None);
        let deacetylation = Modification::new(2, AnyMod::named(&POLYMER_DB, "DeAc").unwrap());
        assert_eq!(deacetylation.average_mz(), None);
        let calcium = Modification::new(1, AnyMod::named(&POLYMER_DB, "Ca").unwrap());
        assert_eq!(calcium.average_mz(), Some(dec!(39.069648884578600)));

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
