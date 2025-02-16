use crate::ModificationInfo;
use std::fmt::{self, Display, Formatter};

impl Display for ModificationInfo<'_, '_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                ModificationInfo::Named(..) => "a named",
                ModificationInfo::Offset(..) => "an offset",
                ModificationInfo::Unlocalized(..) => "an unlocalized",
            }
        )
    }
}

#[cfg(test)]
mod tests {
    use std::{panic, sync::LazyLock};

    use super::*;

    use crate::{
        AtomicDatabase, ChemicalComposition, Count, FunctionalGroup, Modification, NamedMod,
        OffsetKind, OffsetMod, PolymerDatabase, ResidueGroup, ResidueId,
    };

    static ATOMIC_DB: LazyLock<AtomicDatabase> = LazyLock::new(AtomicDatabase::default);
    static POLYMER_DB: LazyLock<PolymerDatabase> = LazyLock::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "test_polymer_database.kdl",
            include_str!("../../tests/data/polymer_database.kdl"),
        )
        .unwrap()
    });

    static NAMED: LazyLock<ModificationInfo> = LazyLock::new(|| {
        let modification = NamedMod::new(&POLYMER_DB, "Am").unwrap();
        let residue_group =
            ResidueGroup(ResidueId(0), FunctionalGroup::new("Carboxyl", "Sidechain"));
        ModificationInfo::Named(modification, residue_group)
    });

    static OFFSET: LazyLock<ModificationInfo> = LazyLock::new(|| {
        let modification = Modification::new(
            Count::new(2).unwrap(),
            OffsetMod::new(
                OffsetKind::Remove,
                ChemicalComposition::new(&ATOMIC_DB, "H2O").unwrap(),
            ),
        );
        let residue = ResidueId(0);
        ModificationInfo::Offset(modification, residue)
    });

    static UNLOCALIZED: LazyLock<[ModificationInfo; 2]> = LazyLock::new(|| {
        let ModificationInfo::Named(named, ..) = NAMED.clone() else {
            unreachable!()
        };
        let ModificationInfo::Offset(offset, ..) = OFFSET.clone() else {
            unreachable!()
        };
        [named.into(), offset.into()].map(ModificationInfo::Unlocalized)
    });

    #[test]
    fn derive_more_is_variant() {
        assert!(NAMED.is_named());
        assert!(!NAMED.is_offset());
        assert!(!NAMED.is_unlocalized());

        assert!(!OFFSET.is_named());
        assert!(OFFSET.is_offset());
        assert!(!OFFSET.is_unlocalized());

        for unlocalized in &*UNLOCALIZED {
            assert!(!unlocalized.is_named());
            assert!(!unlocalized.is_offset());
            assert!(unlocalized.is_unlocalized());
        }
    }

    #[test]
    fn derive_more_unwrap() {
        macro_rules! assert_panics {
            ($expr:expr) => {
                assert!(panic::catch_unwind(|| $expr).is_err());
            };
        }

        NAMED.clone().unwrap_named();
        assert_panics!(NAMED.clone().unwrap_offset());
        assert_panics!(NAMED.clone().unwrap_unlocalized());

        assert_panics!(OFFSET.clone().unwrap_named());
        OFFSET.clone().unwrap_offset();
        assert_panics!(OFFSET.clone().unwrap_unlocalized());

        for unlocalized in &*UNLOCALIZED {
            assert_panics!(unlocalized.clone().unwrap_named());
            assert_panics!(unlocalized.clone().unwrap_offset());
            unlocalized.clone().unwrap_unlocalized();
        }
    }

    #[test]
    fn display_modification_info() {
        assert_eq!(NAMED.to_string(), "a named".to_owned());

        assert_eq!(OFFSET.to_string(), "an offset".to_owned());

        for unlocalized in &*UNLOCALIZED {
            assert_eq!(unlocalized.to_string(), "an unlocalized".to_owned());
        }
    }
}
