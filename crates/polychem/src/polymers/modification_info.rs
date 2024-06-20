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
    use super::*;
    use once_cell::sync::Lazy;

    use crate::{
        AtomicDatabase, ChemicalComposition, Count, FunctionalGroup, Modification, NamedMod,
        OffsetKind, OffsetMod, PolymerDatabase, ResidueGroup, ResidueId,
    };

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);
    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "test_polymer_database.kdl",
            include_str!("../../tests/data/polymer_database.kdl"),
        )
        .unwrap()
    });

    static NAMED: Lazy<ModificationInfo> = Lazy::new(|| {
        let modification = NamedMod::new(&POLYMER_DB, "Am").unwrap();
        let residue_group =
            ResidueGroup(ResidueId(0), FunctionalGroup::new("Carboxyl", "Sidechain"));
        ModificationInfo::Named(modification, residue_group)
    });

    static OFFSET: Lazy<ModificationInfo> = Lazy::new(|| {
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

    static UNLOCALIZED: Lazy<[ModificationInfo; 2]> = Lazy::new(|| {
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
    fn display_modification_info() {
        assert_eq!(NAMED.to_string(), "a named".to_owned());

        assert_eq!(OFFSET.to_string(), "an offset".to_owned());

        for unlocalized in &*UNLOCALIZED {
            assert_eq!(unlocalized.to_string(), "an unlocalized".to_owned());
        }
    }
}
