use std::iter;

use ahash::{HashSet, HashSetExt};

use crate::{
    errors::PolychemError, Charge, Charged, Count, FunctionalGroup, GroupState, Id, Mass, Massive,
    Modification, ModificationId, OffsetMod, Residue, Result,
};

use super::{
    errors::OffsetMultiplierError,
    polymer_database::{PolymerDatabase, ResidueDescription},
};

impl<'a, 'p> Residue<'a, 'p> {
    pub fn new(db: &'p PolymerDatabase<'a>, abbr: impl AsRef<str>, id: Id) -> Result<Self> {
        let abbr = abbr.as_ref();
        let (
            abbr,
            ResidueDescription {
                name,
                composition,
                functional_groups,
            },
        ) = db
            .residues
            .get_key_value(abbr)
            .ok_or_else(|| PolychemError::residue_lookup(abbr))?;
        let functional_groups = functional_groups
            .iter()
            .map(|fg| (fg.into(), GroupState::default()))
            .collect();
        Ok(Self {
            abbr,
            name,
            composition,
            functional_groups,
            offset_modifications: HashSet::new(),
        })
    }

    #[must_use]
    pub const fn abbr(&self) -> &'p str {
        self.abbr
    }

    #[must_use]
    pub const fn name(&self) -> &'p str {
        self.name
    }

    pub fn group_state(&self, functional_group: &FunctionalGroup<'p>) -> Result<&GroupState> {
        self.functional_groups.get(functional_group).ok_or_else(|| {
            PolychemError::group_lookup(*functional_group, self.name, self.abbr).into()
        })
    }

    // FIXME: This cannot be public — it breaks the polymerizer free-group index — think about if residue should have
    // any public methods... I think they should, since they can be looked up by ID in the polymerizer!
    pub(crate) fn group_state_mut(
        &mut self,
        functional_group: &FunctionalGroup<'p>,
    ) -> Result<&mut GroupState> {
        self.functional_groups
            .get_mut(functional_group)
            .ok_or_else(|| {
                PolychemError::group_lookup(*functional_group, self.name, self.abbr).into()
            })
    }

    pub fn offset_modifications(&'a self) -> impl Iterator<Item = ModificationId> + 'a {
        self.offset_modifications.iter().copied()
    }

    // TODO: Write a `named_modifications()` equivalent!
}

impl Massive for Residue<'_, '_> {
    fn monoisotopic_mass(&self) -> Mass {
        self.composition.monoisotopic_mass()
    }

    fn average_mass(&self) -> Mass {
        self.composition.average_mass()
    }
}

impl Charged for Residue<'_, '_> {
    fn charge(&self) -> Charge {
        self.composition.charge()
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::{
        testing_tools::assert_miette_snapshot, AtomicDatabase, Bond, BondTarget, ChargedParticle,
        NamedMod, OffsetKind,
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
    fn errors() {
        let sucrose = Residue::new(&POLYMER_DB, "s", 0);
        assert_miette_snapshot!(sucrose);
        let super_amino = Residue::new(&POLYMER_DB, "Sa", 0);
        assert_miette_snapshot!(super_amino);
    }

    #[test]
    fn id() {
        let alanine = Residue::new(&POLYMER_DB, "A", 0).unwrap();
        assert_eq!(alanine.id(), 0);
        let alanine = Residue::new(&POLYMER_DB, "A", 42).unwrap();
        assert_eq!(alanine.id(), 42);
        let alanine = Residue::new(&POLYMER_DB, "A", 1).unwrap();
        assert_eq!(alanine.id(), 1);
    }

    static N_TERMINAL: FunctionalGroup = FunctionalGroup::new("Amino", "N-Terminal");

    static C_TERMINAL: FunctionalGroup = FunctionalGroup::new("Carboxyl", "C-Terminal");

    #[test]
    fn group_state() {
        let lysine = Residue::new(&POLYMER_DB, "K", 0).unwrap();
        let alanine = Residue::new(&POLYMER_DB, "A", 0).unwrap();
        let sidechain_amino = FunctionalGroup::new("Amino", "Sidechain");

        // Sucessfully lookup functional groups that exist
        assert!(lysine.group_state(&N_TERMINAL).is_ok());
        assert!(lysine.group_state(&C_TERMINAL).is_ok());
        assert!(lysine.group_state(&sidechain_amino).is_ok());

        assert!(alanine.group_state(&N_TERMINAL).is_ok());
        assert!(alanine.group_state(&C_TERMINAL).is_ok());

        // Fail to lookup functional groups that don't exist
        assert_miette_snapshot!(alanine.group_state(&sidechain_amino));
    }

    #[test]
    fn group_state_mut() {
        let mut lysine = Residue::new(&POLYMER_DB, "K", 0).unwrap();
        let mut alanine = Residue::new(&POLYMER_DB, "A", 0).unwrap();
        let sidechain_amino = FunctionalGroup::new("Amino", "Sidechain");

        // Sucessfully lookup functional groups that exist
        assert!(lysine.group_state_mut(&N_TERMINAL).is_ok());
        assert!(lysine.group_state_mut(&C_TERMINAL).is_ok());
        assert!(lysine.group_state_mut(&sidechain_amino).is_ok());

        assert!(alanine.group_state_mut(&N_TERMINAL).is_ok());
        assert!(alanine.group_state_mut(&C_TERMINAL).is_ok());

        // Fail to lookup functional groups that don't exist
        assert_miette_snapshot!(alanine.group_state_mut(&sidechain_amino));
    }

    #[test]
    fn add_offset() {
        let mut alanine = Residue::new(&POLYMER_DB, "A", 0).unwrap();
        let water_loss = OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap();
        let proton = |kind| OffsetMod::new(&ATOMIC_DB, kind, "p").unwrap();
        macro_rules! assert_offset_names_and_counts {
            ($input:ident, $output:expr) => {
                let mut sorted_offsets: Vec<_> = $input
                    .offset_modifications()
                    .map(|m| m.to_string())
                    .collect();
                sorted_offsets.sort_unstable();
                assert_eq!(sorted_offsets, $output);
            };
        }

        assert_eq!(alanine.add_offset(water_loss.clone()).unwrap(), -1);
        assert_offset_names_and_counts!(alanine, vec!["-H2O"]);
        assert_eq!(alanine.add_offset(water_loss.clone()).unwrap(), -2);
        assert_offset_names_and_counts!(alanine, vec!["-2xH2O"]);

        assert_eq!(alanine.add_offsets(proton(OffsetKind::Add), 2).unwrap(), 2);
        assert_offset_names_and_counts!(alanine, vec!["+2xp", "-2xH2O"]);
        assert_eq!(alanine.add_offsets(proton(OffsetKind::Add), 2).unwrap(), 4);
        assert_offset_names_and_counts!(alanine, vec!["+4xp", "-2xH2O"]);

        assert_eq!(
            alanine.add_offsets(proton(OffsetKind::Remove), 2).unwrap(),
            2
        );
        assert_offset_names_and_counts!(alanine, vec!["+2xp", "-2xH2O"]);
        assert_eq!(
            alanine.add_offsets(proton(OffsetKind::Remove), 2).unwrap(),
            0
        );
        assert_offset_names_and_counts!(alanine, vec!["-2xH2O"]);

        assert_eq!(
            alanine
                .add_offsets(proton(OffsetKind::Add), Count::MAX)
                .unwrap(),
            Count::MAX.into()
        );
        assert_offset_names_and_counts!(alanine, vec!["+4294967295xp", "-2xH2O"]);
        assert_miette_snapshot!(alanine.add_offsets(proton(OffsetKind::Add), Count::MAX));

        assert_miette_snapshot!(alanine.add_offsets(water_loss, 0));
    }

    static RESIDUE_SERIES: Lazy<Vec<Residue<'static, 'static>>> = Lazy::new(|| {
        let mut snapshots = Vec::new();

        let mut alanine = Residue::new(&POLYMER_DB, "A", 0).unwrap();
        snapshots.push(alanine.clone());

        // Add a water-loss offset modification
        let water_loss = OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap();
        alanine.add_offset(water_loss).unwrap();
        snapshots.push(alanine.clone());

        // Add an amidation named modification to the C-terminal
        let amidation = GroupState::Modified(NamedMod::new(&POLYMER_DB, "Am").unwrap());
        *alanine.group_state_mut(&C_TERMINAL).unwrap() = amidation;
        snapshots.push(alanine.clone());

        // Add an amidation named modification to the N-terminal (ignoring that that's impossible)
        *alanine.group_state_mut(&N_TERMINAL).unwrap() = amidation;
        snapshots.push(alanine.clone());

        // Out of functional groups, so adding more amidations changes nothing
        *alanine.group_state_mut(&N_TERMINAL).unwrap() = amidation;
        *alanine.group_state_mut(&C_TERMINAL).unwrap() = amidation;
        snapshots.push(alanine.clone());

        // But they can be replaced with bonds
        let peptide_bond =
            Bond::new(&POLYMER_DB, "Peptide", BondTarget::new(0, C_TERMINAL)).unwrap();
        *alanine.group_state_mut(&N_TERMINAL).unwrap() = GroupState::Donor(peptide_bond);
        snapshots.push(alanine.clone());

        // Residues can be protonated
        let proton = OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "p").unwrap();
        alanine.add_offsets(proton, 2).unwrap();
        snapshots.push(alanine.clone());

        // Or can form other adducts
        let ca = OffsetMod::new(&ATOMIC_DB, OffsetKind::Add, "Ca-2e").unwrap();
        alanine.add_offset(ca).unwrap();
        snapshots.push(alanine.clone());

        // Removing the two protons...
        let anti_proton = OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "p").unwrap();
        alanine.add_offsets(anti_proton, 2).unwrap();
        snapshots.push(alanine.clone());

        snapshots
    });

    #[test]
    fn monoisotopic_mass() {
        assert_eq!(
            RESIDUE_SERIES
                .iter()
                .map(Massive::monoisotopic_mass)
                .collect_vec(),
            vec![
                dec!(89.04767846918),
                dec!(71.03711378515),
                dec!(70.05309820224),
                dec!(69.06908261933),
                dec!(69.06908261933),
                dec!(52.04253351821),
                dec!(54.057086451452),
                dec!(94.018580154633870),
                dec!(92.004027221391870),
            ]
        );
    }

    #[test]
    fn average_mass() {
        assert_eq!(
            RESIDUE_SERIES
                .iter()
                .map(Massive::average_mass)
                .collect_vec(),
            vec![
                dec!(89.09330602867854225),
                dec!(71.07801959624870965),
                dec!(70.09325863743200710),
                dec!(69.10849767861530455),
                dec!(69.10849767861530455),
                dec!(52.07797220500217450),
                dec!(54.09252513824417450),
                dec!(94.16945048944377450),
                dec!(92.15489755620177450),
            ]
        );
    }

    #[test]
    fn charge() {
        assert_eq!(
            RESIDUE_SERIES.iter().map(Charged::charge).collect_vec(),
            vec![0, 0, 0, 0, 0, 0, 2, 4, 2]
        );
    }

    #[test]
    fn monoisotopic_mz() {
        assert_eq!(
            RESIDUE_SERIES
                .iter()
                .map(ChargedParticle::monoisotopic_mz)
                .collect_vec(),
            vec![
                None,
                None,
                None,
                None,
                None,
                None,
                Some(dec!(27.028543225726)),
                Some(dec!(23.50464503865846750)),
                Some(dec!(46.002013610695935))
            ]
        );
    }

    #[test]
    fn average_mz() {
        assert_eq!(
            RESIDUE_SERIES
                .iter()
                .map(ChargedParticle::average_mz)
                .collect_vec(),
            vec![
                None,
                None,
                None,
                None,
                None,
                None,
                Some(dec!(27.04626256912208725)),
                Some(dec!(23.5423626223609436250)),
                Some(dec!(46.07744877810088725))
            ]
        );
    }
}
