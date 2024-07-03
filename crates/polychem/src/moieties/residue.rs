use ahash::{HashSet, HashSetExt};

use crate::{
    errors::PolychemError, AverageMass, Charge, Charged, FunctionalGroup, GroupState, Massive,
    ModificationId, MonoisotopicMass, Residue, Result,
};

use super::polymer_database::{PolymerDatabase, ResidueDescription};

impl<'a, 'p> Residue<'a, 'p> {
    pub(crate) fn new(db: &'p PolymerDatabase<'a>, abbr: impl AsRef<str>) -> Result<Self> {
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
        let offset_modifications = HashSet::new();
        Ok(Self {
            abbr,
            name,
            composition,
            functional_groups,
            offset_modifications,
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

    pub fn functional_groups(&self) -> impl Iterator<Item = (&FunctionalGroup<'p>, &GroupState)> {
        self.functional_groups.iter()
    }

    // FIXME: Remove the '_ after Rust 2024 is released
    pub fn offset_modifications(&self) -> impl Iterator<Item = ModificationId> + '_ {
        self.offset_modifications.iter().copied()
    }

    pub fn group_state(&self, functional_group: &FunctionalGroup<'p>) -> Result<&GroupState> {
        self.functional_groups.get(functional_group).ok_or_else(|| {
            PolychemError::group_lookup(*functional_group, self.name, self.abbr).into()
        })
    }

    // NOTE: This cannot be public API — it would let users invalidate the polymerizer group-state index / cache
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

    // NOTE: This cannot be public API — it would let users invalidate `Polymer`s mapping of `ModificationId`s
    pub(crate) fn offset_modifications_mut(&mut self) -> &mut HashSet<ModificationId> {
        &mut self.offset_modifications
    }
}

impl Massive for Residue<'_, '_> {
    fn monoisotopic_mass(&self) -> MonoisotopicMass {
        self.composition.monoisotopic_mass()
    }

    fn average_mass(&self) -> AverageMass {
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
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::{
        testing_tools::assert_miette_snapshot, AtomicDatabase, AverageMz, BondId, ChargedParticle,
        MonoisotopicMz,
    };

    use super::*;

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);

    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "test_polymer_database.kdl",
            include_str!("../../tests/data/polymer_database.kdl"),
        )
        .unwrap()
    });

    static RESIDUES: Lazy<[Residue; 4]> =
        Lazy::new(|| ["A", "m", "X", "K2+"].map(|abbr| Residue::new(&POLYMER_DB, abbr).unwrap()));

    #[test]
    fn errors() {
        let sucrose = Residue::new(&POLYMER_DB, "s");
        assert_miette_snapshot!(sucrose);
        let super_amino = Residue::new(&POLYMER_DB, "Sa");
        assert_miette_snapshot!(super_amino);
    }

    #[test]
    fn names_and_abbrs() {
        assert_eq!(
            RESIDUES.each_ref().map(|m| (m.abbr(), m.name())),
            [
                ("A", "Alanine"),
                ("m", "N-Acetylmuramic Acid"),
                ("X", "Unknown Amino Acid"),
                ("K2+", "Lysine 2+")
            ]
        );
    }

    static N_TERMINAL: FunctionalGroup = FunctionalGroup::new("Amino", "N-Terminal");

    static C_TERMINAL: FunctionalGroup = FunctionalGroup::new("Carboxyl", "C-Terminal");

    #[test]
    fn functional_groups() {
        let mut lysine = Residue::new(&POLYMER_DB, "K").unwrap();
        lysine
            .functional_groups
            .insert(N_TERMINAL, GroupState::Acceptor(BondId(0)));
        lysine
            .functional_groups
            .insert(C_TERMINAL, GroupState::Modified(ModificationId(0)));

        // Overwrite an earlier assignment
        lysine
            .functional_groups
            .insert(N_TERMINAL, GroupState::Donor(BondId(0)));

        let mut groups: Vec<_> = lysine.functional_groups().collect();
        groups.sort_unstable();
        assert_eq!(
            groups,
            vec![
                (
                    &FunctionalGroup {
                        name: "Amino",
                        location: "N-Terminal"
                    },
                    &GroupState::Donor(BondId(0))
                ),
                (
                    &FunctionalGroup {
                        name: "Amino",
                        location: "Sidechain"
                    },
                    &GroupState::Free
                ),
                (
                    &FunctionalGroup {
                        name: "Carboxyl",
                        location: "C-Terminal"
                    },
                    &GroupState::Modified(ModificationId(0))
                ),
            ]
        );
    }

    #[test]
    fn offset_modifications() {
        let mut alanine = Residue::new(&POLYMER_DB, "A").unwrap();
        alanine.offset_modifications.insert(ModificationId(0));
        alanine.offset_modifications.insert(ModificationId(1));

        // Try to insert a duplicate modification
        alanine.offset_modifications.insert(ModificationId(0));

        let mut mods: Vec<_> = alanine.offset_modifications().collect();
        mods.sort_unstable();
        assert_eq!(mods, vec![ModificationId(0), ModificationId(1)]);
    }

    #[test]
    fn group_state() {
        let lysine = Residue::new(&POLYMER_DB, "K").unwrap();
        let alanine = Residue::new(&POLYMER_DB, "A").unwrap();
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
    fn monoisotopic_mass() {
        assert_eq!(
            RESIDUES.each_ref().map(Massive::monoisotopic_mass),
            [
                MonoisotopicMass(dec!(89.04767846918)),
                MonoisotopicMass(dec!(293.11106657336)),
                MonoisotopicMass(dec!(0)),
                MonoisotopicMass(dec!(146.104430568002)),
            ]
        );
    }

    #[test]
    fn average_mass() {
        assert_eq!(
            RESIDUES.each_ref().map(Massive::average_mass),
            [
                AverageMass(dec!(89.09330602867854225)),
                AverageMass(dec!(293.27091179713952985)),
                AverageMass(dec!(0)),
                AverageMass(dec!(146.18647363385097400)),
            ]
        );
    }

    #[test]
    fn charge() {
        assert_eq!(
            RESIDUES.each_ref().map(Charged::charge),
            [Charge(0), Charge(0), Charge(0), Charge(2)]
        );
    }

    #[test]
    fn monoisotopic_mz() {
        assert_eq!(
            RESIDUES.each_ref().map(ChargedParticle::monoisotopic_mz),
            [
                None,
                None,
                None,
                Some(MonoisotopicMz(dec!(73.052215284001)))
            ]
        );
    }

    #[test]
    fn average_mz() {
        assert_eq!(
            RESIDUES.each_ref().map(ChargedParticle::average_mz),
            [
                None,
                None,
                None,
                Some(AverageMz(dec!(73.09323681692548700)))
            ]
        );
    }
}
