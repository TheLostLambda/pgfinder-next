use ahash::{HashSet, HashSetExt};

use crate::{
    errors::PolychemError, AverageMass, Charge, Charged, FunctionalGroup, GroupState, Massive,
    ModificationId, MonoisotopicMass, Residue, Result,
};

use super::polymer_database::{PolymerDatabase, ResidueDescription};

impl<'a, 'p> Residue<'a, 'p> {
    pub fn new(db: &'p PolymerDatabase<'a>, abbr: impl AsRef<str>) -> Result<Self> {
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

    use crate::{testing_tools::assert_miette_snapshot, AtomicDatabase, ChargedParticle};

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

    static RESIDUES: Lazy<[Residue; 3]> =
        Lazy::new(|| ["A", "m", "X"].map(|abbr| Residue::new(&POLYMER_DB, abbr).unwrap()));

    #[test]
    fn errors() {
        let sucrose = Residue::new(&POLYMER_DB, "s");
        assert_miette_snapshot!(sucrose);
        let super_amino = Residue::new(&POLYMER_DB, "Sa");
        assert_miette_snapshot!(super_amino);
    }

    static N_TERMINAL: FunctionalGroup = FunctionalGroup::new("Amino", "N-Terminal");

    static C_TERMINAL: FunctionalGroup = FunctionalGroup::new("Carboxyl", "C-Terminal");

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
    fn group_state_mut() {
        let mut lysine = Residue::new(&POLYMER_DB, "K").unwrap();
        let mut alanine = Residue::new(&POLYMER_DB, "A").unwrap();
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
    fn monoisotopic_mass() {
        assert_eq!(
            RESIDUES.each_ref().map(Massive::monoisotopic_mass),
            [
                MonoisotopicMass(dec!(89.04767846918)),
                MonoisotopicMass(dec!(293.11106657336)),
                MonoisotopicMass(dec!(0))
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
                AverageMass(dec!(0))
            ]
        );
    }

    #[test]
    fn charge() {
        assert_eq!(
            RESIDUES.each_ref().map(Charged::charge),
            [Charge(0), Charge(0), Charge(0)]
        );
    }

    #[test]
    fn monoisotopic_mz() {
        assert_eq!(
            RESIDUES.each_ref().map(ChargedParticle::monoisotopic_mz),
            [None, None, None]
        );
    }

    #[test]
    fn average_mz() {
        assert_eq!(
            RESIDUES.each_ref().map(ChargedParticle::average_mz),
            [None, None, None]
        );
    }
}
