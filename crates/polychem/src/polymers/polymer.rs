use core::slice;

use crate::{
    errors::PolychemError,
    moieties::{polymer_database::BondDescription, target::Target},
    AtomicDatabase, AverageMass, Bond, BondId, BondInfo, Charge, Charged, FunctionalGroup,
    GroupState, Massive, ModificationInfo, MonoisotopicMass, Polymer, PolymerDatabase, Residue,
    ResidueId, Result,
};

use super::{errors::FindFreeGroupsError, polymerizer_state::PolymerizerState};

impl<'a, 'p> Polymer<'a, 'p> {
    #[must_use]
    pub const fn atomic_db(&self) -> &'a AtomicDatabase {
        self.polymerizer_state.polymerizer.atomic_db()
    }

    #[must_use]
    pub const fn polymer_db(&self) -> &'p PolymerDatabase<'a> {
        self.polymerizer_state.polymerizer.polymer_db()
    }

    pub fn new_residue(&mut self, abbr: impl AsRef<str>) -> Result<ResidueId> {
        let residue = Residue::new(self.polymer_db(), abbr)?;
        let id = ResidueId(self.polymerizer_state.next_id());
        self.polymerizer_state.index_residue_groups(id, &residue);
        self.residues.insert(id, residue);
        Ok(id)
    }

    #[must_use]
    pub fn residue(&self, id: ResidueId) -> Option<&Residue<'a, 'p>> {
        self.residues.get(&id)
    }

    pub fn bond_residues(
        &mut self,
        abbr: impl AsRef<str>,
        donor: ResidueId,
        acceptor: ResidueId,
    ) -> Result<BondId> {
        // NOTE: Perhaps this explicit type annotation for `state` won't be necessary after the compiler improves its
        // type inference...
        let accessor = |state: &PolymerizerState<'a, 'p>, target, residue| {
            state.find_any_free_groups(&[target], residue, 1)
        };
        self.bond_with_accessors(abbr, donor, accessor, acceptor, accessor)
    }

    pub fn bond_groups(
        &mut self,
        abbr: impl AsRef<str>,
        donor: ResidueId,
        donor_group: &FunctionalGroup<'p>,
        acceptor: ResidueId,
        acceptor_group: &FunctionalGroup<'p>,
    ) -> Result<BondId> {
        // NOTE: Perhaps this explicit type annotation for `state` won't be necessary after the compiler improves its
        // type inference...
        let accessor = |group| {
            move |state: &PolymerizerState<'a, 'p>, target, residue| {
                state.find_these_free_groups(&[target], residue, slice::from_ref(group))
            }
        };
        self.bond_with_accessors(
            abbr,
            donor,
            accessor(donor_group),
            acceptor,
            accessor(acceptor_group),
        )
    }
}

// Private Methods =====================================================================================================

impl<'a, 'p> Polymer<'a, 'p> {
    fn bond_with_accessors<A, I>(
        &mut self,
        abbr: impl AsRef<str>,
        donor: ResidueId,
        donor_accessor: A,
        acceptor: ResidueId,
        acceptor_accessor: A,
    ) -> Result<BondId>
    where
        A: Fn(&PolymerizerState<'a, 'p>, &'p Target, ResidueId) -> Result<I, FindFreeGroupsError>,
        I: Iterator<Item = FunctionalGroup<'p>>,
    {
        let (
            abbr,
            BondDescription {
                name,
                lost,
                gained,
                from,
                to,
            },
        ) = Bond::lookup_description(self.polymer_db(), abbr)?;

        let residue_error = |id| Err(PolychemError::residue_not_in_polymer(id).into());
        let Some(donor_residue) = self.residue(donor) else {
            return residue_error(donor);
        };
        let Some(acceptor_residue) = self.residue(acceptor) else {
            return residue_error(acceptor);
        };

        let find_free_group = |accessor: A, target, residue, donor_or_acceptor| {
            accessor(&self.polymerizer_state, target, residue)
                // SAFETY: All accessors yield at least one group if they return `Ok`, so this unwrap shouldn't fail!
                .map(|mut gs| gs.next().unwrap())
                .map_err(|e| {
                    PolychemError::bond(
                        abbr,
                        donor,
                        donor_residue.name(),
                        acceptor,
                        acceptor_residue.name(),
                        donor_or_acceptor,
                        e,
                    )
                })
        };

        // Avoid partial updates by performing validation of both group updates *before* updating either group
        let donor_group = find_free_group(donor_accessor, from, donor, "donor")?;
        let acceptor_group = find_free_group(acceptor_accessor, to, acceptor, "acceptor")?;

        let bond = Bond {
            abbr,
            name,
            lost,
            gained,
        };
        let id = BondId(self.polymerizer_state.next_id());
        let bond_info = BondInfo(donor, bond, acceptor);
        self.bonds.insert(id, bond_info);

        self.update_group(donor, &donor_group, GroupState::Donor(id));
        self.update_group(acceptor, &acceptor_group, GroupState::Acceptor(id));

        Ok(id)
    }

    fn update_group(
        &mut self,
        residue_id: ResidueId,
        group: &FunctionalGroup<'p>,
        state: GroupState,
    ) {
        // SAFETY: `update_group()` is only called if the `residue_id` was found in the `group_index` of
        // `self.polymerizer_state`, and if a residue is present in the `group_index`, it should be present in the
        // `Polymer` as well!
        let residue = self.residues.get_mut(&residue_id).unwrap();
        let target = Target::from_residue_and_group(residue, group);

        // SAFETY: As long as this particular residue-group combination has been vetted by `find_*_free_groups()`, this
        // `unwrap()` shouldn't panic!
        self.polymerizer_state
            .update_group_index(target, residue_id, state)
            .unwrap();

        // SAFETY: Same as above, getting the `group` from `find_*_free_groups()` means this shouldn't ever panic!
        let group_state = residue.group_state_mut(group).unwrap();
        *group_state = state;
    }
}

// NOTE: Since we'll be calling the `$accessor` method on different structs (`Residue`s, `Modification`s, and `Bond`s),
// `$accessor` can't have a concrete type. You could solve this with trait objects, but that sacrifices performance —
// using a macro here keeps things DRY and overhead-free!
macro_rules! sum_parts {
    ($self:expr, $accessor:path) => {{
        let residues = $self.residues.values().map($accessor);
        let modifications = $self.modifications.values().map(|mod_info| match mod_info {
            ModificationInfo::Named(m, ..) => $accessor(m),
            ModificationInfo::Offset(m, ..) => $accessor(m),
            ModificationInfo::Unlocalized(m) => $accessor(m),
        });
        let bonds = $self.bonds.values().map(|BondInfo(_, b, _)| $accessor(b));

        residues.chain(modifications).chain(bonds).sum()
    }};
}

impl Massive for Polymer<'_, '_> {
    fn monoisotopic_mass(&self) -> MonoisotopicMass {
        sum_parts!(self, Massive::monoisotopic_mass)
    }

    fn average_mass(&self) -> AverageMass {
        sum_parts!(self, Massive::average_mass)
    }
}

impl Charged for Polymer<'_, '_> {
    fn charge(&self) -> Charge {
        sum_parts!(self, Charged::charge)
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_ron_snapshot;
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::{
        polymers::polymerizer::Polymerizer, testing_tools::assert_miette_snapshot, ChargedParticle,
        FunctionalGroup, ModificationId, ResidueId,
    };

    use super::*;

    const STEM_RESIDUES: [&str; 4] = ["A", "E", "J", "A"];

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);
    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "test_polymer_database.kdl",
            include_str!("../../tests/data/polymer_database.kdl"),
        )
        .unwrap()
    });

    static POLYMERIZER: Lazy<Polymerizer> = Lazy::new(|| Polymerizer::new(&ATOMIC_DB, &POLYMER_DB));

    macro_rules! assert_polymer {
        ($polymer:expr, $mono_mass:literal, $avg_mass:literal) => {
            assert_eq!(
                $polymer.monoisotopic_mass(),
                MonoisotopicMass(dec!($mono_mass))
            );
            assert_eq!($polymer.average_mass(), AverageMass(dec!($avg_mass)));
            assert_eq!($polymer.charge(), Charge(0));
            assert_eq!($polymer.monoisotopic_mz(), None);
            assert_eq!($polymer.average_mz(), None);
        };

        ($polymer:expr, $mono_mass:literal, $avg_mass:literal, $charge:literal, $mono_mz:literal, $avg_mz:literal) => {
            assert_eq!(
                $polymer.monoisotopic_mass(),
                MonoisotopicMass(dec!($mono_mass))
            );
            assert_eq!($polymer.average_mass(), AverageMass(dec!($avg_mass)));
            assert_eq!($polymer.charge(), Charge($charge));
            assert_eq!(
                $polymer.monoisotopic_mz(),
                Some(MonoisotopicMz(dec!($mono_mz)))
            );
            assert_eq!($polymer.average_mz(), Some(AverageMz(dec!($avg_mz))));
        };
    }

    #[test]
    fn recover_databases() {
        let polymer = POLYMERIZER.new_polymer();
        assert_eq!(polymer.atomic_db(), &*ATOMIC_DB);
        assert_eq!(polymer.polymer_db(), &*POLYMER_DB);
    }

    #[test]
    fn update_group() {
        let mut polymer = POLYMERIZER.new_polymer();
        let first_alanine = polymer.new_residue("A").unwrap();
        let second_alanine = polymer.new_residue("A").unwrap();
        let lysine = polymer.new_residue("K").unwrap();

        macro_rules! assert_residue_and_index_state {
            ($group:expr, $residues_and_states:expr) => {
                for (residue_id, state) in $residues_and_states {
                    let residue = polymer.residue(residue_id).unwrap();
                    let target = Target::from_residue_and_group(residue, $group);
                    let residue_state = *residue.group_state($group).unwrap();
                    let index_state = polymer
                        .polymerizer_state
                        .residue_groups(&[target], residue_id)
                        .next()
                        .unwrap()
                        .1;

                    assert_eq!(residue_state, state);
                    assert_eq!(index_state, state);
                }
            };
        }

        let n_terminal = FunctionalGroup::new("Amino", "N-Terminal");
        let n_terminal_groups = [
            (first_alanine, GroupState::Free),
            (second_alanine, GroupState::Free),
            (lysine, GroupState::Free),
        ];
        assert_residue_and_index_state!(&n_terminal, n_terminal_groups);

        polymer.update_group(
            first_alanine,
            &n_terminal,
            GroupState::Modified(ModificationId(42)),
        );
        let n_terminal_groups = [
            (first_alanine, GroupState::Modified(ModificationId(42))),
            (second_alanine, GroupState::Free),
            (lysine, GroupState::Free),
        ];
        assert_residue_and_index_state!(&n_terminal, n_terminal_groups);

        polymer.update_group(lysine, &n_terminal, GroupState::Acceptor(BondId(314)));
        let n_terminal_groups = [
            (first_alanine, GroupState::Modified(ModificationId(42))),
            (second_alanine, GroupState::Free),
            (lysine, GroupState::Acceptor(BondId(314))),
        ];
        assert_residue_and_index_state!(&n_terminal, n_terminal_groups);

        polymer.update_group(first_alanine, &n_terminal, GroupState::Donor(BondId(42)));
        let n_terminal_groups = [
            (first_alanine, GroupState::Donor(BondId(42))),
            (second_alanine, GroupState::Free),
            (lysine, GroupState::Acceptor(BondId(314))),
        ];
        assert_residue_and_index_state!(&n_terminal, n_terminal_groups);

        let c_terminal = FunctionalGroup::new("Carboxyl", "C-Terminal");
        polymer.update_group(second_alanine, &c_terminal, GroupState::Donor(BondId(1337)));
        let c_terminal_groups = [
            (ResidueId(0), GroupState::Free),
            (ResidueId(1), GroupState::Donor(BondId(1337))),
            (ResidueId(2), GroupState::Free),
        ];
        assert_residue_and_index_state!(&c_terminal, c_terminal_groups);
    }

    #[test]
    fn new_residue() {
        let mut polymer = POLYMERIZER.new_polymer();
        let residues = STEM_RESIDUES.map(|abbr| polymer.new_residue(abbr).unwrap());
        assert_eq!(
            residues,
            [ResidueId(0), ResidueId(1), ResidueId(2), ResidueId(3)]
        );
        assert_polymer!(polymer, 515.24387164950, 515.51342919034656875);

        let more_residues = STEM_RESIDUES.map(|abbr| polymer.new_residue(abbr).unwrap());
        assert_eq!(
            more_residues,
            [ResidueId(4), ResidueId(5), ResidueId(6), ResidueId(7)]
        );
        assert_polymer!(polymer, 1030.48774329900, 1031.02685838069313750);

        let residue_refs = residues.map(|id| polymer.residue(id).unwrap());
        assert_ron_snapshot!(residue_refs, {
            ".**.isotopes, .**.functional_groups" => insta::sorted_redaction()
        });

        let nonexistent_residue = polymer.new_residue("?");
        assert_miette_snapshot!(nonexistent_residue);

        let missing_residue = polymer.residue(ResidueId(8));
        assert_eq!(missing_residue, None);
    }

    #[test]
    fn bond_residues() {
        let mut polymer = POLYMERIZER.new_polymer();
        let murnac = polymer.new_residue("m").unwrap();
        let alanine = polymer.new_residue("A").unwrap();
        assert_polymer!(polymer, 382.15874504254, 382.36421782581807210);

        polymer.bond_residues("Stem", murnac, alanine).unwrap();
        assert_polymer!(polymer, 364.14818035851, 364.34893139338823950);

        let lactyl = FunctionalGroup::new("Carboxyl", "Lactyl Ether");
        let n_terminal = FunctionalGroup::new("Amino", "N-Terminal");
        assert!(polymer
            .residue(murnac)
            .unwrap()
            .group_state(&lactyl)
            .unwrap()
            .is_donor());
        assert!(polymer
            .residue(alanine)
            .unwrap()
            .group_state(&n_terminal)
            .unwrap()
            .is_acceptor());

        let all_groups_occupied = polymer.bond_residues("Stem", murnac, alanine);
        assert_miette_snapshot!(all_groups_occupied);

        // Try looking for residues that aren't present in this polymer
        let donor_not_in_polymer = polymer.bond_residues("Stem", ResidueId(2), alanine);
        assert_miette_snapshot!(donor_not_in_polymer);
        let acceptor_not_in_polymer = polymer.bond_residues("Stem", murnac, ResidueId(3));
        assert_miette_snapshot!(acceptor_not_in_polymer);

        let glcnac = polymer.new_residue("g").unwrap();
        let no_matching_groups = polymer.bond_residues("Pep", alanine, glcnac);
        assert_miette_snapshot!(no_matching_groups);
        // When bonding fails due to the acceptor, make sure that the donor remains untouched
        let c_terminal = FunctionalGroup::new("Carboxyl", "C-Terminal");
        assert!(polymer
            .residue(alanine)
            .unwrap()
            .group_state(&c_terminal)
            .unwrap()
            .is_free());

        let nonexistent_bond = polymer.bond_residues("Super", murnac, alanine);
        assert_miette_snapshot!(nonexistent_bond);
    }

    #[test]
    fn bond_groups() {
        let mut polymer = POLYMERIZER.new_polymer();
        let murnac = polymer.new_residue("m").unwrap();
        let alanine = polymer.new_residue("A").unwrap();
        assert_polymer!(polymer, 382.15874504254, 382.36421782581807210);

        let lactyl = FunctionalGroup::new("Carboxyl", "Lactyl Ether");
        let n_terminal = FunctionalGroup::new("Amino", "N-Terminal");
        polymer
            .bond_groups("Stem", murnac, &lactyl, alanine, &n_terminal)
            .unwrap();
        assert_polymer!(polymer, 364.14818035851, 364.34893139338823950);
        assert!(polymer
            .residue(murnac)
            .unwrap()
            .group_state(&lactyl)
            .unwrap()
            .is_donor());
        assert!(polymer
            .residue(alanine)
            .unwrap()
            .group_state(&n_terminal)
            .unwrap()
            .is_acceptor());

        // Use *`bond_residues()`* to check that an A-K link is indeed ambiguous
        let lysine = polymer.new_residue("K").unwrap();
        let ambigious_bond = polymer.bond_residues("Link", alanine, lysine);
        assert_miette_snapshot!(ambigious_bond);

        // But that ambiguity can be resolved using `bond_groups()`:
        let c_terminal = FunctionalGroup::new("Carboxyl", "C-Terminal");
        let sidechain = FunctionalGroup::new("Amino", "Sidechain");
        polymer
            .bond_groups("Link", alanine, &c_terminal, lysine, &sidechain)
            .unwrap();
        assert_polymer!(polymer, 492.24314337370, 492.52144716967893560);
        assert!(polymer
            .residue(alanine)
            .unwrap()
            .group_state(&c_terminal)
            .unwrap()
            .is_donor());
        assert!(polymer
            .residue(lysine)
            .unwrap()
            .group_state(&sidechain)
            .unwrap()
            .is_acceptor());
        // But the n_terminal of lysine remains untouched:
        assert!(polymer
            .residue(lysine)
            .unwrap()
            .group_state(&n_terminal)
            .unwrap()
            .is_free());

        let groups_not_free = polymer.bond_groups("Stem", murnac, &lactyl, alanine, &n_terminal);
        assert_miette_snapshot!(groups_not_free);

        let glcnac = polymer.new_residue("g").unwrap();
        let invalid_bond = polymer.bond_groups("Pep", lysine, &c_terminal, glcnac, &n_terminal);
        assert_miette_snapshot!(invalid_bond);
        // When bonding fails due to the acceptor, make sure that the donor remains untouched
        assert!(polymer
            .residue(lysine)
            .unwrap()
            .group_state(&c_terminal)
            .unwrap()
            .is_free());

        let nonexistent_bond = polymer.bond_groups("Super", murnac, &lactyl, alanine, &n_terminal);
        assert_miette_snapshot!(nonexistent_bond);
    }
}
