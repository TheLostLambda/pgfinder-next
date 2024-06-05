use core::slice;

use itertools::Itertools;

use crate::{
    errors::PolychemError,
    moieties::{
        polymer_database::{BondDescription, ModificationDescription},
        target::Target,
    },
    AnyModification, AtomicDatabase, AverageMass, Bond, BondId, BondInfo, Charge, Charged,
    FunctionalGroup, GroupState, Massive, ModificationId, ModificationInfo, MonoisotopicMass,
    NamedMod, Polymer, PolymerDatabase, Residue, ResidueGroup, ResidueId, Result,
};

use super::{errors::FindFreeGroupsError, polymerizer_state::PolymerizerState};

// FIXME: A horrible hack that's needed to specify the lifetimes captured by `-> impl Trait` correctly. Once Rust 2024
// is stabilized, however, this hack can be removed. Keep an eye on: https://github.com/rust-lang/rust/issues/117587
pub trait Captures<U> {}
impl<T: ?Sized, U> Captures<U> for T {}

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

    pub fn remove_residue(&mut self, id: ResidueId) -> Option<Residue<'a, 'p>> {
        let residue = self.residues.remove(&id);

        if let Some(residue) = &residue {
            self.polymerizer_state.unindex_residue_groups(id, residue);

            for (_, state) in residue.functional_groups() {
                match state {
                    GroupState::Free => (),
                    GroupState::Modified(id) => {
                        self.modifications.remove(&id);
                    }
                    GroupState::Donor(id) => {
                        if let Some(BondInfo(.., acceptor)) = self.bonds.remove(&id) {
                            self.update_group(&acceptor, GroupState::Free);
                        }
                    }
                    GroupState::Acceptor(id) => {
                        if let Some(BondInfo(donor, ..)) = self.bonds.remove(&id) {
                            self.update_group(&donor, GroupState::Free);
                        }
                    }
                }
            }
        }

        residue
    }

    #[must_use]
    pub fn residue(&self, id: ResidueId) -> Option<&Residue<'a, 'p>> {
        self.residues.get(&id)
    }

    pub fn new_chain(
        &mut self,
        abbr: impl AsRef<str>,
        residues: &[impl AsRef<str>],
    ) -> Result<(
        impl Iterator<Item = ResidueId>,
        impl Iterator<Item = BondId>,
    )> {
        let residue_ids: Vec<_> = residues
            .iter()
            .map(|abbr| self.new_residue(abbr))
            .try_collect()?;
        let bond_ids = self.bond_chain(abbr, &residue_ids)?;
        Ok((residue_ids.into_iter(), bond_ids))
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

    pub fn bond_chain(
        &mut self,
        abbr: impl AsRef<str>,
        residues: &[ResidueId],
    ) -> Result<impl Iterator<Item = BondId>> {
        let abbr = abbr.as_ref();
        // NOTE: Once `array_windows()` is stabilised, it should be possible to pattern-match the `residue_pair` and
        // eliminate the need for any manual indexing. Keep an eye on the standard library!
        let bond_ids: Vec<_> = residues
            .windows(2)
            // SAFETY: The indexing of `residue_pair` should never fail, since `windows(2)` will always return a slice
            // containing at least two elements
            .map(|residue_pair| self.bond_residues(abbr, residue_pair[0], residue_pair[1]))
            .try_collect()?;

        Ok(bond_ids.into_iter())
    }

    pub fn modify(
        &mut self,
        // FIXME: Test that this can take NamedMod, OffsetMod, AnyMod, and Modification<_> with any of the preceeding —
        // that's six (6) test-cases in all!
        _modification: impl Into<AnyModification<'a, 'p>>,
        _target: ResidueId,
    ) -> Result<()> {
        // let Modification { multiplier, kind } = modification.into();
        // match kind {
        //     AnyMod::Named(kind) => {
        //         self.modify_with_optional_groups(kind.abbr(), target, multiplier)
        //     }
        //     AnyMod::Offset(kind) => {
        //         // FIXME: Once `Polymerizer` is refactored to store residues, then turn this into a method on `self`
        //         // FIXME: And when you do that, get rid of this nasty discard hack...
        //         target.add_offsets(kind, multiplier).map(|_| ())
        //     }
        // }
        todo!()
    }

    pub fn modify_only_group(
        &mut self,
        abbr: impl AsRef<str>,
        target: ResidueId,
    ) -> Result<ModificationId> {
        self.modify_only_groups(abbr, target, 1)
            // SAFETY: The `self.modify_only_groups()` method should yield at least one group if it returns `Ok`, so
            // this unwrap shouldn't fail!
            .map(|mut ids| ids.next().unwrap())
    }

    pub fn modify_only_groups(
        &mut self,
        abbr: impl AsRef<str>,
        target: ResidueId,
        number: usize,
    ) -> Result<impl Iterator<Item = ModificationId> + Captures<(&'a (), &'p ())>> {
        let accessor = move |state: &PolymerizerState<'a, 'p>, targets, residue| {
            state.find_any_free_groups(targets, residue, number)
        };
        self.modify_with_accessor(abbr, target, accessor)
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
                        name,
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

        let donor = ResidueGroup(donor, donor_group);
        let acceptor = ResidueGroup(acceptor, acceptor_group);
        let id = BondId(self.polymerizer_state.next_id());

        self.update_group(&donor, GroupState::Donor(id));
        self.update_group(&acceptor, GroupState::Acceptor(id));

        let bond = Bond {
            abbr,
            name,
            lost,
            gained,
        };
        let bond_info = BondInfo(donor, bond, acceptor);
        self.bonds.insert(id, bond_info);

        Ok(id)
    }

    fn modify_with_accessor<A, I>(
        &mut self,
        abbr: impl AsRef<str>,
        residue: ResidueId,
        accessor: A,
    ) -> Result<impl Iterator<Item = ModificationId>>
    where
        A: Fn(&PolymerizerState<'a, 'p>, &'p [Target], ResidueId) -> Result<I, FindFreeGroupsError>,
        I: Iterator<Item = FunctionalGroup<'p>>,
    {
        let (
            abbr,
            ModificationDescription {
                name,
                lost,
                gained,
                targets,
            },
        ) = NamedMod::lookup_description(self.polymer_db(), abbr)?;

        let Some(residue_ref) = self.residue(residue) else {
            return Err(PolychemError::residue_not_in_polymer(residue).into());
        };

        let target_groups = accessor(&self.polymerizer_state, targets, residue)
            .map_err(|e| PolychemError::named_modification(name, residue, residue_ref.name(), e))?;

        // NOTE: The reason this is `collect()`ed just to be turned back into an iterator is to ensure that this lazy
        // `map()` code runs before the function returns. Otherwise the requested modifications would only be added to
        // the `Polymer` after the caller has consumed the returned iterator of `ModificationId`s.
        let modification_ids: Vec<_> = target_groups
            .map(|group| {
                let residue_group = ResidueGroup(residue, group);
                let id = ModificationId(self.polymerizer_state.next_id());
                self.update_group(&residue_group, GroupState::Modified(id));

                let modification = NamedMod {
                    abbr,
                    name,
                    lost,
                    gained,
                };
                let modification_info = ModificationInfo::Named(modification, residue_group);
                self.modifications.insert(id, modification_info);

                id
            })
            .collect();

        Ok(modification_ids.into_iter())
    }

    fn update_group(&mut self, residue_group: &ResidueGroup<'p>, state: GroupState) {
        let &ResidueGroup(residue_id, ref group) = residue_group;
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
        polymers::polymerizer::Polymerizer, testing_tools::assert_miette_snapshot, AverageMz,
        ChargedParticle, FunctionalGroup, ModificationId, MonoisotopicMz, ResidueId,
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
            &ResidueGroup(first_alanine, n_terminal),
            GroupState::Modified(ModificationId(42)),
        );
        let n_terminal_groups = [
            (first_alanine, GroupState::Modified(ModificationId(42))),
            (second_alanine, GroupState::Free),
            (lysine, GroupState::Free),
        ];
        assert_residue_and_index_state!(&n_terminal, n_terminal_groups);

        polymer.update_group(
            &ResidueGroup(lysine, n_terminal),
            GroupState::Acceptor(BondId(314)),
        );
        let n_terminal_groups = [
            (first_alanine, GroupState::Modified(ModificationId(42))),
            (second_alanine, GroupState::Free),
            (lysine, GroupState::Acceptor(BondId(314))),
        ];
        assert_residue_and_index_state!(&n_terminal, n_terminal_groups);

        polymer.update_group(
            &ResidueGroup(first_alanine, n_terminal),
            GroupState::Donor(BondId(42)),
        );
        let n_terminal_groups = [
            (first_alanine, GroupState::Donor(BondId(42))),
            (second_alanine, GroupState::Free),
            (lysine, GroupState::Acceptor(BondId(314))),
        ];
        assert_residue_and_index_state!(&n_terminal, n_terminal_groups);

        let c_terminal = FunctionalGroup::new("Carboxyl", "C-Terminal");
        polymer.update_group(
            &ResidueGroup(second_alanine, c_terminal),
            GroupState::Donor(BondId(1337)),
        );
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
    fn remove_residue() {
        let mut polymer = POLYMERIZER.new_polymer();
        let murnac = polymer.new_residue("m").unwrap();
        let alanine = polymer.new_residue("A").unwrap();

        polymer.bond_residues("Stem", murnac, alanine).unwrap();
        assert_polymer!(polymer, 364.14818035851, 364.34893139338823950);

        polymer.modify_only_group("Red", murnac).unwrap();
        assert_polymer!(polymer, 366.16383042297, 366.36481290149979420);

        polymer.modify_only_group("Ca", alanine).unwrap();
        assert_polymer!(
            polymer,
            405.118047659530870,
            405.43446178607839420,
            1,
            405.118047659530870,
            405.43446178607839420
        );

        // Remove the MurNAc residue (and its modifications / bonds)
        polymer.remove_residue(murnac);
        assert_polymer!(
            polymer,
            128.001895705740870,
            128.16295491325714225,
            1,
            128.001895705740870,
            128.16295491325714225
        );

        // Clear the rest of the polymer
        polymer.remove_residue(alanine);
        assert_polymer!(polymer, 0.0, 0.0);
    }

    #[test]
    fn new_chain() {
        // NOTE: A macro since using a closure leads to borrow-checker issues...
        macro_rules! bond_chain {
            ($abbr:literal, $residues:expr) => {{
                let mut polymer = POLYMERIZER.new_polymer();
                polymer
                    .new_chain($abbr, $residues)
                    .map(|(residue_ids, bond_ids)| {
                        (
                            polymer,
                            residue_ids.collect::<Vec<_>>(),
                            bond_ids.collect::<Vec<_>>(),
                        )
                    })
            }};
        }

        // Building chains of 1 or 0 residues results in no bonds
        let empty: &[&str] = &[];
        let (polymer, residues, bonds) = bond_chain!("Pep", empty).unwrap();
        assert!(residues.is_empty());
        assert!(bonds.is_empty());
        assert_polymer!(polymer, 0.0, 0.0);

        let (polymer, residues, bonds) = bond_chain!("Pep", &["A"]).unwrap();
        assert_eq!(residues, vec![ResidueId(0)]);
        assert!(bonds.is_empty());
        assert_polymer!(polymer, 89.04767846918, 89.09330602867854225);

        // Then actually build a full chain (3 bonds for 4 residues)
        let (polymer, residues, bonds) = bond_chain!("Pep", &STEM_RESIDUES).unwrap();
        assert_eq!(
            residues,
            vec![ResidueId(0), ResidueId(1), ResidueId(2), ResidueId(3)]
        );
        assert_eq!(bonds, vec![BondId(4), BondId(5), BondId(6)]);
        assert_polymer!(polymer, 461.21217759741, 461.46756989305707095);
        assert_ron_snapshot!(polymer, {
            ".**.composition, .**.lost, .**.gained" => "<FORMULA>",
            ".**.residues, .**.bonds" => insta::sorted_redaction(),
            ".**.functional_groups" => insta::sorted_redaction()
        });

        let nonexistent_bond = bond_chain!("?", &STEM_RESIDUES);
        assert_miette_snapshot!(nonexistent_bond);

        let nonexistent_residue = bond_chain!("Pep", &["A", "yE", "J"]);
        assert_miette_snapshot!(nonexistent_residue);
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

    #[test]
    fn bond_chain() {
        let mut polymer = POLYMERIZER.new_polymer();
        let residues = STEM_RESIDUES.map(|abbr| polymer.new_residue(abbr).unwrap());
        let polymer_without_bonds = polymer.clone();
        assert_polymer!(polymer, 515.24387164950, 515.51342919034656875);

        // NOTE: A macro since using a closure leads to borrow-checker issues...
        macro_rules! bond_chain {
            ($abbr:literal, $residues:expr) => {
                polymer
                    .bond_chain($abbr, $residues)
                    .map(|ids| ids.collect::<Vec<_>>())
            };
        }

        // Bonding chains of 1 or 0 residues is a no-op
        let new_bonds = bond_chain!("Pep", &[]).unwrap();
        assert!(new_bonds.is_empty());
        assert_eq!(polymer, polymer_without_bonds);
        let new_bonds = bond_chain!("Pep", &residues[..1]).unwrap();
        assert!(new_bonds.is_empty());
        assert_eq!(polymer, polymer_without_bonds);

        // Then actually bond the full chain (3 bonds for 4 residues)
        let bond_ids: Vec<_> = bond_chain!("Pep", &residues).unwrap();
        assert_eq!(bond_ids, vec![BondId(4), BondId(5), BondId(6)]);
        assert_polymer!(polymer, 461.21217759741, 461.46756989305707095);
        assert_ron_snapshot!(polymer, {
            ".**.composition, .**.lost, .**.gained" => "<FORMULA>",
            ".**.residues, .**.bonds" => insta::sorted_redaction(),
            ".**.functional_groups" => insta::sorted_redaction()
        });

        let nonexistent_bond = bond_chain!("?", &residues);
        assert_miette_snapshot!(nonexistent_bond);

        let nonexistent_residue = bond_chain!("Pep", &[ResidueId(3), ResidueId(4)]);
        assert_miette_snapshot!(nonexistent_residue);
    }

    #[ignore]
    #[test]
    fn modify() {
        // TODO: Restore from git
        todo!()
    }

    #[test]
    fn modify_only_group() {
        let mut polymer = POLYMERIZER.new_polymer();
        let murnac = polymer.new_residue("m").unwrap();
        assert_polymer!(polymer, 293.11106657336, 293.27091179713952985);

        let modification_id = polymer.modify_only_group("Anh", murnac).unwrap();
        assert_eq!(modification_id, ModificationId(1));
        assert_polymer!(polymer, 275.10050188933, 275.25562536470969725);
        let reducing_end = FunctionalGroup::new("Hydroxyl", "Reducing End");
        assert!(polymer
            .residue(murnac)
            .unwrap()
            .group_state(&reducing_end)
            .unwrap()
            .is_modified());

        let all_groups_occupied = polymer.modify_only_group("Anh", murnac);
        assert_miette_snapshot!(all_groups_occupied);

        let residue_not_in_polymer = polymer.modify_only_group("Anh", ResidueId(1));
        assert_miette_snapshot!(residue_not_in_polymer);

        let no_matching_groups = polymer.modify_only_group("Am", murnac);
        assert_miette_snapshot!(no_matching_groups);

        let nonexistent_modification = polymer.modify_only_group("Arg", murnac);
        assert_miette_snapshot!(nonexistent_modification);
    }

    #[test]
    fn modify_only_groups() {
        let mut polymer = POLYMERIZER.new_polymer();
        let murnac = polymer.new_residue("m").unwrap();
        assert_polymer!(polymer, 293.11106657336, 293.27091179713952985);

        // NOTE: A macro since using a closure leads to borrow-checker issues...
        macro_rules! modify_only_groups {
            ($abbr:literal, $target:expr, $number:expr) => {
                polymer
                    .modify_only_groups($abbr, $target, $number)
                    .map(Itertools::collect_vec)
            };
        }

        let modification_ids = modify_only_groups!("Met", murnac, 3).unwrap();
        assert_eq!(
            modification_ids,
            vec![ModificationId(1), ModificationId(2), ModificationId(3)]
        );
        assert_polymer!(polymer, 335.15801676674, 335.35076401167994095);
        let hydroxyl_groups = [
            FunctionalGroup::new("Hydroxyl", "Reducing End"),
            FunctionalGroup::new("Hydroxyl", "Nonreducing End"),
            FunctionalGroup::new("Hydroxyl", "6-Position"),
        ];
        for hydroxyl_group in hydroxyl_groups {
            assert!(polymer
                .residue(murnac)
                .unwrap()
                .group_state(&hydroxyl_group)
                .unwrap()
                .is_modified());
        }

        polymer.remove_residue(murnac);
        let murnac = polymer.new_residue("m").unwrap();
        assert_polymer!(polymer, 293.11106657336, 293.27091179713952985);

        let modification_ids = modify_only_groups!("Ca", murnac, 4).unwrap();
        assert_eq!(
            modification_ids,
            vec![
                ModificationId(5),
                ModificationId(6),
                ModificationId(7),
                ModificationId(8)
            ]
        );
        assert_polymer!(
            polymer,
            448.927935519603480,
            449.54950733545392985,
            4,
            112.231983879900870,
            112.3873768338634824625
        );

        polymer.remove_residue(murnac);
        let alanine = polymer.new_residue("A").unwrap();
        assert_polymer!(polymer, 89.04767846918, 89.09330602867854225);
        let modification_ids = modify_only_groups!("Ca", alanine, 2).unwrap();
        assert_eq!(
            modification_ids,
            vec![ModificationId(10), ModificationId(11)]
        );
        assert_polymer!(
            polymer,
            166.956112942301740,
            167.23260379783574225,
            2,
            83.478056471150870,
            83.6163018989178711250
        );

        polymer.remove_residue(alanine);
        let murnac = polymer.new_residue("m").unwrap();
        let modification_ids = modify_only_groups!("Met", murnac, 3).unwrap();
        assert_eq!(
            modification_ids,
            vec![ModificationId(13), ModificationId(14), ModificationId(15)]
        );
        let all_groups_occupied = modify_only_groups!("Ca", murnac, 4);
        assert_miette_snapshot!(all_groups_occupied);
        let still_all_groups_occupied = modify_only_groups!("Ca", murnac, 2);
        assert_miette_snapshot!(still_all_groups_occupied);

        assert_polymer!(polymer, 335.15801676674, 335.35076401167994095);
        let modification_ids = modify_only_groups!("Ca", murnac, 1).unwrap();
        assert_eq!(modification_ids, vec![ModificationId(16)]);
        assert_polymer!(
            polymer,
            374.112234003300870,
            374.42041289625854095,
            1,
            374.112234003300870,
            374.42041289625854095
        );

        let residue_not_in_polymer = modify_only_groups!("Ca", ResidueId(42), 2);
        assert_miette_snapshot!(residue_not_in_polymer);

        let murnac = polymer.new_residue("m").unwrap();
        let no_matching_groups = modify_only_groups!("Anh", murnac, 2);
        assert_miette_snapshot!(no_matching_groups);
    }
}
