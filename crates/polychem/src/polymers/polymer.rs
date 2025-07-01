use core::slice;
use std::borrow::Borrow;

use ahash::HashSet;
use itertools::Itertools;

use crate::{
    AnyMod, AnyModification, AtomicDatabase, AverageMass, Bond, BondId, BondInfo, Charge, Charged,
    ChemicalComposition, Count, FunctionalGroup, GroupState, Massive, Modification, ModificationId,
    ModificationInfo, MonoisotopicMass, NamedMod, OffsetKind, OffsetMod, Polymer, PolymerDatabase,
    Residue, ResidueGroup, ResidueId, Result,
    errors::PolychemError,
    moieties::{
        polymer_database::{BondDescription, ModificationDescription},
        target::Target,
    },
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

    pub fn remove_residue(&mut self, id: ResidueId) -> Option<Residue<'a, 'p>> {
        let residue = self.residues.remove(&id);

        if let Some(residue) = &residue {
            self.polymerizer_state.unindex_residue_groups(id, residue);

            for (_, state) in residue.functional_groups() {
                match state {
                    GroupState::Free => (),
                    GroupState::Modified(id) => {
                        self.modifications.remove(id);
                    }
                    GroupState::Donor(id) => {
                        if let Some(BondInfo(.., acceptor)) = self.bonds.remove(id) {
                            self.update_group(&acceptor, GroupState::Free);
                        }
                    }
                    GroupState::Acceptor(id) => {
                        if let Some(BondInfo(donor, ..)) = self.bonds.remove(id) {
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

    // FIXME: This is part of an incomplete set of methods — be sure to add in methods for iterating over pairs and
    // also values? Needs more design thought in general! Maybe I should *only* have the pair version?
    // FIXME: Also needs testing!
    pub fn residue_ids(&self) -> impl Iterator<Item = ResidueId> {
        self.residues.keys().copied()
    }

    // FIXME: This feels a bit inconsistent with the `named_mod` that I've been using elsewhere? Maybe instead of
    // `named_mod`s and `offset_mod`s, I should just go with `modification`s and `offset`s?
    pub fn new_modification(
        &mut self,
        multiplier: impl TryInto<Count>,
        abbr: impl AsRef<str>,
    ) -> Result<ModificationId> {
        let multiplier = multiplier
            .try_into()
            .map_err(|_| PolychemError::ZeroMultiplier)?;
        let modification = Modification::new(multiplier, NamedMod::new(self.polymer_db(), abbr)?);

        Ok(self.unlocalized_modification(modification))
    }

    // FIXME: Again a bit unsure about the naming here...
    pub fn new_offset(
        &mut self,
        kind: OffsetKind,
        multiplier: impl TryInto<Count>,
        formula: impl AsRef<str>,
    ) -> Result<ModificationId> {
        let modification = self.build_offset_mod(kind, multiplier, formula)?;
        Ok(self.unlocalized_modification(modification))
    }

    pub fn new_offset_with_composition(
        &mut self,
        kind: OffsetKind,
        multiplier: impl TryInto<Count>,
        composition: ChemicalComposition<'a>,
    ) -> Result<ModificationId> {
        let modification = Self::build_offset_mod_with_composition(kind, multiplier, composition)?;
        Ok(self.unlocalized_modification(modification))
    }

    // TODO: Add `remove_modification` and be sure to clear any groups the modification was attached to!

    #[must_use]
    pub fn modification(&self, id: ModificationId) -> Option<&ModificationInfo<'_, '_>> {
        self.modifications.get(&id)
    }

    // FIXME: This is part of an incomplete set of methods — be sure to add in methods for iterating over IDs and
    // also just values? Needs more design thought in general! For example, should I have a method that just returns
    // unlocalized modifications?
    // FIXME: Also needs testing!
    pub fn modifications(
        &self,
    ) -> impl Iterator<Item = (ModificationId, &ModificationInfo<'a, 'p>)> {
        self.modifications.iter().map(|(&id, info)| (id, info))
    }

    // FIXME: This is part of an incomplete set of methods — be sure to add in methods for iterating over IDs and
    // also just values? Needs more design thought in general! For example, should I have a method that just returns
    // unlocalized modifications?
    // FIXME: Also needs testing!
    pub fn modification_refs(&self) -> impl Iterator<Item = &ModificationInfo<'a, 'p>> {
        self.modifications.values()
    }

    pub fn new_chain(
        &mut self,
        abbr: impl AsRef<str>,
        residues: impl IntoIterator<Item: AsRef<str>>,
    ) -> Result<(Vec<ResidueId>, Vec<BondId>)> {
        let residue_ids: Vec<_> = residues
            .into_iter()
            .map(|abbr| self.new_residue(abbr))
            .try_collect()?;
        let bond_ids = self.bond_chain(abbr, &residue_ids)?;
        Ok((residue_ids, bond_ids))
    }

    pub fn bond_residues(
        &mut self,
        abbr: impl AsRef<str>,
        donor: ResidueId,
        acceptor: ResidueId,
    ) -> Result<BondId> {
        // NOTE: Perhaps this explicit type annotation for `state` won't be necessary after the compiler improves its
        // type inference...
        let accessor = |residue| {
            move |state: &PolymerizerState<'a, 'p>, target| {
                state.find_any_free_groups(&[target], residue, 1)
            }
        };
        self.bond_with_accessors(abbr, donor, accessor(donor), acceptor, accessor(acceptor))
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
        let accessor = |residue, group| {
            move |state: &PolymerizerState<'a, 'p>, target| {
                state.find_these_free_groups(&[target], residue, slice::from_ref(group))
            }
        };
        self.bond_with_accessors(
            abbr,
            donor,
            accessor(donor, donor_group),
            acceptor,
            accessor(acceptor, acceptor_group),
        )
    }

    pub fn bond_chain(
        &mut self,
        abbr: impl AsRef<str>,
        residues: impl IntoIterator<Item: Borrow<ResidueId>>,
    ) -> Result<Vec<BondId>> {
        let abbr = abbr.as_ref();
        let bond_ids = residues
            .into_iter()
            .map(|r| r.borrow().to_owned())
            .tuple_windows()
            .map(|(donor, acceptor)| self.bond_residues(abbr, donor, acceptor))
            .try_collect()?;

        Ok(bond_ids)
    }

    // FIXME: Also needs testing!
    pub fn remove_bond(&mut self, id: BondId) -> Option<BondInfo<'a, 'p>> {
        let bond = self.bonds.remove(&id);

        if let Some(BondInfo(donor, _, acceptor)) = &bond {
            self.update_group(donor, GroupState::Free);
            self.update_group(acceptor, GroupState::Free);
        }

        bond
    }
    // TODO: Need to add `bond()` accessor / lookup function like I have for `residue()` and `modification()`

    // FIXME: This is part of an incomplete set of methods — be sure to add in methods for iterating over IDs and
    // also pairs? Needs more design thought in general!
    // FIXME: Also needs testing!
    pub fn bond_refs(&self) -> impl Iterator<Item = &BondInfo<'a, 'p>> {
        self.bonds.values()
    }

    // FIXME: This is part of an incomplete set of methods — be sure to add in methods for iterating over IDs and
    // also pairs? Needs more design thought in general!
    // FIXME: Also needs testing!
    pub fn bonds(&self) -> impl Iterator<Item = (BondId, &BondInfo<'a, 'p>)> {
        self.bonds.iter().map(|(&id, info)| (id, info))
    }

    pub fn localize_modification(
        &mut self,
        modification: ModificationId,
        residue: ResidueId,
    ) -> Result<Vec<ModificationId>> {
        // NOTE: We're `remove()`-ing the modification here with the presumption that the user is trying to localize a
        // currently unlocalized modification. In the error case where this is not true, the modification is reinserted
        // into the map before returning an error. The avoids an extra lookup in the "happy path"
        let Modification { multiplier, kind } = match self.modifications.remove(&modification) {
            Some(ModificationInfo::Unlocalized(m)) => Ok(m),
            Some(info) => {
                let error = PolychemError::modification_already_localized(modification, &info);
                self.modifications.insert(modification, info);
                Err(error)
            }
            None => Err(PolychemError::modification_not_in_polymer(modification)),
        }?;

        match kind {
            // NOTE: We're only extracting the `abbr` since we'll need to perform another lookup to check for legal
            // targets anyways, and we may also need to create more than one modification if the original, unlocalized
            // modification had a multiplier greater than one!
            AnyMod::Named(kind) => self.modify_only_groups(kind.abbr(), residue, multiplier),
            AnyMod::Offset(kind) => {
                let modification = Modification::new(multiplier, kind);
                self.offset_with_modification(modification, residue)
                    .map(|id| vec![id])
            }
        }
    }

    // FIXME: Look into deduplicating with `modify_with_accessors()` — could use DRYing...
    pub fn modify_polymer(&mut self, abbr: impl AsRef<str>) -> Result<Vec<ModificationId>> {
        let (
            abbr,
            ModificationDescription {
                name,
                lost,
                gained,
                targets,
            },
        ) = NamedMod::lookup_description(self.polymer_db(), abbr)?;

        // NOTE: Without `.collect()`ing here, I'm holding on to a reference to `self.polymerizer_state` which prevents
        // me from calling `self.group_modification()` later, since that needs a unique mutable reference to the
        // `polymerizer_state` (to mark functional groups as modified)
        #[allow(clippy::needless_collect)]
        let residue_groups: Vec<_> = self
            .polymerizer_state
            .free_polymer_groups(targets)
            .collect();

        Ok(residue_groups
            .into_iter()
            .map(|residue_group| {
                let modification = NamedMod {
                    abbr,
                    name,
                    lost,
                    gained,
                };
                self.group_modification(modification, residue_group)
            })
            .collect())
    }

    pub fn modify_only_group(
        &mut self,
        abbr: impl AsRef<str>,
        residue: ResidueId,
    ) -> Result<ModificationId> {
        self.modify_only_groups(abbr, residue, 1)
            // SAFETY: The `self.modify_only_groups()` method should yield at least one group if it returns `Ok`, so
            // this indexing should never panic!
            .map(|ids| ids[0])
    }

    // FIXME: Think about changing the argument order so that `number` comes first? Is that more consistent with
    // `offset_residue()` and friends?
    pub fn modify_only_groups(
        &mut self,
        abbr: impl AsRef<str>,
        residue: ResidueId,
        number: impl TryInto<Count> + Copy,
    ) -> Result<Vec<ModificationId>> {
        let accessor = |state: &PolymerizerState<'a, 'p>, targets| {
            state.find_any_free_groups(targets, residue, number)
        };
        self.modify_with_accessor(abbr, residue, accessor)
    }

    pub fn modify_group(
        &mut self,
        abbr: impl AsRef<str>,
        residue: ResidueId,
        group: &FunctionalGroup<'p>,
    ) -> Result<ModificationId> {
        self.modify_groups(abbr, residue, slice::from_ref(group))
            // SAFETY: The `self.modify_groups()` method should yield at least one group if it returns `Ok`, so this
            // indexing should never panic!
            .map(|ids| ids[0])
    }

    pub fn modify_groups(
        &mut self,
        abbr: impl AsRef<str>,
        residue: ResidueId,
        groups: impl IntoIterator<Item: Borrow<FunctionalGroup<'p>>> + Copy,
    ) -> Result<Vec<ModificationId>> {
        let accessor = |state: &PolymerizerState<'a, 'p>, targets| {
            state.find_these_free_groups(targets, residue, groups)
        };
        self.modify_with_accessor(abbr, residue, accessor)
    }

    pub fn offset_residue(
        &mut self,
        kind: OffsetKind,
        multiplier: impl TryInto<Count>,
        formula: impl AsRef<str>,
        residue: ResidueId,
    ) -> Result<ModificationId> {
        let modification = self.build_offset_mod(kind, multiplier, formula)?;
        self.offset_with_modification(modification, residue)
    }

    pub fn offset_residue_with_composition(
        &mut self,
        kind: OffsetKind,
        multiplier: impl TryInto<Count>,
        composition: ChemicalComposition<'a>,
        residue: ResidueId,
    ) -> Result<ModificationId> {
        let modification = Self::build_offset_mod_with_composition(kind, multiplier, composition)?;
        self.offset_with_modification(modification, residue)
    }
}

// Private Methods =====================================================================================================

impl<'a, 'p> Polymer<'a, 'p> {
    fn build_offset_mod(
        &self,
        kind: OffsetKind,
        multiplier: impl TryInto<Count>,
        formula: impl AsRef<str>,
    ) -> Result<Modification<OffsetMod<'a>>> {
        let composition = ChemicalComposition::new(self.atomic_db(), formula)?;

        Self::build_offset_mod_with_composition(kind, multiplier, composition)
    }

    fn build_offset_mod_with_composition(
        kind: OffsetKind,
        multiplier: impl TryInto<Count>,
        composition: ChemicalComposition<'a>,
    ) -> Result<Modification<OffsetMod<'a>>> {
        let multiplier = multiplier
            .try_into()
            .map_err(|_| PolychemError::ZeroMultiplier)?;
        let modification = Modification::new(multiplier, OffsetMod::new(kind, composition));

        Ok(modification)
    }

    fn unlocalized_modification(
        &mut self,
        modification: impl Into<AnyModification<'a, 'p>>,
    ) -> ModificationId {
        let id = ModificationId(self.polymerizer_state.next_id());
        let modification_info = ModificationInfo::Unlocalized(modification.into());
        self.modifications.insert(id, modification_info);

        id
    }

    fn group_modification(
        &mut self,
        modification: NamedMod<'a, 'p>,
        residue_group: ResidueGroup<'p>,
    ) -> ModificationId {
        let id = ModificationId(self.polymerizer_state.next_id());
        self.update_group(&residue_group, GroupState::Modified(id));

        let modification_info = ModificationInfo::Named(modification, residue_group);
        self.modifications.insert(id, modification_info);

        id
    }

    fn bond_with_accessors<A>(
        &mut self,
        abbr: impl AsRef<str>,
        donor: ResidueId,
        donor_accessor: A,
        acceptor: ResidueId,
        acceptor_accessor: A,
    ) -> Result<BondId>
    where
        A: Fn(
            &PolymerizerState<'a, 'p>,
            &'p Target,
        ) -> Result<HashSet<FunctionalGroup<'p>>, FindFreeGroupsError>,
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

        let find_free_group = |accessor: A, target, donor_or_acceptor| {
            accessor(&self.polymerizer_state, target)
                // SAFETY: All accessors yield at least one group if they return `Ok`, so this unwrap shouldn't fail!
                .map(|gs| gs.into_iter().next().unwrap())
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
        let donor_group = find_free_group(donor_accessor, from, "donor")?;
        let acceptor_group = find_free_group(acceptor_accessor, to, "acceptor")?;

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

    fn modify_with_accessor<A>(
        &mut self,
        abbr: impl AsRef<str>,
        residue: ResidueId,
        accessor: A,
    ) -> Result<Vec<ModificationId>>
    where
        A: Fn(
            &PolymerizerState<'a, 'p>,
            &'p [Target],
        ) -> Result<HashSet<FunctionalGroup<'p>>, FindFreeGroupsError>,
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

        let target_groups = accessor(&self.polymerizer_state, targets)
            .map_err(|e| PolychemError::named_modification(name, residue, residue_ref.name(), e))?;

        let modification_ids = target_groups
            .into_iter()
            .map(|group| {
                let modification = NamedMod {
                    abbr,
                    name,
                    lost,
                    gained,
                };
                let residue_group = ResidueGroup(residue, group);

                self.group_modification(modification, residue_group)
            })
            .collect();

        Ok(modification_ids)
    }

    fn offset_with_modification(
        &mut self,
        modification: Modification<OffsetMod<'a>>,
        residue: ResidueId,
    ) -> Result<ModificationId> {
        let Some(residue_ref) = self.residues.get_mut(&residue) else {
            return Err(PolychemError::residue_not_in_polymer(residue).into());
        };

        let id = ModificationId(self.polymerizer_state.next_id());
        residue_ref.offset_modifications_mut().insert(id);

        let modification_info = ModificationInfo::Offset(modification, residue);
        self.modifications.insert(id, modification_info);

        Ok(id)
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
    use rust_decimal_macros::dec;
    use std::sync::LazyLock;

    use crate::{
        AverageMz, ChargedParticle, FunctionalGroup, ModificationId, MonoisotopicMz, OffsetKind,
        ResidueId, polymers::polymerizer::Polymerizer, testing_tools::assert_miette_snapshot,
    };

    use super::*;

    const STEM_RESIDUES: [&str; 4] = ["A", "E", "J", "A"];

    static ATOMIC_DB: LazyLock<AtomicDatabase> = LazyLock::new(AtomicDatabase::default);
    static POLYMER_DB: LazyLock<PolymerDatabase> = LazyLock::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "test_polymer_database.kdl",
            include_str!("../../tests/data/polymer_database.kdl"),
        )
        .unwrap()
    });

    static POLYMERIZER: LazyLock<Polymerizer> =
        LazyLock::new(|| Polymerizer::new(&ATOMIC_DB, &POLYMER_DB));

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
    #[allow(clippy::cognitive_complexity)]
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
    fn new_modification() {
        let mut polymer = POLYMERIZER.new_polymer();

        let amidation = polymer.new_modification(1, "Am").unwrap();
        assert_eq!(amidation, ModificationId(0));
        assert_polymer!(polymer, -0.98401558291, -0.98476095881670255);

        let double_amidation = polymer.new_modification(2, "Am").unwrap();
        assert_eq!(double_amidation, ModificationId(1));
        assert_polymer!(polymer, -2.95204674873, -2.95428287645010765);

        let calcium = polymer.new_modification(2, "Ca").unwrap();
        assert_eq!(calcium, ModificationId(2));
        assert_polymer!(
            polymer,
            74.956387724391740,
            75.18501489270709235,
            2,
            37.478193862195870,
            37.5925074463535461750
        );

        let modification_refs =
            [amidation, double_amidation, calcium].map(|id| polymer.modification(id).unwrap());
        assert_ron_snapshot!(modification_refs, {
            ".**.lost, .**.gained" => "<FORMULA>",
        });

        let zero_multiplier = polymer.new_modification(0, "Ac");
        assert_miette_snapshot!(zero_multiplier);

        let nonexistent_modification = polymer.new_modification(1, "?");
        assert_miette_snapshot!(nonexistent_modification);

        let missing_modification = polymer.modification(ModificationId(8));
        assert_eq!(missing_modification, None);
    }

    #[test]
    fn new_offset() {
        let mut polymer = POLYMERIZER.new_polymer();

        let water_loss = polymer.new_offset(OffsetKind::Remove, 1, "H2O").unwrap();
        assert_eq!(water_loss, ModificationId(0));
        assert_polymer!(polymer, -18.01056468403, -18.01528643242983260);

        let double_water_gain = polymer.new_offset(OffsetKind::Add, 2, "H2O").unwrap();
        assert_eq!(double_water_gain, ModificationId(1));
        assert_polymer!(polymer, 18.01056468403, 18.01528643242983260);

        let triple_ammonium_loss = polymer.new_offset(OffsetKind::Remove, 3, "NH3+p").unwrap();
        assert_eq!(triple_ammonium_loss, ModificationId(2));
        assert_polymer!(
            polymer,
            -36.090912019193,
            -36.09811938827255755,
            -3,
            -12.030304006397666666666666667,
            -12.032706462757519183333333333
        );

        let modification_refs = [water_loss, double_water_gain, triple_ammonium_loss]
            .map(|id| polymer.modification(id).unwrap());
        assert_ron_snapshot!(modification_refs, {
            ".**.composition" => "<FORMULA>",
        });

        let zero_multiplier = polymer.new_offset(OffsetKind::Remove, 0, "H");
        assert_miette_snapshot!(zero_multiplier);

        let invalid_composition = polymer.new_offset(OffsetKind::Add, 1, "H[2]O");
        assert_miette_snapshot!(invalid_composition);

        let missing_modification = polymer.modification(ModificationId(8));
        assert_eq!(missing_modification, None);
    }

    #[test]
    fn new_offset_with_composition() {
        let mut polymer_a = POLYMERIZER.new_polymer();
        let mut polymer_b = POLYMERIZER.new_polymer();
        let args = [
            (OffsetKind::Remove, 1, "H2O"),
            (OffsetKind::Add, 2, "H2O"),
            (OffsetKind::Remove, 3, "NH3+p"),
        ];

        for (kind, multiplier, formula) in args {
            let composition = ChemicalComposition::new(&ATOMIC_DB, formula).unwrap();

            let offset_a = polymer_a.new_offset(kind, multiplier, formula);
            let offset_b = polymer_b.new_offset_with_composition(kind, multiplier, composition);

            assert_eq!(offset_a, offset_b);
            assert_eq!(polymer_a, polymer_b);
        }
    }

    #[test]
    fn new_chain() {
        // NOTE: A macro since using a closure leads to borrow-checker issues...
        macro_rules! bond_chain {
            ($abbr:literal, $residues:expr) => {{
                let mut polymer = POLYMERIZER.new_polymer();
                polymer
                    .new_chain($abbr, $residues)
                    .map(|(residue_ids, bond_ids)| (polymer, residue_ids, bond_ids))
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
        assert!(
            polymer
                .residue(murnac)
                .unwrap()
                .group_state(&lactyl)
                .unwrap()
                .is_donor()
        );
        assert!(
            polymer
                .residue(alanine)
                .unwrap()
                .group_state(&n_terminal)
                .unwrap()
                .is_acceptor()
        );

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
        assert!(
            polymer
                .residue(alanine)
                .unwrap()
                .group_state(&c_terminal)
                .unwrap()
                .is_free()
        );

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
        assert!(
            polymer
                .residue(murnac)
                .unwrap()
                .group_state(&lactyl)
                .unwrap()
                .is_donor()
        );
        assert!(
            polymer
                .residue(alanine)
                .unwrap()
                .group_state(&n_terminal)
                .unwrap()
                .is_acceptor()
        );

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
        assert!(
            polymer
                .residue(alanine)
                .unwrap()
                .group_state(&c_terminal)
                .unwrap()
                .is_donor()
        );
        assert!(
            polymer
                .residue(lysine)
                .unwrap()
                .group_state(&sidechain)
                .unwrap()
                .is_acceptor()
        );
        // But the n_terminal of lysine remains untouched:
        assert!(
            polymer
                .residue(lysine)
                .unwrap()
                .group_state(&n_terminal)
                .unwrap()
                .is_free()
        );

        let groups_not_free = polymer.bond_groups("Stem", murnac, &lactyl, alanine, &n_terminal);
        assert_miette_snapshot!(groups_not_free);

        let glcnac = polymer.new_residue("g").unwrap();
        let invalid_bond = polymer.bond_groups("Pep", lysine, &c_terminal, glcnac, &n_terminal);
        assert_miette_snapshot!(invalid_bond);
        // When bonding fails due to the acceptor, make sure that the donor remains untouched
        assert!(
            polymer
                .residue(lysine)
                .unwrap()
                .group_state(&c_terminal)
                .unwrap()
                .is_free()
        );

        let nonexistent_bond = polymer.bond_groups("Super", murnac, &lactyl, alanine, &n_terminal);
        assert_miette_snapshot!(nonexistent_bond);
    }

    // NOTE: Needed to help type-inference figure out what `[]` is supposed to be
    const NO_RESIDUES: [ResidueId; 0] = [];

    #[test]
    fn bond_chain() {
        let mut polymer = POLYMERIZER.new_polymer();
        let residues = STEM_RESIDUES.map(|abbr| polymer.new_residue(abbr).unwrap());
        let polymer_without_bonds = polymer.clone();
        assert_polymer!(polymer, 515.24387164950, 515.51342919034656875);

        // Bonding chains of 1 or 0 residues is a no-op
        let new_bonds = polymer.bond_chain("Pep", NO_RESIDUES).unwrap();
        assert!(new_bonds.is_empty());
        assert_eq!(polymer, polymer_without_bonds);
        let new_bonds = polymer.bond_chain("Pep", &residues[..1]).unwrap();
        assert!(new_bonds.is_empty());
        assert_eq!(polymer, polymer_without_bonds);

        // Then actually bond the full chain (3 bonds for 4 residues)
        let bond_ids: Vec<_> = polymer.bond_chain("Pep", residues).unwrap();
        assert_eq!(bond_ids, vec![BondId(4), BondId(5), BondId(6)]);
        assert_polymer!(polymer, 461.21217759741, 461.46756989305707095);
        assert_ron_snapshot!(polymer, {
            ".**.composition, .**.lost, .**.gained" => "<FORMULA>",
            ".**.residues, .**.bonds" => insta::sorted_redaction(),
            ".**.functional_groups" => insta::sorted_redaction()
        });

        let nonexistent_bond = polymer.bond_chain("?", residues);
        assert_miette_snapshot!(nonexistent_bond);

        let nonexistent_residue = polymer.bond_chain("Pep", [ResidueId(3), ResidueId(4)]);
        assert_miette_snapshot!(nonexistent_residue);
    }

    #[test]
    #[allow(clippy::cognitive_complexity, clippy::too_many_lines)]
    fn localize_modification() {
        let mut polymer = POLYMERIZER.new_polymer();
        let murnac = polymer.new_residue("m").unwrap();
        assert_polymer!(polymer, 293.11106657336, 293.27091179713952985);

        let anhydro = polymer.new_modification(1, "Anh").unwrap();
        assert_eq!(anhydro, ModificationId(1));
        assert_polymer!(polymer, 275.10050188933, 275.25562536470969725);
        for (_, group_state) in polymer.residue(murnac).unwrap().functional_groups() {
            assert!(group_state.is_free());
        }

        let moved_anhydro = polymer.localize_modification(anhydro, murnac).unwrap();
        assert_eq!(moved_anhydro, vec![ModificationId(2)]);
        assert_eq!(polymer.modification(anhydro), None);
        assert_polymer!(polymer, 275.10050188933, 275.25562536470969725);
        let reducing_end = FunctionalGroup::new("Hydroxyl", "Reducing End");
        assert!(
            polymer
                .residue(murnac)
                .unwrap()
                .group_state(&reducing_end)
                .unwrap()
                .is_modified()
        );

        let triple_calcium = polymer.new_modification(3, "Ca").unwrap();
        assert_eq!(triple_calcium, ModificationId(3));
        assert_polymer!(
            polymer,
            391.963153599012610,
            392.46457201844549725,
            3,
            130.65438453300420333333333333,
            130.82152400614849908333333333
        );
        for (&group, group_state) in polymer.residue(murnac).unwrap().functional_groups() {
            if group != reducing_end {
                assert!(group_state.is_free());
            }
        }

        let moved_triple_calcium = polymer
            .localize_modification(triple_calcium, murnac)
            .unwrap();
        assert_eq!(
            moved_triple_calcium,
            vec![ModificationId(4), ModificationId(5), ModificationId(6)]
        );
        assert_eq!(polymer.modification(triple_calcium), None);
        assert_polymer!(
            polymer,
            391.963153599012610,
            392.46457201844549725,
            3,
            130.65438453300420333333333333,
            130.82152400614849908333333333
        );
        let n_acetyl = FunctionalGroup::new("Acetyl", "Secondary Amide");
        for (&group, group_state) in polymer.residue(murnac).unwrap().functional_groups() {
            if group != n_acetyl {
                assert!(group_state.is_modified());
            }
        }

        let phosphate = polymer.new_offset(OffsetKind::Add, 1, "PO4+3e").unwrap();
        assert_eq!(phosphate, ModificationId(7));
        assert_polymer!(polymer, 486.918219815439805, 487.43759945386580385);
        assert_eq!(
            polymer
                .residue(murnac)
                .unwrap()
                .offset_modifications()
                .count(),
            0
        );

        let moved_phosphate = polymer.localize_modification(phosphate, murnac).unwrap();
        assert_eq!(moved_phosphate, vec![ModificationId(8)]);
        assert_eq!(polymer.modification(phosphate), None);
        assert_polymer!(polymer, 486.918219815439805, 487.43759945386580385);
        assert_eq!(
            polymer
                .residue(murnac)
                .unwrap()
                .offset_modifications()
                .count(),
            1
        );

        let quadruple_hydrogen_loss = polymer.new_offset(OffsetKind::Remove, 4, "H").unwrap();
        assert_eq!(quadruple_hydrogen_loss, ModificationId(9));
        assert_polymer!(polymer, 482.886919686519805, 483.40583643764269445);
        assert_eq!(
            polymer
                .residue(murnac)
                .unwrap()
                .offset_modifications()
                .count(),
            1
        );

        let moved_quadruple_hydrogen_loss = polymer
            .localize_modification(quadruple_hydrogen_loss, murnac)
            .unwrap();
        assert_eq!(moved_quadruple_hydrogen_loss, vec![ModificationId(10)]);
        assert_eq!(polymer.modification(quadruple_hydrogen_loss), None);
        assert_polymer!(polymer, 482.886919686519805, 483.40583643764269445);
        assert_eq!(
            polymer
                .residue(murnac)
                .unwrap()
                .offset_modifications()
                .count(),
            2
        );

        let deleted_modification = polymer.localize_modification(phosphate, murnac);
        assert_miette_snapshot!(deleted_modification);
        // Make sure the polymer remains untouched!
        assert_polymer!(polymer, 482.886919686519805, 483.40583643764269445);

        let already_localized_offset = polymer.localize_modification(moved_phosphate[0], murnac);
        assert_miette_snapshot!(already_localized_offset);
        // Make sure the polymer remains untouched!
        assert_polymer!(polymer, 482.886919686519805, 483.40583643764269445);

        let already_localized_modification =
            polymer.localize_modification(moved_anhydro[0], murnac);
        assert_miette_snapshot!(already_localized_modification);
        // Make sure the polymer remains untouched!
        assert_polymer!(polymer, 482.886919686519805, 483.40583643764269445);
    }

    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn modify_polymer() {
        let mut no_targets = POLYMERIZER.new_polymer();
        for abbr in STEM_RESIDUES {
            no_targets.new_residue(abbr).unwrap();
        }
        assert_polymer!(no_targets, 515.24387164950, 515.51342919034656875);

        let modification_ids = no_targets.modify_polymer("Met").unwrap();
        assert_eq!(modification_ids.len(), 0);
        assert_polymer!(no_targets, 515.24387164950, 515.51342919034656875);

        let mut one_residue = POLYMERIZER.new_polymer();
        one_residue.new_residue("m").unwrap();
        assert_polymer!(one_residue, 293.11106657336, 293.27091179713952985);

        let modification_ids = one_residue.modify_polymer("Met").unwrap();
        assert_eq!(modification_ids.len(), 3);
        assert_polymer!(one_residue, 335.15801676674, 335.35076401167994095);

        let mut two_residues = POLYMERIZER.new_polymer();
        let glcnac = two_residues.new_residue("g").unwrap();
        let murnac = two_residues.new_residue("m").unwrap();
        assert_polymer!(two_residues, 514.20100377866, 514.47904303921364750);

        let bond_id = two_residues.bond_residues("Gly", glcnac, murnac).unwrap();
        assert_polymer!(two_residues, 496.19043909463, 496.46375660678381490);

        let modification_ids = two_residues.modify_polymer("Met").unwrap();
        assert_eq!(modification_ids.len(), 4);
        assert_polymer!(two_residues, 552.25303935247, 552.57022622617102970);

        let groups = |id| two_residues.residue(id).unwrap().functional_groups();

        for (fg, &gs) in groups(glcnac) {
            match gs {
                GroupState::Free => assert_eq!(fg.name, "Acetyl"),
                GroupState::Modified(_) => assert_eq!(fg.name, "Hydroxyl"),
                GroupState::Donor(id) => {
                    assert_eq!(id, bond_id);
                    assert_eq!(fg.location, "Reducing End");
                }
                GroupState::Acceptor(_) => panic!(),
            }
        }

        for (fg, &gs) in groups(murnac) {
            match gs {
                GroupState::Free => assert_ne!(fg.name, "Hydroxyl"),
                GroupState::Modified(_) => assert_eq!(fg.name, "Hydroxyl"),
                GroupState::Donor(_) => panic!(),
                GroupState::Acceptor(id) => {
                    assert_eq!(id, bond_id);
                    assert_eq!(fg.location, "Nonreducing End");
                }
            }
        }

        let nonexistent_modification = two_residues.modify_polymer("Arg");
        assert_miette_snapshot!(nonexistent_modification);
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
        assert!(
            polymer
                .residue(murnac)
                .unwrap()
                .group_state(&reducing_end)
                .unwrap()
                .is_modified()
        );

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
    #[allow(clippy::cognitive_complexity)]
    fn modify_only_groups() {
        let mut polymer = POLYMERIZER.new_polymer();
        let murnac = polymer.new_residue("m").unwrap();
        assert_polymer!(polymer, 293.11106657336, 293.27091179713952985);

        let modification_ids = polymer.modify_only_groups("Met", murnac, 3).unwrap();
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
            assert!(
                polymer
                    .residue(murnac)
                    .unwrap()
                    .group_state(&hydroxyl_group)
                    .unwrap()
                    .is_modified()
            );
        }

        polymer.remove_residue(murnac);
        let murnac = polymer.new_residue("m").unwrap();
        assert_polymer!(polymer, 293.11106657336, 293.27091179713952985);

        let modification_ids = polymer.modify_only_groups("Ca", murnac, 4).unwrap();
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
        let modification_ids = polymer.modify_only_groups("Ca", alanine, 2).unwrap();
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
        let modification_ids = polymer.modify_only_groups("Met", murnac, 3).unwrap();
        assert_eq!(
            modification_ids,
            vec![ModificationId(13), ModificationId(14), ModificationId(15)]
        );
        let all_groups_occupied = polymer.modify_only_groups("Ca", murnac, 4);
        assert_miette_snapshot!(all_groups_occupied);
        let still_all_groups_occupied = polymer.modify_only_groups("Ca", murnac, 2);
        assert_miette_snapshot!(still_all_groups_occupied);

        assert_polymer!(polymer, 335.15801676674, 335.35076401167994095);
        let modification_ids = polymer.modify_only_groups("Ca", murnac, 1).unwrap();
        assert_eq!(modification_ids, vec![ModificationId(16)]);
        assert_polymer!(
            polymer,
            374.112234003300870,
            374.42041289625854095,
            1,
            374.112234003300870,
            374.42041289625854095
        );

        let residue_not_in_polymer = polymer.modify_only_groups("Ca", ResidueId(42), 2);
        assert_miette_snapshot!(residue_not_in_polymer);

        let murnac = polymer.new_residue("m").unwrap();
        let no_matching_groups = polymer.modify_only_groups("Anh", murnac, 2);
        assert_miette_snapshot!(no_matching_groups);
    }

    #[test]
    fn modify_group() {
        let mut polymer = POLYMERIZER.new_polymer();
        let murnac = polymer.new_residue("m").unwrap();
        assert_polymer!(polymer, 293.11106657336, 293.27091179713952985);

        let reducing_end = FunctionalGroup::new("Hydroxyl", "Reducing End");
        let modification_id = polymer.modify_group("Met", murnac, &reducing_end).unwrap();
        assert_eq!(modification_id, ModificationId(1));
        assert_polymer!(polymer, 307.12671663782, 307.29752920198633355);
        assert!(
            polymer
                .residue(murnac)
                .unwrap()
                .group_state(&reducing_end)
                .unwrap()
                .is_modified()
        );

        let modify_non_free_group = polymer.modify_group("Anh", murnac, &reducing_end);
        assert_miette_snapshot!(modify_non_free_group);

        let residue_from_wrong_polymer = polymer.modify_group("Anh", ResidueId(1), &reducing_end);
        assert_miette_snapshot!(residue_from_wrong_polymer);

        let invalid_group = polymer.modify_group("Ac", murnac, &reducing_end);
        assert_miette_snapshot!(invalid_group);

        let alanine = polymer.new_residue("A").unwrap();
        let nonexistent_group = polymer.modify_group("Red", alanine, &reducing_end);
        assert_miette_snapshot!(nonexistent_group);

        let nonexistent_modification = polymer.modify_group("Arg", murnac, &reducing_end);
        assert_miette_snapshot!(nonexistent_modification);
    }

    #[test]
    fn modify_groups() {
        let mut polymer = POLYMERIZER.new_polymer();
        let murnac = polymer.new_residue("m").unwrap();
        assert_polymer!(polymer, 293.11106657336, 293.27091179713952985);

        let reducing_end = FunctionalGroup::new("Hydroxyl", "Reducing End");
        let nonreducing_end = FunctionalGroup::new("Hydroxyl", "Nonreducing End");
        let six_position = FunctionalGroup::new("Hydroxyl", "6-Position");

        let modification_ids = polymer
            .modify_groups("Met", murnac, [nonreducing_end, six_position])
            .unwrap();
        assert_eq!(modification_ids, vec![ModificationId(1), ModificationId(2)]);
        assert_polymer!(polymer, 321.14236670228, 321.32414660683313725);
        assert!(
            polymer
                .residue(murnac)
                .unwrap()
                .group_state(&reducing_end)
                .unwrap()
                .is_free()
        );
        assert!(
            polymer
                .residue(murnac)
                .unwrap()
                .group_state(&nonreducing_end)
                .unwrap()
                .is_modified()
        );
        assert!(
            polymer
                .residue(murnac)
                .unwrap()
                .group_state(&six_position)
                .unwrap()
                .is_modified()
        );

        let modify_non_free_group =
            polymer.modify_groups("Ca", murnac, [reducing_end, six_position]);
        assert_miette_snapshot!(modify_non_free_group);

        let modify_invalid_group =
            polymer.modify_groups("Anh", murnac, [reducing_end, six_position]);
        assert_miette_snapshot!(modify_invalid_group);

        let modification_ids = polymer
            .modify_groups("Anh", murnac, [reducing_end])
            .unwrap();
        assert_eq!(modification_ids, vec![ModificationId(3)]);
        assert_polymer!(polymer, 303.13180201825, 303.30886017440330465);
        assert!(
            polymer
                .residue(murnac)
                .unwrap()
                .group_state(&reducing_end)
                .unwrap()
                .is_modified()
        );

        let modify_multiple_errors =
            polymer.modify_groups("Anh", murnac, [reducing_end, six_position]);
        assert_miette_snapshot!(modify_multiple_errors);

        let residue_from_wrong_polymer =
            polymer.modify_groups("Anh", ResidueId(1337), [reducing_end, six_position]);
        assert_miette_snapshot!(residue_from_wrong_polymer);

        let alanine = polymer.new_residue("A").unwrap();
        let nonexistent_group = polymer.modify_groups("Red", alanine, [reducing_end, six_position]);
        assert_miette_snapshot!(nonexistent_group);

        let nonexistent_modification = polymer.modify_groups("Arg", murnac, [reducing_end]);
        assert_miette_snapshot!(nonexistent_modification);
    }

    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn offset_residue() {
        let mut polymer = POLYMERIZER.new_polymer();
        let alanine = polymer.new_residue("A").unwrap();
        assert_polymer!(polymer, 89.04767846918, 89.09330602867854225);

        let modification_id = polymer
            .offset_residue(OffsetKind::Add, 1, "p", alanine)
            .unwrap();
        assert_eq!(modification_id, ModificationId(1));
        assert_polymer!(
            polymer,
            90.054954935801,
            90.10058249529954225,
            1,
            90.054954935801,
            90.10058249529954225
        );

        let modification_id = polymer
            .offset_residue(OffsetKind::Add, 2, "e", alanine)
            .unwrap();
        assert_eq!(modification_id, ModificationId(2));
        assert_polymer!(
            polymer,
            90.056052095619130,
            90.10167965511767225,
            -1,
            90.056052095619130,
            90.10167965511767225
        );

        let murnac = polymer.new_residue("m").unwrap();
        let modification_id = polymer
            .offset_residue(OffsetKind::Add, 1, "H2O+p", murnac)
            .unwrap();
        assert_eq!(modification_id, ModificationId(4));
        assert_polymer!(polymer, 402.184959819630130, 402.39515435130803470);

        let modification_id = polymer
            .offset_residue(OffsetKind::Remove, 4, "CO", alanine)
            .unwrap();
        assert_eq!(modification_id, ModificationId(5));
        assert_polymer!(polymer, 290.205301341350130, 290.35459106709392710);

        let mut alanine_offsets: Vec<_> = polymer
            .residue(alanine)
            .unwrap()
            .offset_modifications()
            .collect();
        alanine_offsets.sort_unstable();
        assert_eq!(
            alanine_offsets,
            vec![ModificationId(1), ModificationId(2), ModificationId(5)]
        );

        let mut murnac_offsets: Vec<_> = polymer
            .residue(murnac)
            .unwrap()
            .offset_modifications()
            .collect();
        murnac_offsets.sort_unstable();
        assert_eq!(murnac_offsets, vec![ModificationId(4)]);

        let zero_multiplier = polymer.offset_residue(OffsetKind::Remove, 0, "H", alanine);
        assert_miette_snapshot!(zero_multiplier);

        let invalid_composition = polymer.offset_residue(OffsetKind::Add, 1, "H[2]O", alanine);
        assert_miette_snapshot!(invalid_composition);

        let residue_from_wrong_polymer =
            polymer.offset_residue(OffsetKind::Add, 1, "H2O", ResidueId(4));
        assert_miette_snapshot!(residue_from_wrong_polymer);
    }

    #[test]
    fn offset_residue_with_composition() {
        let mut polymer_a = POLYMERIZER.new_polymer();
        let mut polymer_b = POLYMERIZER.new_polymer();
        let args = [
            (OffsetKind::Add, 1, "p"),
            (OffsetKind::Add, 2, "e"),
            (OffsetKind::Add, 1, "H2O+p"),
            (OffsetKind::Remove, 4, "CO"),
        ];

        let alanine_a = polymer_a.new_residue("A").unwrap();
        let alanine_b = polymer_b.new_residue("A").unwrap();

        for (kind, multiplier, formula) in args {
            let composition = ChemicalComposition::new(&ATOMIC_DB, formula).unwrap();

            let offset_a = polymer_a.offset_residue(kind, multiplier, formula, alanine_a);
            let offset_b =
                polymer_b.offset_residue_with_composition(kind, multiplier, composition, alanine_b);

            assert_eq!(offset_a, offset_b);
            assert_eq!(polymer_a, polymer_b);
        }
    }
}
