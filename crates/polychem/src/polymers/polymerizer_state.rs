use std::{
    borrow::{Borrow, Cow},
    cmp::Ordering,
    collections::hash_map::Entry,
    mem,
};

use ahash::{HashMap, HashMapExt, HashSet, HashSetExt};

use crate::{
    Count, FunctionalGroup, GroupState, Id, Residue, ResidueGroup, ResidueId,
    moieties::target::{Index, Target},
};

use super::{errors::FindFreeGroupsError, polymerizer::Polymerizer};

#[derive(Clone, Eq, PartialEq, Debug)]
// NOTE: This whole module is `pub(crate)`, so that's how we can get away with `pub` here without leaking this struct
// into the public API!
pub struct PolymerizerState<'a, 'p> {
    // NOTE: Currently `pub` because this struct is already private, and there isn't any fancy logic I need to attach
    // to the `polymerizer` field. Definitely private this field if `PolymerizerState` is ever made public!
    pub polymerizer: Polymerizer<'a, 'p>,
    next_id: Id,
    // PERF: `GroupState` is currently 16-times larger than a `bool`, and this index is only really "essential" for
    // tracking free groups, so some space could be saved by moving back to a `bool`. With that being said, tracking the
    // full `GroupState` makes the error reporting code a bit cleaner, and storing a `GroupState` is both clearer than a
    // `bool` and less messy than a second `IsFree` enum, so this seems like a worthy compromise until benchmarks can
    // prove otherwise!
    group_index: Index<'p, HashMap<ResidueId, GroupState>>,
}

impl<'a, 'p> PolymerizerState<'a, 'p> {
    #[must_use]
    pub fn new(polymerizer: &Polymerizer<'a, 'p>) -> Self {
        Self {
            polymerizer: *polymerizer,
            next_id: Id::default(),
            group_index: Index::new(),
        }
    }

    #[must_use]
    // NOTE: I don't want to promise that this will always be `const` in the API yet...
    #[expect(clippy::missing_const_for_fn)]
    pub fn next_id(&mut self) -> Id {
        let current_id = self.next_id;
        self.next_id += 1;
        current_id
    }

    pub fn index_residue_groups(&mut self, id: ResidueId, residue: &Residue<'a, 'p>) {
        for (group, &state) in residue.functional_groups() {
            let target = Target::from_residue_and_group(residue, group);
            self.group_index
                .entry(target)
                .or_default()
                .insert(id, state);
        }
    }

    pub fn unindex_residue_groups(&mut self, id: ResidueId, residue: &Residue<'a, 'p>) {
        for (group, _) in residue.functional_groups() {
            let target = Target::from_residue_and_group(residue, group);
            if let Some(id_states) = self.group_index.get_mut(target) {
                id_states.remove(&id);
                if id_states.is_empty() {
                    self.group_index.remove(target);
                }
            }
        }
    }

    pub fn update_group_index(
        &mut self,
        target: impl Into<Target<&'p str>>,
        residue: ResidueId,
        state: GroupState,
    ) -> Option<GroupState> {
        let residue_state = self.group_index.get_mut(target)?.get_mut(&residue)?;
        Some(mem::replace(residue_state, state))
    }

    pub fn find_any_free_groups<'t, T: 'p>(
        &self,
        targets: &'t [T],
        residue: ResidueId,
        number: impl TryInto<Count>,
    ) -> Result<HashSet<FunctionalGroup<'p>>, FindFreeGroupsError>
    where
        &'t T: Into<Target<&'p str>>,
    {
        let number = number
            .try_into()
            .map_err(|_| FindFreeGroupsError::ZeroGroupNumber)?
            .into();

        // PERF: There are likely some small performance gains to be made by moving to `Vec`s instead of `HashSet`s, but
        // using a `HashSet` makes it clearer that the returned groups are unique!
        let free_groups: HashSet<_> = self.free_residue_groups(targets, residue).collect();
        match free_groups.len().cmp(&number) {
            Ordering::Less => {
                Err(self.diagnose_any_missing_groups(targets, residue, number, &free_groups))
            }
            Ordering::Equal => Ok(free_groups),
            Ordering::Greater => Err(FindFreeGroupsError::ambiguous_groups(
                residue,
                number,
                &free_groups,
            )),
        }
    }

    pub fn find_these_free_groups<'t, T: 'p>(
        &self,
        targets: &'t [T],
        residue: ResidueId,
        // FIXME: It's a bit odd that this is a generic, whilst `targets` isn't, but that's just to make calling this
        // from `Polymer` easier... Worth another look...
        groups: impl IntoIterator<Item: Borrow<FunctionalGroup<'p>>>,
    ) -> Result<HashSet<FunctionalGroup<'p>>, FindFreeGroupsError>
    where
        &'t T: Into<Target<&'p str>>,
    {
        let mut target_groups = HashSet::new();
        for group in groups {
            let group = group.borrow().to_owned();
            if !target_groups.insert(group) {
                return Err(FindFreeGroupsError::duplicate_target_group(&group));
            }
        }

        if target_groups.is_empty() {
            return Err(FindFreeGroupsError::NoTargetGroups);
        }

        let free_groups: HashSet<_> = self.free_residue_groups(targets, residue).collect();
        if target_groups.is_subset(&free_groups) {
            Ok(target_groups)
        } else {
            let missing_groups = target_groups.difference(&free_groups).copied();
            Err(self.diagnose_these_missing_groups(targets, residue, missing_groups))
        }
    }

    pub fn polymer_groups<'t, T: 'p>(
        &self,
        targets: &'t [T],
    ) -> impl Iterator<Item = (FunctionalGroup<'p>, Cow<'_, HashMap<ResidueId, GroupState>>)>
    where
        &'t T: Into<Target<&'p str>>,
    {
        let polymer_groups = targets.iter().flat_map(|target| {
            self.group_index.matches(target, |g, l, _, m| {
                // SAFETY: `.matches()` always returns a complete `Target` with no `None` fields, so `.unwrap()` is safe
                // NOTE: The `Cow` is used here since normally we won't need to mutate any of the underlying `HashMap`s,
                // and a reference is just fine, but if there are some duplicate functional groups, then their
                // `HashMap`s will need to be merged, and we'll need to create a copy that's mutable!
                (FunctionalGroup::new(g, l.unwrap()), Cow::Borrowed(m))
            })
        });

        let mut deduplicated_groups = HashMap::new();
        for (group, id_states) in polymer_groups {
            match deduplicated_groups.entry(group) {
                Entry::Occupied(mut e) => {
                    // NOTE: Type inference gets a bit confused here, so I need to explicitly mention that I've got a
                    // `&mut Cow<HashMap<_, _>>` here! Perhaps someday this could be simplified down to something like:
                    // `e.get_mut().to_mut().extend(id_states)`
                    let existing_id_states: &mut Cow<HashMap<_, _>> = e.get_mut();
                    existing_id_states.to_mut().extend(&*id_states);
                }
                Entry::Vacant(e) => {
                    e.insert(id_states);
                }
            }
        }
        deduplicated_groups.into_iter()
    }

    pub fn free_polymer_groups<'t, T: 'p>(
        &self,
        targets: &'t [T],
    ) -> impl Iterator<Item = ResidueGroup<'p>>
    where
        &'t T: Into<Target<&'p str>>,
    {
        self.polymer_groups(targets)
            .flat_map(|(fg, ids)| -> Vec<_> {
                // FIXME: I really wish I didn't need `.collect()` here...
                ids.iter()
                    .filter_map(move |(&id, gs)| gs.is_free().then_some(ResidueGroup(id, fg)))
                    .collect()
            })
    }

    pub fn residue_groups<'t, T: 'p>(
        &self,
        targets: &'t [T],
        residue: ResidueId,
    ) -> impl Iterator<Item = (FunctionalGroup<'p>, GroupState)>
    where
        &'t T: Into<Target<&'p str>>,
    {
        self.polymer_groups(targets)
            .filter_map(move |(fg, ids)| ids.get(&residue).map(|&gs| (fg, gs)))
    }

    pub fn free_residue_groups<'t, T: 'p>(
        &self,
        targets: &'t [T],
        residue: ResidueId,
    ) -> impl Iterator<Item = FunctionalGroup<'p>>
    where
        &'t T: Into<Target<&'p str>>,
    {
        self.residue_groups(targets, residue)
            .filter_map(|(fg, gs)| gs.is_free().then_some(fg))
    }
}

// Private Methods =====================================================================================================

impl<'p> PolymerizerState<'_, 'p> {
    fn diagnose_any_missing_groups<'t, T: 'p>(
        &self,
        targets: &'t [T],
        residue: ResidueId,
        number: usize,
        free_groups: &HashSet<FunctionalGroup<'p>>,
    ) -> FindFreeGroupsError
    where
        &'t T: Into<Target<&'p str>>,
    {
        // MISSING: This method does not and cannot check if the provided `ResidueId` is actually a valid ID in the
        // current polymer! If you pass it an invalid ID, then it will just return the generic "too-few groups" error.
        // It is expected that the validity of `ResidueId`s is checked *before* they are passed to this method!
        let non_free_groups: Vec<_> = self.residue_groups(targets, residue).collect();

        if non_free_groups.len() >= number {
            FindFreeGroupsError::groups_occupied(residue, number, non_free_groups)
        } else {
            FindFreeGroupsError::too_few_matching_groups(residue, targets, number, free_groups)
        }
    }

    fn diagnose_these_missing_groups<'t, T: 'p>(
        &self,
        targets: &'t [T],
        residue: ResidueId,
        missing_groups: impl Iterator<Item = FunctionalGroup<'p>>,
    ) -> FindFreeGroupsError
    where
        &'t T: Into<Target<&'p str>>,
    {
        // MISSING: This method does not and cannot check if the provided `ResidueId` is actually a valid ID in the
        // current polymer! It is expected that the validity of `ResidueId`s is checked *before* they are passed to this
        // method!
        let mut errors: Vec<_> = missing_groups
            .map(|group| {
                let current_target = Target::from(group);

                let theoretically_possible = targets.iter().any(|possible_target| {
                    // NOTE: The `residue` field of each target is cleared since it's impossible for the
                    // `current_target` (which contains only a `group` and `location`) to ever match a "complete"
                    // target that specifies a `residue`
                    let Target {
                        group, location, ..
                    } = possible_target.into();
                    current_target.matches(&Target::new(group, location, None))
                });
                // NOTE: Utterly perplexed as to why I need to add the turbofish — `::<Target<&str>>` — below...
                // Without it, Rust expects `current_target` to be `T`? Perhaps calling a self method means inheriting
                // the caller's generic parameters? Feels like a borderline compiler bug...
                let non_free_group = self
                    .residue_groups::<Target<&str>>(&[current_target], residue)
                    .next();

                if !theoretically_possible {
                    FindFreeGroupsError::invalid_target_group(targets, current_target)
                } else if let Some((group, state)) = non_free_group {
                    FindFreeGroupsError::group_occupied(&group, residue, state)
                } else {
                    FindFreeGroupsError::nonexistent_group(&group, residue)
                }
            })
            .collect();
        if errors.len() == 1 {
            // SAFETY: We've just checked that there is one item is present in the Vec, so this will never panic
            errors.pop().unwrap()
        } else {
            FindFreeGroupsError::multiple_missing_free_groups(errors)
        }
    }
}

#[cfg(test)]
mod tests {
    use std::iter::zip;

    use std::sync::LazyLock;

    use crate::{
        AtomicDatabase, BondId, ModificationId, PolymerDatabase, Residue, moieties::target::Target,
        testing_tools::assert_miette_snapshot,
    };

    use super::*;

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

    #[test]
    fn next_id() {
        let mut state = PolymerizerState::new(&POLYMERIZER);
        assert_eq!(state.next_id(), 0);
        assert_eq!(state.next_id(), 1);
        assert_eq!(state.next_id(), 2);
        assert_eq!(state.next_id(), 3);
        assert_eq!(state.next_id(), 4);
    }

    #[test]
    fn index_residue_groups() {
        let mut state = PolymerizerState::new(&POLYMERIZER);
        let alanine = Residue::new(&POLYMER_DB, "A").unwrap();
        let lysine = Residue::new(&POLYMER_DB, "K").unwrap();
        state.index_residue_groups(ResidueId(0), &alanine);
        state.index_residue_groups(ResidueId(1), &alanine);
        state.index_residue_groups(ResidueId(2), &lysine);

        macro_rules! assert_targeted_groups {
            ($target:expr, $output:expr) => {
                let mut groups: Vec<_> = state
                    .group_index
                    .matches($target, |g, l, r, m| {
                        let mut id_states: Vec<_> = m.iter().map(|(&i, &f)| (i, f)).collect();
                        id_states.sort_unstable();
                        (Target::new(g, l, r), id_states)
                    })
                    .collect();
                groups.sort_unstable();
                assert_eq!(groups, $output)
            };
        }

        let amino = Target::new("Amino", None, None);
        let amino_groups = vec![
            (
                Target::new("Amino", Some("N-Terminal"), Some("Alanine")),
                vec![
                    (ResidueId(0), GroupState::Free),
                    (ResidueId(1), GroupState::Free),
                ],
            ),
            (
                Target::new("Amino", Some("N-Terminal"), Some("Lysine")),
                vec![(ResidueId(2), GroupState::Free)],
            ),
            (
                Target::new("Amino", Some("Sidechain"), Some("Lysine")),
                vec![(ResidueId(2), GroupState::Free)],
            ),
        ];
        assert_targeted_groups!(amino, amino_groups);

        let n_terminal = Target::new("Amino", Some("N-Terminal"), None);
        let n_terminal_groups = vec![
            (
                Target::new("Amino", Some("N-Terminal"), Some("Alanine")),
                vec![
                    (ResidueId(0), GroupState::Free),
                    (ResidueId(1), GroupState::Free),
                ],
            ),
            (
                Target::new("Amino", Some("N-Terminal"), Some("Lysine")),
                vec![(ResidueId(2), GroupState::Free)],
            ),
        ];
        assert_targeted_groups!(n_terminal, n_terminal_groups);

        let carboxyl = Target::new("Carboxyl", None, None);
        let carboxyl_groups = vec![
            (
                Target::new("Carboxyl", Some("C-Terminal"), Some("Alanine")),
                vec![
                    (ResidueId(0), GroupState::Free),
                    (ResidueId(1), GroupState::Free),
                ],
            ),
            (
                Target::new("Carboxyl", Some("C-Terminal"), Some("Lysine")),
                vec![(ResidueId(2), GroupState::Free)],
            ),
        ];
        assert_targeted_groups!(carboxyl, carboxyl_groups);
    }

    #[test]
    fn unindex_residue_groups() {
        let mut state = PolymerizerState::new(&POLYMERIZER);
        let alanine = Residue::new(&POLYMER_DB, "A").unwrap();
        let lysine = Residue::new(&POLYMER_DB, "K").unwrap();
        state.index_residue_groups(ResidueId(0), &alanine);
        state.index_residue_groups(ResidueId(1), &alanine);
        state.index_residue_groups(ResidueId(2), &lysine);

        macro_rules! assert_targeted_groups {
            ($target:expr, $output:expr) => {
                let mut groups: Vec<_> = state
                    .group_index
                    .matches($target, |g, l, r, m| {
                        let mut id_states: Vec<_> = m.iter().map(|(&i, &f)| (i, f)).collect();
                        id_states.sort_unstable();
                        (Target::new(g, l, r), id_states)
                    })
                    .collect();
                groups.sort_unstable();
                assert_eq!(groups, $output)
            };
        }

        let amino = Target::new("Amino", None, None);
        let amino_groups = vec![
            (
                Target::new("Amino", Some("N-Terminal"), Some("Alanine")),
                vec![
                    (ResidueId(0), GroupState::Free),
                    (ResidueId(1), GroupState::Free),
                ],
            ),
            (
                Target::new("Amino", Some("N-Terminal"), Some("Lysine")),
                vec![(ResidueId(2), GroupState::Free)],
            ),
            (
                Target::new("Amino", Some("Sidechain"), Some("Lysine")),
                vec![(ResidueId(2), GroupState::Free)],
            ),
        ];
        assert_targeted_groups!(amino, amino_groups);

        let carboxyl = Target::new("Carboxyl", None, None);
        let carboxyl_groups = vec![
            (
                Target::new("Carboxyl", Some("C-Terminal"), Some("Alanine")),
                vec![
                    (ResidueId(0), GroupState::Free),
                    (ResidueId(1), GroupState::Free),
                ],
            ),
            (
                Target::new("Carboxyl", Some("C-Terminal"), Some("Lysine")),
                vec![(ResidueId(2), GroupState::Free)],
            ),
        ];
        assert_targeted_groups!(carboxyl, carboxyl_groups);

        // And then remove one residue from the index
        state.unindex_residue_groups(ResidueId(2), &lysine);

        let amino_groups = vec![(
            Target::new("Amino", Some("N-Terminal"), Some("Alanine")),
            vec![
                (ResidueId(0), GroupState::Free),
                (ResidueId(1), GroupState::Free),
            ],
        )];
        assert_targeted_groups!(amino, amino_groups);

        let carboxyl_groups = vec![(
            Target::new("Carboxyl", Some("C-Terminal"), Some("Alanine")),
            vec![
                (ResidueId(0), GroupState::Free),
                (ResidueId(1), GroupState::Free),
            ],
        )];
        assert_targeted_groups!(carboxyl, carboxyl_groups);

        // And then remove one more
        state.unindex_residue_groups(ResidueId(0), &alanine);

        let amino_groups = vec![(
            Target::new("Amino", Some("N-Terminal"), Some("Alanine")),
            vec![(ResidueId(1), GroupState::Free)],
        )];
        assert_targeted_groups!(amino, amino_groups);

        let carboxyl_groups = vec![(
            Target::new("Carboxyl", Some("C-Terminal"), Some("Alanine")),
            vec![(ResidueId(1), GroupState::Free)],
        )];
        assert_targeted_groups!(carboxyl, carboxyl_groups);

        // Duplicate removals accomplish nothing
        state.unindex_residue_groups(ResidueId(0), &alanine);
        assert_targeted_groups!(amino, amino_groups);
        assert_targeted_groups!(carboxyl, carboxyl_groups);

        // But we can get rid of the last one
        state.unindex_residue_groups(ResidueId(1), &alanine);
        assert_targeted_groups!(amino, Vec::new());
        assert_targeted_groups!(carboxyl, Vec::new());
    }

    #[test]
    fn update_group_index() {
        let mut state = PolymerizerState::new(&POLYMERIZER);
        let alanine = Residue::new(&POLYMER_DB, "A").unwrap();
        let lysine = Residue::new(&POLYMER_DB, "K").unwrap();
        state.index_residue_groups(ResidueId(0), &alanine);
        state.index_residue_groups(ResidueId(1), &alanine);
        state.index_residue_groups(ResidueId(2), &lysine);

        macro_rules! assert_targeted_groups {
            ($target:expr, $output:expr) => {
                let mut groups: Vec<_> = state
                    .group_index
                    .matches($target, |_, _, _, m| m.iter().map(|(&r, &g)| (r, g)))
                    .flatten()
                    .collect();
                groups.sort_unstable();
                assert_eq!(groups, $output)
            };
        }

        let n_terminal = |residue| Target::new("Amino", Some("N-Terminal"), residue);
        let n_terminal_groups = vec![
            (ResidueId(0), GroupState::Free),
            (ResidueId(1), GroupState::Free),
            (ResidueId(2), GroupState::Free),
        ];
        assert_targeted_groups!(n_terminal(None), n_terminal_groups);

        let previous_state = state.update_group_index(
            n_terminal(Some("Alanine")),
            ResidueId(0),
            GroupState::Modified(ModificationId(42)),
        );
        let n_terminal_groups = vec![
            (ResidueId(0), GroupState::Modified(ModificationId(42))),
            (ResidueId(1), GroupState::Free),
            (ResidueId(2), GroupState::Free),
        ];
        assert_eq!(previous_state, Some(GroupState::Free));
        assert_targeted_groups!(n_terminal(None), n_terminal_groups);

        let previous_state = state.update_group_index(
            n_terminal(Some("Lysine")),
            ResidueId(2),
            GroupState::Acceptor(BondId(314)),
        );
        let n_terminal_groups = vec![
            (ResidueId(0), GroupState::Modified(ModificationId(42))),
            (ResidueId(1), GroupState::Free),
            (ResidueId(2), GroupState::Acceptor(BondId(314))),
        ];
        assert_eq!(previous_state, Some(GroupState::Free));
        assert_targeted_groups!(n_terminal(None), n_terminal_groups);

        let previous_state = state.update_group_index(
            n_terminal(Some("Alanine")),
            ResidueId(0),
            GroupState::Donor(BondId(42)),
        );
        let n_terminal_groups = vec![
            (ResidueId(0), GroupState::Donor(BondId(42))),
            (ResidueId(1), GroupState::Free),
            (ResidueId(2), GroupState::Acceptor(BondId(314))),
        ];
        assert_eq!(
            previous_state,
            Some(GroupState::Modified(ModificationId(42)))
        );
        assert_targeted_groups!(n_terminal(None), n_terminal_groups);

        let c_terminal = |residue| Target::new("Carboxyl", Some("C-Terminal"), residue);
        let previous_state = state.update_group_index(
            c_terminal(Some("Alanine")),
            ResidueId(1),
            GroupState::Donor(BondId(1337)),
        );
        let c_terminal_groups = vec![
            (ResidueId(0), GroupState::Free),
            (ResidueId(1), GroupState::Donor(BondId(1337))),
            (ResidueId(2), GroupState::Free),
        ];
        assert_eq!(previous_state, Some(GroupState::Free));
        assert_targeted_groups!(c_terminal(None), c_terminal_groups);

        let sidechain = Target::new("Amino", Some("Sidechain"), Some("Lysine"));
        let residue_found = state.update_group_index(sidechain, ResidueId(2), GroupState::Free);
        assert!(residue_found.is_some());

        let residue_not_found = state.update_group_index(sidechain, ResidueId(1), GroupState::Free);
        assert!(residue_not_found.is_none());

        let crazy = Target::new("Silly", Some("Goofy"), None);
        let group_not_found = state.update_group_index(crazy, ResidueId(1), GroupState::Free);
        assert!(group_not_found.is_none());

        let group_still_not_found = state.update_group_index(crazy, ResidueId(2), GroupState::Free);
        assert!(group_still_not_found.is_none());
    }

    #[test]
    fn polymer_groups() {
        let mut state = PolymerizerState::new(&POLYMERIZER);
        let alanine = Residue::new(&POLYMER_DB, "A").unwrap();
        let lysine = Residue::new(&POLYMER_DB, "K").unwrap();
        state.index_residue_groups(ResidueId(0), &alanine);
        state.index_residue_groups(ResidueId(1), &alanine);
        state.index_residue_groups(ResidueId(2), &lysine);

        macro_rules! assert_polymer_groups {
            ($targets:expr, $output:expr) => {
                let mut groups: Vec<_> = state
                    .polymer_groups($targets)
                    .map(|(fg, m)| {
                        let mut id_states: Vec<_> = m.iter().map(|(&i, &f)| (i, f)).collect();
                        id_states.sort_unstable();
                        (fg, id_states)
                    })
                    .collect();
                groups.sort_unstable();
                assert_eq!(groups, $output)
            };
        }
        let amino = Target::new("Amino", None, None);
        let amino_groups = vec![
            (
                FunctionalGroup::new("Amino", "N-Terminal"),
                vec![
                    (ResidueId(0), GroupState::Free),
                    (ResidueId(1), GroupState::Free),
                    (ResidueId(2), GroupState::Free),
                ],
            ),
            (
                FunctionalGroup::new("Amino", "Sidechain"),
                vec![(ResidueId(2), GroupState::Free)],
            ),
        ];
        assert_polymer_groups!(&[amino], amino_groups);

        let carboxyl = Target::new("Carboxyl", None, None);
        let carboxyl_groups = vec![(
            FunctionalGroup::new("Carboxyl", "C-Terminal"),
            vec![
                (ResidueId(0), GroupState::Free),
                (ResidueId(1), GroupState::Free),
                (ResidueId(2), GroupState::Free),
            ],
        )];
        assert_polymer_groups!(&[carboxyl], carboxyl_groups);

        // Looking for both targets merges the results
        let both_groups = [amino_groups, carboxyl_groups].concat();
        assert_polymer_groups!(&[amino, carboxyl], both_groups);

        // Adding in an overlapping (subset) target doesn't change the result
        let n_terminal = Target::new("Amino", Some("N-Terminal"), None);
        assert_polymer_groups!(&[amino, carboxyl, n_terminal], both_groups);

        // Looking for a particular type of residue restricts the results
        let alanine_n_terminal = Target::new("Amino", Some("N-Terminal"), Some("Alanine"));
        let alanine_n_terminal_groups = vec![(
            FunctionalGroup::new("Amino", "N-Terminal"),
            vec![
                (ResidueId(0), GroupState::Free),
                (ResidueId(1), GroupState::Free),
            ],
        )];
        assert_polymer_groups!(&[alanine_n_terminal], alanine_n_terminal_groups);
    }

    #[test]
    fn free_polymer_groups() {
        let mut state = PolymerizerState::new(&POLYMERIZER);
        let mut alanine = Residue::new(&POLYMER_DB, "A").unwrap();
        let mut lysine = Residue::new(&POLYMER_DB, "K").unwrap();
        // Residue 0 is indexed before so that all of its groups remain free!
        state.index_residue_groups(ResidueId(0), &alanine);

        // Manually modify a couple of groups so that they are filtered out later
        let c_terminal_carboxyl = FunctionalGroup::new("Carboxyl", "C-Terminal");
        *alanine.group_state_mut(&c_terminal_carboxyl).unwrap() =
            GroupState::Modified(ModificationId(42));
        let sidechain_amino = FunctionalGroup::new("Amino", "Sidechain");
        *lysine.group_state_mut(&sidechain_amino).unwrap() =
            GroupState::Modified(ModificationId(9001));

        state.index_residue_groups(ResidueId(1), &alanine);
        state.index_residue_groups(ResidueId(2), &lysine);

        macro_rules! assert_free_polymer_groups {
            ($targets:expr, $output:expr) => {
                let mut groups: Vec<_> = state.free_polymer_groups($targets).collect();
                groups.sort_unstable();
                assert_eq!(groups, $output)
            };
        }

        let amino = Target::new("Amino", None, None);
        let amino_groups = vec![
            ResidueGroup(ResidueId(0), FunctionalGroup::new("Amino", "N-Terminal")),
            ResidueGroup(ResidueId(1), FunctionalGroup::new("Amino", "N-Terminal")),
            ResidueGroup(ResidueId(2), FunctionalGroup::new("Amino", "N-Terminal")),
        ];
        assert_free_polymer_groups!(&[amino], amino_groups);

        let carboxyl = Target::new("Carboxyl", None, None);
        let carboxyl_groups = vec![
            ResidueGroup(ResidueId(0), FunctionalGroup::new("Carboxyl", "C-Terminal")),
            ResidueGroup(ResidueId(2), FunctionalGroup::new("Carboxyl", "C-Terminal")),
        ];
        assert_free_polymer_groups!(&[carboxyl], carboxyl_groups);

        // Looking for both targets merges the results
        let both_groups = [
            amino_groups[0],
            carboxyl_groups[0],
            amino_groups[1],
            amino_groups[2],
            carboxyl_groups[1],
        ];
        assert_free_polymer_groups!(&[amino, carboxyl], both_groups);

        // Adding in an overlapping (subset) target doesn't change the result
        let n_terminal = Target::new("Amino", Some("N-Terminal"), None);
        assert_free_polymer_groups!(&[amino, carboxyl, n_terminal], both_groups);

        // Looking for a particular type of residue restricts the results
        let alanine_c_terminal = Target::new("Carboxyl", Some("C-Terminal"), Some("Alanine"));
        let alanine_c_terminal_groups = vec![ResidueGroup(
            ResidueId(0),
            FunctionalGroup::new("Carboxyl", "C-Terminal"),
        )];
        assert_free_polymer_groups!(&[alanine_c_terminal], alanine_c_terminal_groups);
    }

    #[test]
    fn residue_groups() {
        let mut state = PolymerizerState::new(&POLYMERIZER);
        let alanine = Residue::new(&POLYMER_DB, "A").unwrap();
        let lysine = Residue::new(&POLYMER_DB, "K").unwrap();
        state.index_residue_groups(ResidueId(0), &alanine);
        state.index_residue_groups(ResidueId(1), &lysine);

        macro_rules! assert_residue_groups {
            ($targets:expr, $residue:expr, $output:expr) => {
                let mut groups: Vec<_> = state.residue_groups($targets, $residue).collect();
                groups.sort_unstable();
                assert_eq!(groups, $output)
            };
        }
        let amino = Target::new("Amino", None, None);
        let ala_amino_groups = vec![(
            FunctionalGroup::new("Amino", "N-Terminal"),
            GroupState::Free,
        )];
        assert_residue_groups!(&[amino], ResidueId(0), ala_amino_groups);

        let lys_amino_groups = vec![
            (
                FunctionalGroup::new("Amino", "N-Terminal"),
                GroupState::Free,
            ),
            (FunctionalGroup::new("Amino", "Sidechain"), GroupState::Free),
        ];
        assert_residue_groups!(&[amino], ResidueId(1), lys_amino_groups);

        // A more specific target will only show the N-Terminal Amino group (regardless of Ala vs Lys)
        let n_terminal = Target::new("Amino", Some("N-Terminal"), None);
        assert_residue_groups!(&[n_terminal], ResidueId(0), ala_amino_groups);
        assert_residue_groups!(&[n_terminal], ResidueId(1), ala_amino_groups);

        // Carboxyl groups are only present at the C-Terminal
        let carboxyl = Target::new("Carboxyl", None, None);
        let carboxyl_groups = vec![(
            FunctionalGroup::new("Carboxyl", "C-Terminal"),
            GroupState::Free,
        )];
        assert_residue_groups!(&[carboxyl], ResidueId(0), carboxyl_groups);
        assert_residue_groups!(&[carboxyl], ResidueId(1), carboxyl_groups);

        // Looking for both Amino and Carboxyl targets merges the results
        let ala_both_groups = [ala_amino_groups, carboxyl_groups.clone()].concat();
        let lys_both_groups = [lys_amino_groups, carboxyl_groups].concat();
        assert_residue_groups!(&[amino, carboxyl], ResidueId(0), ala_both_groups);
        assert_residue_groups!(&[amino, carboxyl], ResidueId(1), lys_both_groups);
    }

    #[test]
    fn free_residue_groups() {
        let mut state = PolymerizerState::new(&POLYMERIZER);
        let mut alanine = Residue::new(&POLYMER_DB, "A").unwrap();
        let mut lysine = Residue::new(&POLYMER_DB, "K").unwrap();

        // Manually modify a couple of groups so that they are filtered out later
        let c_terminal_carboxyl = FunctionalGroup::new("Carboxyl", "C-Terminal");
        *alanine.group_state_mut(&c_terminal_carboxyl).unwrap() =
            GroupState::Modified(ModificationId(42));
        let sidechain_amino = FunctionalGroup::new("Amino", "Sidechain");
        *lysine.group_state_mut(&sidechain_amino).unwrap() =
            GroupState::Modified(ModificationId(9001));

        state.index_residue_groups(ResidueId(0), &alanine);
        state.index_residue_groups(ResidueId(1), &lysine);

        macro_rules! assert_free_residue_groups {
            ($targets:expr, $residue:expr, $output:expr) => {
                let mut groups: Vec<_> = state.free_residue_groups($targets, $residue).collect();
                groups.sort_unstable();
                assert_eq!(groups, $output)
            };
        }
        let amino = Target::new("Amino", None, None);
        let amino_groups = vec![FunctionalGroup::new("Amino", "N-Terminal")];
        assert_free_residue_groups!(&[amino], ResidueId(0), amino_groups);

        // Because the sidechain group is modified, it isn't returned, and the result looks the same as that for alanine
        assert_free_residue_groups!(&[amino], ResidueId(1), amino_groups);

        // A more specific target will only show the N-Terminal Amino group (regardless of Ala vs Lys)
        let n_terminal = Target::new("Amino", Some("N-Terminal"), None);
        assert_free_residue_groups!(&[n_terminal], ResidueId(0), amino_groups);
        assert_free_residue_groups!(&[n_terminal], ResidueId(1), amino_groups);

        // Carboxyl groups are only present at the C-Terminal, and is already modified for alanine!
        let carboxyl = Target::new("Carboxyl", None, None);
        let carboxyl_groups = vec![FunctionalGroup::new("Carboxyl", "C-Terminal")];
        assert_free_residue_groups!(&[carboxyl], ResidueId(0), Vec::new());
        assert_free_residue_groups!(&[carboxyl], ResidueId(1), carboxyl_groups);

        // Looking for both Amino and Carboxyl targets merges the results
        let both_groups = [amino_groups.clone(), carboxyl_groups].concat();
        assert_free_residue_groups!(&[amino, carboxyl], ResidueId(0), amino_groups);
        assert_free_residue_groups!(&[amino, carboxyl], ResidueId(1), both_groups);
    }

    #[test]
    fn find_single_unambiguous_free_groups() {
        let mut state = PolymerizerState::new(&POLYMERIZER);
        let mut murnac = Residue::new(&POLYMER_DB, "m").unwrap();
        state.index_residue_groups(ResidueId(0), &murnac);

        let carboxyl = Target::new("Carboxyl", None, None);
        let carboxyl_group = state
            .find_any_free_groups(&[carboxyl], ResidueId(0), 1)
            .unwrap();
        assert_eq!(carboxyl_group.len(), 1);
        assert!(carboxyl_group.contains(&FunctionalGroup::new("Carboxyl", "Lactyl Ether")));

        let zero_groups = state.find_any_free_groups(&[carboxyl], ResidueId(0), 0);
        assert_miette_snapshot!(zero_groups);

        let hydroxyl = Target::new("Hydroxyl", None, None);
        let ambiguous_group = state.find_any_free_groups(&[hydroxyl], ResidueId(0), 1);
        assert_miette_snapshot!(ambiguous_group);

        let amino = Target::new("Amino", None, None);
        let crazy = Target::new("Crazy", None, None);
        let no_matching_groups = state.find_any_free_groups(&[amino, crazy], ResidueId(0), 1);
        assert_miette_snapshot!(no_matching_groups);

        // Modify all of MurNAc's "Hydroxyl" groups and try to find a free one again
        let mut hydroxyl_groups: Vec<_> = state
            .residue_groups(&[hydroxyl], ResidueId(0))
            .map(|(fg, _)| fg)
            .collect();
        hydroxyl_groups.sort_unstable();
        let group_states = [
            GroupState::Modified(ModificationId(1)),
            GroupState::Donor(BondId(2)),
            GroupState::Acceptor(BondId(3)),
        ];
        for (fg, gs) in zip(hydroxyl_groups, group_states) {
            *murnac.group_state_mut(&fg).unwrap() = gs;
        }
        state.index_residue_groups(ResidueId(0), &murnac);
        let all_groups_occupied = state.find_any_free_groups(&[hydroxyl], ResidueId(0), 1);
        assert_miette_snapshot!(all_groups_occupied);
    }

    #[test]
    fn find_multiple_unambiguous_free_groups() {
        let mut state = PolymerizerState::new(&POLYMERIZER);
        let mut murnac = Residue::new(&POLYMER_DB, "m").unwrap();
        state.index_residue_groups(ResidueId(0), &murnac);

        let hydroxyl = Target::new("Hydroxyl", None, None);
        let hydroxyl_groups = state
            .find_any_free_groups(&[hydroxyl], ResidueId(0), 3)
            .unwrap();
        assert_eq!(hydroxyl_groups.len(), 3);
        for location in ["Reducing End", "Nonreducing End", "6-Position"] {
            assert!(hydroxyl_groups.contains(&FunctionalGroup::new("Hydroxyl", location)));
        }

        let too_few_hydroxyl_groups = state.find_any_free_groups(&[hydroxyl], ResidueId(0), 4);
        assert_miette_snapshot!(too_few_hydroxyl_groups);

        let ambiguous_group = state.find_any_free_groups(&[hydroxyl], ResidueId(0), 2);
        assert_miette_snapshot!(ambiguous_group);

        let nonreducing_end = FunctionalGroup::new("Hydroxyl", "Nonreducing End");
        *murnac.group_state_mut(&nonreducing_end).unwrap() =
            GroupState::Modified(ModificationId(42));
        state.index_residue_groups(ResidueId(0), &murnac);
        let groups_occupied = state.find_any_free_groups(&[hydroxyl], ResidueId(0), 3);
        assert_miette_snapshot!(groups_occupied);

        let remaining_hydroxyl_groups = state
            .find_any_free_groups(&[hydroxyl], ResidueId(0), 2)
            .unwrap();
        assert_eq!(remaining_hydroxyl_groups.len(), 2);
        for location in ["Reducing End", "6-Position"] {
            assert!(
                remaining_hydroxyl_groups.contains(&FunctionalGroup::new("Hydroxyl", location))
            );
        }

        let lysine = Residue::new(&POLYMER_DB, "K").unwrap();
        state.index_residue_groups(ResidueId(1), &lysine);
        let n_terminal = Target::new("Amino", Some("N-Terminal"), None);
        let carboxyl = Target::new("Carboxyl", None, None);
        let terminals = state
            .find_any_free_groups(&[n_terminal, carboxyl], ResidueId(1), 2)
            .unwrap();
        assert_eq!(terminals.len(), 2);
        assert!(terminals.contains(&FunctionalGroup::new("Amino", "N-Terminal")));
        assert!(terminals.contains(&FunctionalGroup::new("Carboxyl", "C-Terminal")));
    }

    // NOTE: Needed to help type-inference figure out what `[]` is supposed to be
    const NO_GROUPS: [FunctionalGroup; 0] = [];

    #[test]
    fn find_single_specific_free_groups() {
        let mut state = PolymerizerState::new(&POLYMERIZER);
        let mut murnac = Residue::new(&POLYMER_DB, "m").unwrap();
        state.index_residue_groups(ResidueId(0), &murnac);

        let hydroxyl = Target::new("Hydroxyl", None, None);
        let nonreducing_end = FunctionalGroup::new("Hydroxyl", "Nonreducing End");
        let hydroxyl_group = state
            .find_these_free_groups(&[hydroxyl], ResidueId(0), [nonreducing_end])
            .unwrap();
        assert_eq!(hydroxyl_group.len(), 1);
        assert!(hydroxyl_group.contains(&FunctionalGroup::new("Hydroxyl", "Nonreducing End")));

        let no_target_groups = state.find_these_free_groups(&[hydroxyl], ResidueId(0), NO_GROUPS);
        assert_miette_snapshot!(no_target_groups);

        let duplicate_target_group = state.find_these_free_groups(
            &[hydroxyl],
            ResidueId(0),
            [nonreducing_end, nonreducing_end],
        );
        assert_miette_snapshot!(duplicate_target_group);

        let crazy = Target::new("Crazy", None, None);
        let invalid_target =
            state.find_these_free_groups(&[crazy], ResidueId(0), [nonreducing_end]);
        assert_miette_snapshot!(invalid_target);

        let alanine = Residue::new(&POLYMER_DB, "A").unwrap();
        state.index_residue_groups(ResidueId(1), &alanine);
        let nonexistent_group =
            state.find_these_free_groups(&[hydroxyl], ResidueId(1), [nonreducing_end]);
        assert_miette_snapshot!(nonexistent_group);

        *murnac.group_state_mut(&nonreducing_end).unwrap() =
            GroupState::Modified(ModificationId(42));
        state.index_residue_groups(ResidueId(0), &murnac);
        let group_occupied =
            state.find_these_free_groups(&[hydroxyl], ResidueId(0), [nonreducing_end]);
        assert_miette_snapshot!(group_occupied);
    }

    #[test]
    fn find_multiple_specific_free_groups() {
        let mut state = PolymerizerState::new(&POLYMERIZER);
        let mut murnac = Residue::new(&POLYMER_DB, "m").unwrap();
        state.index_residue_groups(ResidueId(0), &murnac);

        let hydroxyl = Target::new("Hydroxyl", None, None);
        let reducing_end = FunctionalGroup::new("Hydroxyl", "Reducing End");
        let nonreducing_end = FunctionalGroup::new("Hydroxyl", "Nonreducing End");
        let six_position = FunctionalGroup::new("Hydroxyl", "6-Position");

        let one_hydroxyl_group = state
            .find_these_free_groups(&[hydroxyl], ResidueId(0), [reducing_end])
            .unwrap();
        assert_eq!(one_hydroxyl_group.len(), 1);
        assert!(one_hydroxyl_group.contains(&reducing_end));

        let two_hydroxyl_groups = state
            .find_these_free_groups(&[hydroxyl], ResidueId(0), [reducing_end, nonreducing_end])
            .unwrap();
        assert_eq!(two_hydroxyl_groups.len(), 2);
        assert!(two_hydroxyl_groups.contains(&reducing_end));
        assert!(two_hydroxyl_groups.contains(&nonreducing_end));

        let all_hydroxyl_groups = state
            .find_these_free_groups(
                &[hydroxyl],
                ResidueId(0),
                [reducing_end, nonreducing_end, six_position],
            )
            .unwrap();
        assert_eq!(all_hydroxyl_groups.len(), 3);
        assert!(all_hydroxyl_groups.contains(&reducing_end));
        assert!(all_hydroxyl_groups.contains(&nonreducing_end));
        assert!(all_hydroxyl_groups.contains(&six_position));

        let n_terminal = FunctionalGroup::new("Amino", "N-Terminal");
        let invalid_target = state.find_these_free_groups(
            &[hydroxyl],
            ResidueId(0),
            [reducing_end, nonreducing_end, six_position, n_terminal],
        );
        assert_miette_snapshot!(invalid_target);

        let crazy = Target::new("Crazy", None, None);
        let still_invalid_target = state.find_these_free_groups(
            &[hydroxyl, crazy],
            ResidueId(0),
            [reducing_end, nonreducing_end, six_position, n_terminal],
        );
        assert_miette_snapshot!(still_invalid_target);

        let amino = Target::new("Amino", None, None);
        let nonexistent_group = state.find_these_free_groups(
            &[hydroxyl, crazy, amino],
            ResidueId(0),
            [reducing_end, nonreducing_end, six_position, n_terminal],
        );
        assert_miette_snapshot!(nonexistent_group);

        *murnac.group_state_mut(&nonreducing_end).unwrap() =
            GroupState::Modified(ModificationId(42));
        state.index_residue_groups(ResidueId(0), &murnac);
        let multiple_errors = state.find_these_free_groups(
            &[hydroxyl, crazy, amino],
            ResidueId(0),
            [reducing_end, nonreducing_end, six_position, n_terminal],
        );
        assert_miette_snapshot!(multiple_errors);
    }
}
