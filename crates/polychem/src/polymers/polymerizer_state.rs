use std::{borrow::Cow, collections::hash_map::Entry};

use ahash::{HashMap, HashMapExt};

use crate::{
    moieties::target::{Index, Target},
    FunctionalGroup, GroupState, Id, Residue, ResidueId,
};

use super::polymerizer::Polymerizer;

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
    pub fn next_id(&mut self) -> Id {
        let current_id = self.next_id;
        self.next_id += 1;
        current_id
    }

    pub fn index_residue_groups(&mut self, id: ResidueId, residue: &Residue<'a, 'p>) {
        for (group, state) in residue.functional_groups() {
            let target = Target::from_residue_and_group(residue, &group);
            self.group_index
                .entry(target)
                .or_default()
                .insert(id, state);
        }
    }

    pub fn polymer_groups<T: Into<Target<&'p str>>>(
        &self,
        targets: impl IntoIterator<Item = T>,
    ) -> HashMap<FunctionalGroup<'p>, Cow<'_, HashMap<ResidueId, GroupState>>> {
        let polymer_groups = targets.into_iter().flat_map(|target| {
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
        deduplicated_groups
    }

    pub fn residue_groups<T: Into<Target<&'p str>>>(
        &self,
        targets: impl IntoIterator<Item = T>,
        residue: ResidueId,
    ) -> impl Iterator<Item = (FunctionalGroup<'p>, GroupState)> + '_ {
        self.polymer_groups(targets)
            .into_iter()
            .filter_map(move |(fg, ids)| ids.get(&residue).map(|&gs| (fg, gs)))
    }
}

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;

    use crate::{moieties::target::Target, AtomicDatabase, PolymerDatabase, Residue};

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

    static POLYMERIZER: Lazy<Polymerizer> = Lazy::new(|| Polymerizer::new(&ATOMIC_DB, &POLYMER_DB));

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
                    .into_iter()
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
}
