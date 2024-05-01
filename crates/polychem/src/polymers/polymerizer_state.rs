use ahash::HashMap;

use crate::{
    moieties::target::{Index, Target},
    GroupState, Id, Residue, ResidueId,
};

use super::polymerizer::Polymerizer;

#[derive(Clone, Eq, PartialEq, Debug)]
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
}
