use ahash::{HashMap, HashMapExt};

use crate::{moieties::target::Index, AtomicDatabase, Id, Polymer, PolymerDatabase, ResidueId};

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct Polymerizer<'a, 'p> {
    atomic_db: &'a AtomicDatabase,
    polymer_db: &'p PolymerDatabase<'a>,
}

#[derive(Clone, Eq, PartialEq, Debug)]
pub(crate) struct PolymerizerState<'a, 'p> {
    pub polymerizer: Polymerizer<'a, 'p>,
    pub next_id: Id,
    pub free_groups: Index<'p, HashMap<ResidueId, bool>>,
}

impl<'a, 'p> Polymerizer<'a, 'p> {
    #[must_use]
    pub const fn new(atomic_db: &'a AtomicDatabase, polymer_db: &'p PolymerDatabase<'a>) -> Self {
        Self {
            atomic_db,
            polymer_db,
        }
    }

    #[must_use]
    pub fn new_polymer(&self) -> Polymer {
        let polymerizer_state = PolymerizerState {
            polymerizer: *self,
            next_id: Id::default(),
            free_groups: Index::new(),
        };
        Polymer {
            polymerizer_state,
            residues: HashMap::new(),
            modifications: HashMap::new(),
            bonds: HashMap::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;

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

    #[test]
    fn new_polymer() {
        let polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let polymer: Polymer = polymerizer.new_polymer();
        assert_eq!(polymer.polymerizer_state.polymerizer, polymerizer);
        assert_eq!(polymer.polymerizer_state.next_id, 0);
        assert_eq!(polymer.polymerizer_state.free_groups, Index::new());
    }
}
