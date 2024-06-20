use ahash::{HashMap, HashMapExt};

use crate::{AtomicDatabase, Polymer, PolymerDatabase};

use super::polymerizer_state::PolymerizerState;

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct Polymerizer<'a, 'p> {
    atomic_db: &'a AtomicDatabase,
    polymer_db: &'p PolymerDatabase<'a>,
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
    pub const fn atomic_db(&self) -> &'a AtomicDatabase {
        self.atomic_db
    }

    #[must_use]
    pub const fn polymer_db(&self) -> &'p PolymerDatabase<'a> {
        self.polymer_db
    }

    #[must_use]
    pub fn new_polymer(&self) -> Polymer<'a, 'p> {
        let polymerizer_state = PolymerizerState::new(self);
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
    fn recover_databases() {
        let polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        assert_eq!(polymerizer.atomic_db(), &*ATOMIC_DB);
        assert_eq!(polymerizer.polymer_db(), &*POLYMER_DB);
    }

    #[test]
    fn new_polymer() {
        let polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let polymer: Polymer = polymerizer.new_polymer();
        assert_eq!(polymer.polymerizer_state.polymerizer, polymerizer);
        assert_eq!(
            polymer.polymerizer_state,
            PolymerizerState::new(&polymerizer)
        );
    }
}
