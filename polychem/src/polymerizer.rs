use crate::{atoms::atomic_database::AtomicDatabase, polymers::polymer_database::PolymerDatabase};

struct Polymerizer<'a, 'p> {
    atomic_db: &'a AtomicDatabase,
    polymer_db: &'p PolymerDatabase<'a>,
    residue_idx: usize,
}

impl<'a, 'p> Polymerizer<'a, 'p> {
    pub const fn new(atomic_db: &'a AtomicDatabase, polymer_db: &'p PolymerDatabase<'a>) -> Self {
        Self {
            atomic_db,
            polymer_db,
            residue_idx: 0,
        }
    }
}
