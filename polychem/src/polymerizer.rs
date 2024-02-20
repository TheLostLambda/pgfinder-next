use crate::{atoms::atomic_database::AtomicDatabase, polymers::polymer_database::PolymerDatabase};

struct Polymerizer {
    atomic_db: AtomicDatabase,
    polymer_db: PolymerDatabase,
    residue_idx: usize,
}

impl Polymerizer {
    pub const fn new(atomic_db: AtomicDatabase, polymer_db: PolymerDatabase) -> Self {
        Self {
            atomic_db,
            polymer_db,
            residue_idx: 0,
        }
    }
}
