use crate::{
    atoms::atomic_database::AtomicDatabase, polymers::polymer_chemistry::PolymerChemistry,
};

struct Polymerizer {
    atomic_db: AtomicDatabase,
    chemistry: PolymerChemistry,
    residue_idx: usize,
}

impl Polymerizer {
    pub const fn new(atomic_db: AtomicDatabase, chemistry: PolymerChemistry) -> Self {
        Self {
            atomic_db,
            chemistry,
            residue_idx: 0,
        }
    }
}
