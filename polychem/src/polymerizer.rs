use crate::{
    atoms::atomic_database::AtomicDatabase, polymers::polymer_database::PolymerDatabase, Residue,
    Result,
};

// FIXME: This API in general will need some work... Is it better to take a _guard() like `tracing` does and to let
// the scope of that value determine the range in which numbering is unique? I think it might be... But for now,
// the `Polymerizer` is Copy, so I could just create a new one each time... I think, however, once I get to
// building cache-mappings, I won't want to just blow this whole struct away...
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
    // FIXME: Maybe just call this function `residue`?
    pub fn new_residue(&mut self, abbr: impl AsRef<str>) -> Result<Residue> {
        self.residue_idx += 1;
        Residue::new(self.polymer_db, abbr, self.residue_idx)
    }
}
