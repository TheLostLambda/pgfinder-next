use ahash::HashMap;

use crate::{moieties::target::Index, AtomicDatabase, Id, PolymerDatabase, ResidueId};

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct Polymerizer<'a, 'p> {
    atomic_db: &'a AtomicDatabase,
    polymer_db: &'p PolymerDatabase<'a>,
}

#[derive(Clone, Eq, PartialEq, Debug)]
pub(crate) struct PolymerizerState<'a, 'p> {
    atomic_db: &'a AtomicDatabase,
    polymer_db: &'p PolymerDatabase<'a>,
    next_id: Id,
    free_groups: Index<'p, HashMap<ResidueId, bool>>,
}
