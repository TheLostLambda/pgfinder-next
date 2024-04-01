use crate::{BondTarget, FunctionalGroup, Id};

impl<'p> BondTarget<'p> {
    #[must_use]
    pub const fn new(residue: Id, group: FunctionalGroup<'p>) -> Self {
        Self { residue, group }
    }
}
