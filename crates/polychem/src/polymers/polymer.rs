use crate::{AtomicDatabase, Polymer, PolymerDatabase, Residue, ResidueId, Result};

impl<'a, 'p> Polymer<'a, 'p> {
    #[must_use]
    pub const fn atomic_db(&self) -> &'a AtomicDatabase {
        self.polymerizer_state.polymerizer.atomic_db()
    }

    #[must_use]
    pub const fn polymer_db(&self) -> &'p PolymerDatabase<'a> {
        self.polymerizer_state.polymerizer.polymer_db()
    }

    pub fn new_residue(&mut self, abbr: impl AsRef<str>) -> Result<ResidueId> {
        let residue = Residue::new(self.polymer_db(), abbr)?;
        let id = ResidueId(self.polymerizer_state.next_id());
        self.polymerizer_state.index_residue_groups(id, &residue);
        self.residues.insert(id, residue);
        Ok(id)
    }

    #[must_use]
    pub fn residue(&self, id: ResidueId) -> Option<&Residue<'a, 'p>> {
        self.residues.get(&id)
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_ron_snapshot;
    use once_cell::sync::Lazy;

    use crate::{
        polymers::polymerizer::Polymerizer, testing_tools::assert_miette_snapshot, ResidueId,
    };

    use super::*;

    const STEM_RESIDUES: [&str; 4] = ["A", "E", "J", "A"];

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
    fn recover_databases() {
        let polymer = POLYMERIZER.new_polymer();
        assert_eq!(polymer.atomic_db(), &*ATOMIC_DB);
        assert_eq!(polymer.polymer_db(), &*POLYMER_DB);
    }

    #[test]
    fn new_residue() {
        let mut polymer = POLYMERIZER.new_polymer();
        let residues = STEM_RESIDUES.map(|abbr| polymer.new_residue(abbr).unwrap());
        assert_eq!(
            residues,
            [ResidueId(0), ResidueId(1), ResidueId(2), ResidueId(3)]
        );

        let more_residues = STEM_RESIDUES.map(|abbr| polymer.new_residue(abbr).unwrap());
        assert_eq!(
            more_residues,
            [ResidueId(4), ResidueId(5), ResidueId(6), ResidueId(7)]
        );

        let residue_refs = residues.map(|id| polymer.residue(id).unwrap());
        assert_ron_snapshot!(residue_refs, {
            ".**.isotopes, .**.functional_groups" => insta::sorted_redaction()
        });

        let nonexistent_residue = polymer.new_residue("?");
        assert_miette_snapshot!(nonexistent_residue);

        let missing_residue = polymer.residue(ResidueId(8));
        assert_eq!(missing_residue, None);
    }
}
