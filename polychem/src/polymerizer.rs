use crate::{
    atoms::atomic_database::AtomicDatabase, polymers::polymer_database::PolymerDatabase, Residue,
    Result,
};

// FIXME: This API in general will need some work... Is it better to take a _guard() like `tracing` does and to let
// the scope of that value determine the range in which numbering is unique? I think it might be... But for now,
// the `Polymerizer` is Copy, so I could just create a new one each time... I think, however, once I get to
// building cache-mappings, I won't want to just blow this whole struct away...
#[derive(Copy, Clone)]
pub struct Polymerizer<'a, 'p> {
    atomic_db: &'a AtomicDatabase,
    polymer_db: &'p PolymerDatabase<'a>,
    residue_counter: usize,
}

impl<'a, 'p> Polymerizer<'a, 'p> {
    #[must_use]
    pub const fn new(atomic_db: &'a AtomicDatabase, polymer_db: &'p PolymerDatabase<'a>) -> Self {
        Self {
            atomic_db,
            polymer_db,
            residue_counter: 0,
        }
    }

    // FIXME: Should this just be reset(&mut self)?
    #[must_use]
    pub const fn reset(self) -> Self {
        Self::new(self.atomic_db, self.polymer_db)
    }

    pub fn new_residue(&mut self, abbr: impl AsRef<str>) -> Result<Residue<'a, 'p>> {
        self.residue_counter += 1;
        Residue::new(self.polymer_db, abbr, self.residue_counter)
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_ron_snapshot;
    use once_cell::sync::Lazy;

    use crate::{
        atoms::atomic_database::AtomicDatabase, polymers::polymer_database::PolymerDatabase,
    };

    use super::Polymerizer;

    const STEM_RESIDUES: [&str; 4] = ["A", "E", "J", "A"];

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(|| {
        AtomicDatabase::from_kdl(
            "atomic_database.kdl",
            include_str!("../atomic_database.kdl"),
        )
        .unwrap()
    });

    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::from_kdl(
            &ATOMIC_DB,
            "muropeptide_chemistry.kdl",
            include_str!("../muropeptide_chemistry.kdl"),
        )
        .unwrap()
    });

    #[test]
    fn residue_construction() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let residues = STEM_RESIDUES.map(|abbr| polymerizer.new_residue(abbr).unwrap());
        assert_ron_snapshot!(residues, {
            ".**.isotopes, .**.functional_groups" => insta::sorted_redaction()
        });

        let residues = STEM_RESIDUES.map(|abbr| polymerizer.new_residue(abbr).unwrap().id());
        assert_eq!(residues, [5, 6, 7, 8]);

        let mut polymerizer = polymerizer.reset();
        let residues = STEM_RESIDUES.map(|abbr| polymerizer.new_residue(abbr).unwrap().id());
        assert_eq!(residues, [1, 2, 3, 4]);
    }
}
