use rust_decimal::Decimal;

use crate::{GroupState, Residue, Result};

use super::polymer_database::{PolymerDatabase, ResidueDescription};

impl<'a, 'p> Residue<'a, 'p> {
    // FIXME: Get rid of unwrap
    pub fn new(db: &'p PolymerDatabase<'a>, abbr: impl AsRef<str>, id: usize) -> Self {
        let (
            abbr,
            ResidueDescription {
                name,
                composition,
                functional_groups,
            },
        ) = db.residues.get_key_value(abbr.as_ref()).unwrap();
        let functional_groups = functional_groups
            .iter()
            .map(|fg| (fg, GroupState::default()))
            .collect();
        Self {
            id,
            abbr,
            name,
            composition,
            functional_groups,
            offset_modifications: Vec::new(),
        }
    }

    // FIXME: Should these mass functions be made into a trait? I think they probably should be...
    // FIXME: I also don't like using Report anywhere... (with Result<T>) I need to move to module error enums again
    pub fn monoisotopic_mass(&self) -> Result<Decimal> {
        self.composition.monoisotopic_mass()
    }
}

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::{
        atoms::atomic_database::AtomicDatabase, polymers::polymer_database::PolymerDatabase,
        Residue,
    };

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(|| {
        AtomicDatabase::from_kdl(
            "atomic_database.kdl",
            include_str!("../../atomic_database.kdl"),
        )
        .unwrap()
    });

    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::from_kdl(
            &ATOMIC_DB,
            "muropeptide_chemistry.kdl",
            include_str!("../../muropeptide_chemistry.kdl"),
        )
        .unwrap()
    });

    // FIXME: Rubbish, actually split this into sensible tests!
    #[test]
    fn residue_construction() {
        let alanine = Residue::new(&POLYMER_DB, "A", 0);
        assert_eq!(alanine.monoisotopic_mass().unwrap(), dec!(89.04767846918));
    }
}
