use rust_decimal::Decimal;

use crate::{GroupState, Massive, PolychemError, Residue, Result};

use super::polymer_database::{PolymerDatabase, ResidueDescription};

impl<'a, 'p> Residue<'a, 'p> {
    pub fn new(db: &'p PolymerDatabase<'a>, abbr: impl AsRef<str>, id: usize) -> Result<Self> {
        let abbr = abbr.as_ref();
        let (
            abbr,
            ResidueDescription {
                name,
                composition,
                functional_groups,
            },
        ) = db
            .residues
            .get_key_value(abbr)
            .ok_or_else(|| PolychemError::ResidueLookup(abbr.to_owned()))?;
        let functional_groups = functional_groups
            .iter()
            .map(|fg| (fg, GroupState::default()))
            .collect();
        Ok(Self {
            id,
            abbr,
            name,
            composition,
            functional_groups,
            offset_modifications: Vec::new(),
        })
    }
    // TODO: Write named_modifications and bonds
}

// FIXME: These need to take into account offset modifications and functional groups!
impl Massive for Residue<'_, '_> {
    fn monoisotopic_mass(&self) -> Decimal {
        self.composition.monoisotopic_mass()
    }

    fn average_mass(&self) -> Decimal {
        self.composition.average_mass()
    }
}

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::{
        atoms::atomic_database::AtomicDatabase, polymers::polymer_database::PolymerDatabase,
        Massive, Residue,
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
        let alanine = Residue::new(&POLYMER_DB, "A", 0).unwrap();
        assert_eq!(alanine.monoisotopic_mass(), dec!(89.04767846918));
    }
}
