use rust_decimal::Decimal;

use crate::{GroupState, Residue, Result};

use super::polymer_chemistry::{PolymerChemistry, ResidueDescription};

impl Residue {
    // FIXME: Get rid of unwrap
    pub fn new(db: &PolymerChemistry, abbr: impl AsRef<str>, id: usize) -> Self {
        let abbr = abbr.as_ref();
        let ResidueDescription {
            name,
            composition,
            functional_groups,
        } = db.residues.get(abbr).unwrap();
        let functional_groups = functional_groups
            .iter()
            // FIXME: Another clone I shouldn't need...
            .map(|fg| (fg.clone(), GroupState::default()))
            .collect();
        // FIXME: Need to eradicate these clones
        Self {
            id,
            abbr: abbr.to_owned(),
            name: name.clone(),
            composition: composition.clone(),
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
        atoms::atomic_database::AtomicDatabase, polymers::polymer_chemistry::PolymerChemistry,
        Residue,
    };

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(|| {
        AtomicDatabase::from_kdl(
            "atomic_database.kdl",
            include_str!("../../atomic_database.kdl"),
        )
        .unwrap()
    });

    static POLYMER_DB: Lazy<PolymerChemistry> = Lazy::new(|| {
        PolymerChemistry::from_kdl(
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
