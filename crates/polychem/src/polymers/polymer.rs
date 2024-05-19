use crate::{
    AtomicDatabase, AverageMass, BondInfo, Charge, Charged, Massive, ModificationInfo,
    MonoisotopicMass, Polymer, PolymerDatabase, Residue, ResidueId, Result,
};

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

// NOTE: Since we'll be calling the `$accessor` method on different structs (`Residue`s, `Modification`s, and `Bond`s),
// `$accessor` can't have a concrete type. You could solve this with trait objects, but that sacrifices performance —
// using a macro here keeps things DRY and overhead-free!
macro_rules! sum_parts {
    ($self:expr, $accessor:path) => {{
        let residues = $self.residues.values().map($accessor);
        let modifications = $self.modifications.values().map(|mod_info| match mod_info {
            ModificationInfo::Named(m, ..) => $accessor(m),
            ModificationInfo::Offset(m, ..) => $accessor(m),
            ModificationInfo::Unlocalized(m) => $accessor(m),
        });
        let bonds = $self.bonds.values().map(|BondInfo(_, b, _)| $accessor(b));

        residues.chain(modifications).chain(bonds).sum()
    }};
}

impl Massive for Polymer<'_, '_> {
    fn monoisotopic_mass(&self) -> MonoisotopicMass {
        sum_parts!(self, Massive::monoisotopic_mass)
    }

    fn average_mass(&self) -> AverageMass {
        sum_parts!(self, Massive::average_mass)
    }
}

impl Charged for Polymer<'_, '_> {
    fn charge(&self) -> Charge {
        sum_parts!(self, Charged::charge)
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_ron_snapshot;
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::{
        polymers::polymerizer::Polymerizer, testing_tools::assert_miette_snapshot, ChargedParticle,
        ResidueId,
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

    macro_rules! assert_polymer {
        ($polymer:expr, $mono_mass:literal, $avg_mass:literal) => {
            assert_eq!(
                $polymer.monoisotopic_mass(),
                MonoisotopicMass(dec!($mono_mass))
            );
            assert_eq!($polymer.average_mass(), AverageMass(dec!($avg_mass)));
            assert_eq!($polymer.charge(), Charge(0));
            assert_eq!($polymer.monoisotopic_mz(), None);
            assert_eq!($polymer.average_mz(), None);
        };

        ($polymer:expr, $mono_mass:literal, $avg_mass:literal, $charge:literal, $mono_mz:literal, $avg_mz:literal) => {
            assert_eq!(
                $polymer.monoisotopic_mass(),
                MonoisotopicMass(dec!($mono_mass))
            );
            assert_eq!($polymer.average_mass(), AverageMass(dec!($avg_mass)));
            assert_eq!($polymer.charge(), Charge($charge));
            assert_eq!(
                $polymer.monoisotopic_mz(),
                Some(MonoisotopicMz(dec!($mono_mz)))
            );
            assert_eq!($polymer.average_mz(), Some(AverageMz(dec!($avg_mz))));
        };
    }

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
        assert_polymer!(polymer, 515.24387164950, 515.51342919034656875);

        let more_residues = STEM_RESIDUES.map(|abbr| polymer.new_residue(abbr).unwrap());
        assert_eq!(
            more_residues,
            [ResidueId(4), ResidueId(5), ResidueId(6), ResidueId(7)]
        );
        assert_polymer!(polymer, 1030.48774329900, 1031.02685838069313750);

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
