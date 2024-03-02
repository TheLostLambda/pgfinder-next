use std::collections::HashSet;

use crate::{
    atoms::atomic_database::AtomicDatabase,
    polymers::{
        polymer_database::{BondDescription, ModificationDescription, PolymerDatabase},
        target::{Target, TargetIndex},
    },
    Bond, FunctionalGroup, Id, NamedMod, PolychemError, Residue, Result,
};

#[derive(Clone)]
pub struct Polymerizer<'a, 'p> {
    atomic_db: &'a AtomicDatabase,
    polymer_db: &'p PolymerDatabase<'a>,
    residue_counter: usize,
    free_group_index: TargetIndex<'p, HashSet<Id>>,
}

impl<'a, 'p> Polymerizer<'a, 'p> {
    #[must_use]
    pub fn new(atomic_db: &'a AtomicDatabase, polymer_db: &'p PolymerDatabase<'a>) -> Self {
        Self {
            atomic_db,
            polymer_db,
            residue_counter: 0,
            free_group_index: TargetIndex::new(),
        }
    }

    // FIXME: Should this just be reset(&mut self)?
    #[must_use]
    pub fn reset(self) -> Self {
        Self::new(self.atomic_db, self.polymer_db)
    }

    pub fn new_residue(&mut self, abbr: impl AsRef<str>) -> Result<Residue<'a, 'p>> {
        self.residue_counter += 1;
        let residue = Residue::new(self.polymer_db, abbr, self.residue_counter)?;

        let free_groups = residue
            .functional_groups
            .iter()
            .filter_map(|(&fg, gs)| gs.is_free().then_some(fg));

        for group in free_groups {
            let target = Target::from_residue_and_group(&residue, group);
            self.free_group_index
                .entry(target)
                .or_default()
                .insert(residue.id);
        }

        Ok(residue)
    }

    // FIXME: TODO: Add update_group method that handles adding and removing residue ids to the free-groups index and
    // is also in charge of the actual mutation. This avoids me needing to update that free-groups index everywhere!
    pub fn modify_group(
        &self,
        abbr: impl AsRef<str>,
        target: &mut Residue<'a, 'p>,
        target_group: &FunctionalGroup,
    ) -> Result<()> {
        let target_state = target.group_state_mut(target_group)?;
        let (
            abbr,
            ModificationDescription {
                name,
                lost,
                gained,
                targets,
            },
        ) = NamedMod::lookup_description(self.polymer_db, abbr)?;

        Ok(())
    }

    // FIXME: Finish this!
    pub fn bond_groups(
        &self,
        kind: impl AsRef<str>,
        donor: &mut Residue<'a, 'p>,
        donor_group: &FunctionalGroup,
        acceptor: &mut Residue<'a, 'p>,
        acceptor_group: &FunctionalGroup,
    ) -> Result<()> {
        let donor_state = donor.group_state_mut(donor_group)?;
        let acceptor_state = acceptor.group_state_mut(acceptor_group)?;

        if !donor_state.is_free() {
            return Err(PolychemError::bond_group_occupied(donor_group, donor).into());
        }
        if !acceptor_state.is_free() {
            return Err(PolychemError::bond_group_occupied(acceptor_group, acceptor).into());
        }

        let (kind, BondDescription { from, to, lost }) =
            Bond::lookup_description(self.polymer_db, kind)?;
        Ok(())
    }
}

// FIXME: Oh boy... Where to I belong... Maybe replace with a From impl for tuples?
impl<'a> Target<&'a str> {
    pub(crate) fn from_residue_and_group(
        residue: &Residue<'_, 'a>,
        group: &'a FunctionalGroup,
    ) -> Self {
        Self::new(
            group.name.as_str(),
            Some(group.location.as_str()),
            Some(residue.name),
        )
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
