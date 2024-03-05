use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use miette::Diagnostic;
use thiserror::Error;

use crate::{
    atoms::atomic_database::AtomicDatabase,
    polymers::{
        polymer_database::{BondDescription, ModificationDescription, PolymerDatabase},
        target::{Target, TargetIndex},
    },
    Bond, FunctionalGroup, GroupState, Id, NamedMod, PolychemError, Residue, Result,
};

#[derive(Clone)]
pub struct Polymerizer<'a, 'p> {
    atomic_db: &'a AtomicDatabase,
    polymer_db: &'p PolymerDatabase<'a>,
    residue_counter: usize,
    free_group_index: TargetIndex<'p, HashMap<Id, bool>>,
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
                .insert(residue.id, true);
        }

        Ok(residue)
    }

    // FIXME: Move to bottom of impl block
    fn update_group(
        &mut self,
        valid_targets: &[Target],
        target: &mut Residue<'a, 'p>,
        target_group: &'p FunctionalGroup,
        group_state: GroupState<'a, 'p>,
    ) -> Result<(), PolymerizerError> {
        let current_target = Target::from_residue_and_group(target, target_group);
        let valid_groups: Vec<_> = valid_targets
            .iter()
            .flat_map(|t| self.free_group_index.matches_with_targets(t))
            .collect();

        let target_is_valid = valid_groups.iter().find_map(|&(possible_target, ids)| {
            if current_target == possible_target {
                ids.get(&target.id())
            } else {
                None
            }
        });

        if let Some(target_is_free) = target_is_valid {
            if !target_is_free {
                return Err(PolymerizerError::group_occupied(target_group, target));
            }
        } else {
            let theoretically_possible = valid_targets
                .iter()
                .any(|possible_target| current_target.matches(&possible_target.into()));
            if theoretically_possible {
                return Err(PolymerizerError::residue_not_in_polymer(target));
            }
            return Err(PolymerizerError::invalid_target(
                &current_target,
                valid_targets,
            ));
        }

        self.free_group_index
            .get_mut(current_target)
            .unwrap()
            .insert(target.id(), false);
        // SAFETY: The TargetIndex told us that this group exists and is free, so unwrap() is safe
        let target_state = target.group_state_mut(target_group).unwrap();

        *target_state = group_state;
        Ok(())
    }

    pub fn modify_group(
        &mut self,
        abbr: impl AsRef<str>,
        target: &mut Residue<'a, 'p>,
        target_group: &'p FunctionalGroup,
    ) -> Result<()> {
        let (
            abbr,
            ModificationDescription {
                name,
                lost,
                gained,
                targets,
            },
        ) = NamedMod::lookup_description(self.polymer_db, abbr)?;

        let modified_state = GroupState::Modified(NamedMod {
            abbr,
            name,
            lost,
            gained,
        });

        self.update_group(targets, target, target_group, modified_state)
            .map_err(|e| PolychemError::modification(name, abbr, target, e).into())
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

        // if !donor_state.is_free() {
        //     return Err(PolychemError::bond_group_occupied(donor_group, donor).into());
        // }
        // if !acceptor_state.is_free() {
        //     return Err(PolychemError::bond_group_occupied(acceptor_group, acceptor).into());
        // }

        // let (kind, BondDescription { from, to, lost }) =
        //     Bond::lookup_description(self.polymer_db, kind)?;
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

#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
pub(crate) enum PolymerizerError {
    #[error("the functional group {0} of residue {1} was already {2}, but must be free")]
    GroupOccupied(FunctionalGroup, Id, String),

    #[error("residue {0} does not belong to the current polymer")]
    #[diagnostic(help(
        "the referenced residue was likely created by a different Polymerizer instance"
    ))]
    ResidueNotInPolymer(Id),

    #[error("expected a target matching {0}, got {1}")]
    InvalidTarget(String, Target),
}

impl PolymerizerError {
    fn group_occupied(group: &FunctionalGroup, residue: &Residue) -> Self {
        Self::GroupOccupied(
            group.clone(),
            residue.id,
            residue.group_state(group).unwrap().to_string(),
        )
    }

    fn residue_not_in_polymer(residue: &Residue) -> Self {
        Self::ResidueNotInPolymer(residue.id)
    }

    fn invalid_target(target: impl Into<Target>, valid_targets: &[Target]) -> Self {
        let valid_targets = valid_targets.iter().join(", or");
        Self::InvalidTarget(valid_targets, target.into())
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_ron_snapshot;
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::{
        atoms::atomic_database::AtomicDatabase, polymers::polymer_database::PolymerDatabase,
        testing_tools::assert_miette_snapshot, FunctionalGroup, GroupState, Massive,
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

    #[test]
    fn modify_group() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut murnac = polymerizer.new_residue("m").unwrap();
        assert_eq!(murnac.monoisotopic_mass(), dec!(293.11106657336));

        let reducing_end = FunctionalGroup::new("Hydroxyl", "Reducing End");
        polymerizer
            .modify_group("Anh", &mut murnac, &reducing_end)
            .unwrap();
        assert_eq!(murnac.monoisotopic_mass(), dec!(275.10050188933));
        assert!(matches!(
            murnac.group_state(&reducing_end).unwrap(),
            GroupState::Modified(_)
        ));

        let modify_non_free_group = polymerizer.modify_group("Anh", &mut murnac, &reducing_end);
        assert_miette_snapshot!(modify_non_free_group);

        // Start a new polymer by resetting the polymerizer
        let mut polymerizer = polymerizer.reset();
        let residue_from_wrong_polymer =
            polymerizer.modify_group("Anh", &mut murnac, &reducing_end);
        assert_miette_snapshot!(residue_from_wrong_polymer);

        let invalid_group = polymerizer.modify_group("Ac", &mut murnac, &reducing_end);
        assert_miette_snapshot!(invalid_group);
    }
}
