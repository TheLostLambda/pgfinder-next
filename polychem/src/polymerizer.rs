use std::{collections::HashMap, slice};

use itertools::Itertools;
use miette::Diagnostic;
use thiserror::Error;

use crate::{
    atoms::atomic_database::AtomicDatabase,
    polymers::{
        polymer_database::{BondDescription, ModificationDescription, PolymerDatabase},
        target::{Index, Target},
    },
    Bond, BondTarget, FunctionalGroup, GroupState, Id, NamedMod, PolychemError, Residue, Result,
};

#[derive(Clone)]
pub struct Polymerizer<'a, 'p> {
    atomic_db: &'a AtomicDatabase,
    polymer_db: &'p PolymerDatabase<'a>,
    residue_counter: usize,
    free_group_index: Index<'p, HashMap<Id, bool>>,
}

impl<'a, 'p> Polymerizer<'a, 'p> {
    #[must_use]
    pub fn new(atomic_db: &'a AtomicDatabase, polymer_db: &'p PolymerDatabase<'a>) -> Self {
        Self {
            atomic_db,
            polymer_db,
            residue_counter: 0,
            free_group_index: Index::new(),
        }
    }

    #[must_use]
    pub fn reset(self) -> Self {
        Self::new(self.atomic_db, self.polymer_db)
    }

    pub fn new_residue(&mut self, abbr: impl AsRef<str>) -> Result<Residue<'a, 'p>> {
        self.residue_counter += 1;
        let residue = Residue::new(self.polymer_db, abbr, self.residue_counter)?;

        // FIXME: Maybe get rid of this... All of the groups should start free...
        let free_groups = residue
            .functional_groups
            .iter()
            .filter_map(|(&fg, gs)| gs.is_free().then_some(fg));

        // FIXME: And maybe move this to it's own function? But it could be fine here...
        for group in free_groups {
            let target = Target::from_residue_and_group(&residue, group);
            self.free_group_index
                .entry(target)
                .or_default()
                .insert(residue.id(), true);
        }

        Ok(residue)
    }

    pub fn modify(&mut self, abbr: impl AsRef<str>, target: &mut Residue<'a, 'p>) -> Result<()> {
        self.modify_with_optional_group(abbr, target, None)
    }

    // PERF: Could create an `_unchecked` version for when you've already called `self.free_*_groups()` — skip straight
    // to `self.update_group()`!
    pub fn modify_group(
        &mut self,
        abbr: impl AsRef<str>,
        target: &mut Residue<'a, 'p>,
        target_group: FunctionalGroup<'p>,
    ) -> Result<()> {
        self.modify_with_optional_group(abbr, target, Some(target_group))
    }

    fn modify_with_optional_group(
        &mut self,
        abbr: impl AsRef<str>,
        target: &mut Residue<'a, 'p>,
        target_group: Option<FunctionalGroup<'p>>,
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

        let target_group = self
            .find_targeted_group(&targets, target, target_group)
            .map_err(|e| PolychemError::modification(name, abbr, target, e))?;

        let modified_state = GroupState::Modified(NamedMod {
            abbr,
            name,
            lost,
            gained,
        });
        self.update_group(target, target_group, modified_state);

        Ok(())
    }

    pub fn bond(
        &mut self,
        kind: impl AsRef<str>,
        donor: &mut Residue<'a, 'p>,
        acceptor: &mut Residue<'a, 'p>,
    ) -> Result<()> {
        self.bond_with_optional_groups(kind, donor, None, acceptor, None)
    }

    // PERF: Could create an `_unchecked` version for when you've already called `self.free_*_groups()` — skip straight
    // to `self.update_group()`!
    pub fn bond_groups(
        &mut self,
        kind: impl AsRef<str>,
        donor: &mut Residue<'a, 'p>,
        donor_group: FunctionalGroup<'p>,
        acceptor: &mut Residue<'a, 'p>,
        acceptor_group: FunctionalGroup<'p>,
    ) -> Result<()> {
        self.bond_with_optional_groups(
            kind,
            donor,
            Some(donor_group),
            acceptor,
            Some(acceptor_group),
        )
    }

    fn bond_with_optional_groups(
        &mut self,
        kind: impl AsRef<str>,
        donor: &mut Residue<'a, 'p>,
        donor_group: Option<FunctionalGroup<'p>>,
        acceptor: &mut Residue<'a, 'p>,
        acceptor_group: Option<FunctionalGroup<'p>>,
    ) -> Result<()> {
        let (kind, BondDescription { from, to, lost }) =
            Bond::lookup_description(self.polymer_db, kind)?;

        // Avoid partial updates by performing validation of both group updates *before* updating either group
        let donor_group = self
            .find_targeted_group(&slice::from_ref(from), donor, donor_group)
            .map_err(|e| PolychemError::bond(kind, donor, acceptor, "donor", e))?;
        let acceptor_group = self
            .find_targeted_group(&slice::from_ref(to), acceptor, acceptor_group)
            .map_err(|e| PolychemError::bond(kind, donor, acceptor, "acceptor", e))?;

        let donor_state = GroupState::Donor(Bond {
            kind,
            lost,
            acceptor: BondTarget {
                residue: acceptor.id(),
                group: acceptor_group,
            },
        });
        self.update_group(donor, donor_group, donor_state);
        self.update_group(acceptor, acceptor_group, GroupState::Acceptor);

        Ok(())
    }

    fn molecule_groups<'s, T: Into<Target<&'p str>>>(
        &'s self,
        targets: &'s (impl IntoIterator<Item = T> + Copy),
    ) -> impl Iterator<Item = (FunctionalGroup<'p>, &HashMap<Id, bool>)> + 's {
        targets
            .into_iter()
            .flat_map(|target| self.free_group_index.matches_with_targets(target))
            .map(|(target, ids)| {
                (
                    // SAFETY: `.matches_with_targets()` aways returns a complete `Target` with no `None` fields, so
                    // `.unwrap()` is safe
                    FunctionalGroup::new(target.group, target.location.unwrap()),
                    ids,
                )
            })
    }

    // TODO: Write `free_molecule_groups`

    fn residue_groups<'s, T: Into<Target<&'p str>>>(
        &'s self,
        targets: &'s (impl IntoIterator<Item = T> + Copy),
        residue: &Residue<'a, 'p>,
    ) -> impl Iterator<Item = (FunctionalGroup<'p>, bool)> + 's {
        let residue_id = residue.id();
        self.molecule_groups(targets).filter_map(move |(fg, ids)| {
            if let Some(&is_free) = ids.get(&residue_id) {
                Some((fg, is_free))
            } else {
                None
            }
        })
    }

    fn free_residue_groups<'s, T: Into<Target<&'p str>>>(
        &'s self,
        targets: &'s (impl IntoIterator<Item = T> + Copy),
        residue: &Residue<'a, 'p>,
    ) -> impl Iterator<Item = FunctionalGroup<'p>> + 's {
        self.residue_groups(targets, residue)
            .filter_map(|(fg, is_free)| is_free.then_some(fg))
    }

    fn find_targeted_group<T: Into<Target<&'p str>>>(
        &self,
        targets: &(impl IntoIterator<Item = T> + Copy),
        residue: &Residue<'a, 'p>,
        group: Option<FunctionalGroup<'p>>,
    ) -> Result<FunctionalGroup<'p>, PolymerizerError> {
        let free_groups: Vec<_> = self.free_residue_groups(targets, residue).collect();
        match (&free_groups[..], group) {
            ([], None) => {
                let non_free_groups: Vec<_> = self.residue_groups(targets, residue).collect();
                todo!()
            }
            ([], Some(target_group)) => {
                // FIXME: Refactor out into another function!
                let current_target = Target::from_residue_and_group(residue, target_group);

                let group_in_index = self
                    .residue_groups(&[current_target], residue)
                    .next()
                    .is_some();
                let residue_has_group = residue.functional_groups.contains_key(&target_group);
                let theoretically_possible = targets
                    .into_iter()
                    .any(|possible_target| current_target.matches(&possible_target.into()));

                if group_in_index {
                    Err(PolymerizerError::group_occupied(target_group, residue))
                } else if !residue_has_group {
                    Err(PolymerizerError::nonexistent_group(target_group, residue))
                } else if !theoretically_possible {
                    Err(PolymerizerError::invalid_target(&current_target, targets))
                } else {
                    Err(PolymerizerError::residue_not_in_polymer(residue))
                }
            }
            (&[found_group], None) => Ok(found_group),
            (_, None) => panic!("Ambiguous target"),
            (found_groups, Some(target_group)) => {
                if found_groups.contains(&target_group) {
                    Ok(target_group)
                } else {
                    panic!("Shit")
                }
            }
        }
    }

    fn update_group(
        &mut self,
        target: &mut Residue<'a, 'p>,
        target_group: FunctionalGroup<'p>,
        group_state: GroupState<'a, 'p>,
    ) {
        let current_target = Target::from_residue_and_group(target, target_group);

        // FIXME: Update comment
        // SAFETY: This .unwrap() might panic if this update hasn't been pre-validated by self.validate_group_update()!
        self.free_group_index
            .get_mut(current_target)
            .unwrap()
            .insert(target.id(), false);

        // FIXME: Update comment
        // SAFETY: This .unwrap() might panic if this update hasn't been pre-validated by self.validate_group_update()!
        let target_state = target.group_state_mut(&target_group).unwrap();
        *target_state = group_state;
    }
}

// FIXME: Oh boy... Where to I belong... Maybe replace with a From impl for tuples?
impl<'p> Target<&'p str> {
    pub(crate) const fn from_residue_and_group(
        residue: &Residue<'_, 'p>,
        // FIXME: Should I be passing this by reference?!
        group: FunctionalGroup<'p>,
    ) -> Self {
        Self::new(group.name, Some(group.location), Some(residue.name))
    }
}

#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
pub(crate) enum PolymerizerError {
    #[error("the functional group {0} of residue {1} was already {2}, but must be free")]
    GroupOccupied(String, Id, String),

    #[error("the functional group {0} does not exist on residue {1}")]
    NonexistentGroup(String, Id),

    #[error("expected a target matching {0}, got {1}")]
    InvalidTarget(String, Target),

    #[error("residue {0} does not belong to the current polymer")]
    #[diagnostic(help(
        "the referenced residue was likely created by a different Polymerizer instance"
    ))]
    ResidueNotInPolymer(Id),
}

impl PolymerizerError {
    fn group_occupied(group: FunctionalGroup, residue: &Residue) -> Self {
        Self::GroupOccupied(
            group.to_string(),
            residue.id(),
            residue.group_state(&group).unwrap().to_string(),
        )
    }

    fn nonexistent_group(group: FunctionalGroup, residue: &Residue) -> Self {
        Self::NonexistentGroup(group.to_string(), residue.id())
    }

    fn invalid_target<'a, T: Into<Target<&'a str>>>(
        target: impl Into<Target>,
        valid_targets: &(impl IntoIterator<Item = T> + Copy),
    ) -> Self {
        let valid_targets = valid_targets.into_iter().map(Into::into).join(", or");
        Self::InvalidTarget(valid_targets, target.into())
    }

    const fn residue_not_in_polymer(residue: &Residue) -> Self {
        Self::ResidueNotInPolymer(residue.id())
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

        let nonexistent_residue = polymerizer.new_residue("?");
        assert_miette_snapshot!(nonexistent_residue);
    }

    #[test]
    fn modify_group() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut murnac = polymerizer.new_residue("m").unwrap();
        assert_eq!(murnac.monoisotopic_mass(), dec!(293.11106657336));

        let reducing_end = FunctionalGroup::new("Hydroxyl", "Reducing End");
        polymerizer
            .modify_group("Anh", &mut murnac, reducing_end)
            .unwrap();
        assert_eq!(murnac.monoisotopic_mass(), dec!(275.10050188933));
        assert!(matches!(
            murnac.group_state(&reducing_end).unwrap(),
            GroupState::Modified(_)
        ));

        let modify_non_free_group = polymerizer.modify_group("Anh", &mut murnac, reducing_end);
        assert_miette_snapshot!(modify_non_free_group);

        // Start a new polymer by resetting the polymerizer
        let mut polymerizer = polymerizer.reset();
        let residue_from_wrong_polymer = polymerizer.modify_group("Anh", &mut murnac, reducing_end);
        assert_miette_snapshot!(residue_from_wrong_polymer);

        let invalid_group = polymerizer.modify_group("Ac", &mut murnac, reducing_end);
        assert_miette_snapshot!(invalid_group);

        let n_terminal = FunctionalGroup::new("Amino", "N-Terminal");
        let nonexistent_group = polymerizer.modify_group("Ac", &mut murnac, n_terminal);
        assert_miette_snapshot!(nonexistent_group);

        let nonexistent_modification = polymerizer.modify_group("Arg", &mut murnac, reducing_end);
        assert_miette_snapshot!(nonexistent_modification);
    }

    #[test]
    fn bond_groups() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut murnac = polymerizer.new_residue("m").unwrap();
        let mut alanine = polymerizer.new_residue("A").unwrap();
        assert_eq!(
            murnac.monoisotopic_mass() + alanine.monoisotopic_mass(),
            dec!(382.15874504254)
        );

        let lactyl = FunctionalGroup::new("Carboxyl", "Lactyl Ether");
        let n_terminal = FunctionalGroup::new("Amino", "N-Terminal");
        polymerizer
            .bond_groups("Stem", &mut murnac, lactyl, &mut alanine, n_terminal)
            .unwrap();
        assert_eq!(
            murnac.monoisotopic_mass() + alanine.monoisotopic_mass(),
            dec!(364.14818035851)
        );
        assert!(matches!(
            murnac.group_state(&lactyl).unwrap(),
            GroupState::Donor(_)
        ));
        assert!(matches!(
            alanine.group_state(&n_terminal).unwrap(),
            GroupState::Acceptor
        ));

        let groups_not_free =
            polymerizer.bond_groups("Stem", &mut murnac, lactyl, &mut alanine, n_terminal);
        assert_miette_snapshot!(groups_not_free);

        let c_terminal = FunctionalGroup::new("Carboxyl", "C-Terminal");
        let mut glcnac = polymerizer.new_residue("g").unwrap();
        let invalid_bond =
            polymerizer.bond_groups("Peptide", &mut alanine, c_terminal, &mut glcnac, n_terminal);
        assert_miette_snapshot!(invalid_bond);
        // When bonding fails due to the acceptor, make sure that the donor remains untouched
        assert!(matches!(
            alanine.group_state(&c_terminal).unwrap(),
            GroupState::Free
        ));

        let nonexistent_bond =
            polymerizer.bond_groups("Super", &mut murnac, lactyl, &mut alanine, n_terminal);
        assert_miette_snapshot!(nonexistent_bond);
    }

    // FIXME: Needs finishing!
    #[test]
    fn modify() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut murnac = polymerizer.new_residue("m").unwrap();
        assert_eq!(murnac.monoisotopic_mass(), dec!(293.11106657336));

        polymerizer.modify("Anh", &mut murnac).unwrap();
        assert_eq!(murnac.monoisotopic_mass(), dec!(275.10050188933));
        let reducing_end = FunctionalGroup::new("Hydroxyl", "Reducing End");
        assert!(matches!(
            murnac.group_state(&reducing_end).unwrap(),
            GroupState::Modified(_)
        ));

        // let modify_non_free_group = polymerizer.modify("Anh", &mut murnac);
        // assert_miette_snapshot!(modify_non_free_group);

        // // FIXME: When no targets are found, simply report the functional groups and states on the target

        // // Start a new polymer by resetting the polymerizer
        // // FIXME: I can make this work, by enumerating the functional groups from the residue provided, looking them
        // // up in the index, and seeing if they contain the right ID
        // let mut polymerizer = polymerizer.reset();
        // let residue_from_wrong_polymer = polymerizer.modify("Anh", &mut murnac);
        // assert_miette_snapshot!(residue_from_wrong_polymer);

        // let nonexistent_modification = polymerizer.modify("Arg", &mut murnac);
        // assert_miette_snapshot!(nonexistent_modification);
    }

    // FIXME: Needs finishing!
    #[test]
    fn bond() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut murnac = polymerizer.new_residue("m").unwrap();
        let mut alanine = polymerizer.new_residue("A").unwrap();
        assert_eq!(
            murnac.monoisotopic_mass() + alanine.monoisotopic_mass(),
            dec!(382.15874504254)
        );

        polymerizer.bond("Stem", &mut murnac, &mut alanine).unwrap();
        assert_eq!(
            murnac.monoisotopic_mass() + alanine.monoisotopic_mass(),
            dec!(364.14818035851)
        );
        let lactyl = FunctionalGroup::new("Carboxyl", "Lactyl Ether");
        let n_terminal = FunctionalGroup::new("Amino", "N-Terminal");
        assert!(matches!(
            murnac.group_state(&lactyl).unwrap(),
            GroupState::Donor(_)
        ));
        assert!(matches!(
            alanine.group_state(&n_terminal).unwrap(),
            GroupState::Acceptor
        ));

        // let groups_not_free =
        //     polymerizer.bond_groups("Stem", &mut murnac, lactyl, &mut alanine, n_terminal);
        // assert_miette_snapshot!(groups_not_free);

        // let c_terminal = FunctionalGroup::new("Carboxyl", "C-Terminal");
        // let mut glcnac = polymerizer.new_residue("g").unwrap();
        // let invalid_bond =
        //     polymerizer.bond_groups("Peptide", &mut alanine, c_terminal, &mut glcnac, n_terminal);
        // assert_miette_snapshot!(invalid_bond);
        // // When bonding fails due to the acceptor, make sure that the donor remains untouched
        // assert!(matches!(
        //     alanine.group_state(&c_terminal).unwrap(),
        //     GroupState::Free
        // ));

        // let nonexistent_bond =
        //     polymerizer.bond_groups("Super", &mut murnac, lactyl, &mut alanine, n_terminal);
        // assert_miette_snapshot!(nonexistent_bond);
    }
}
