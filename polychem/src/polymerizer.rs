use std::{collections::HashMap, fmt::Display, slice};

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

    pub fn residue(&mut self, abbr: impl AsRef<str>) -> Result<Residue<'a, 'p>> {
        self.residue_counter += 1;
        let residue = Residue::new(self.polymer_db, abbr, self.residue_counter)?;

        // NOTE: This assumes that all functional groups returned by `Residue::new()` start free!
        for &group in residue.functional_groups.keys() {
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
}

// FIXME: Add header for private section!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
impl<'a, 'p> Polymerizer<'a, 'p> {
    fn find_free_group<T: Into<Target<&'p str>>>(
        &self,
        targets: &(impl IntoIterator<Item = T> + Copy),
        residue: &Residue<'a, 'p>,
        group: Option<FunctionalGroup<'p>>,
    ) -> Result<FunctionalGroup<'p>, PolymerizerError> {
        let free_groups: Vec<_> = self.free_residue_groups(targets, residue).collect();
        match (group, &free_groups[..]) {
            (None, []) => Err(self.diagnose_missing_target(targets, residue)),
            (None, &[found_group]) => Ok(found_group),
            (None, found_groups) => Err(PolymerizerError::ambiguous_group(residue, found_groups)),
            (Some(target_group), found_groups) if found_groups.contains(&target_group) => {
                Ok(target_group)
            }
            (Some(target_group), _) => {
                Err(self.diagnose_missing_target_group(targets, residue, target_group))
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

        // SAFETY: These `.unwrap()`s might panic if the target hasn't first been validated by `self.find_free_group()`!
        self.free_group_index
            .get_mut(current_target)
            .unwrap()
            .insert(target.id(), group_state.is_free());

        let target_state = target.group_state_mut(&target_group).unwrap();
        *target_state = group_state;
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
            .find_free_group(&targets, target, target_group)
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
            .find_free_group(&slice::from_ref(from), donor, donor_group)
            .map_err(|e| PolychemError::bond(kind, donor, acceptor, "donor", e))?;
        let acceptor_group = self
            .find_free_group(&slice::from_ref(to), acceptor, acceptor_group)
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

    // FIXME: Of these group-fetching methods, `*_groups()`, some should be made public!
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

    fn diagnose_missing_target<T: Into<Target<&'p str>>>(
        &self,
        targets: &(impl IntoIterator<Item = T> + Copy),
        residue: &Residue<'a, 'p>,
    ) -> PolymerizerError {
        let non_free_groups: Vec<_> = self.residue_groups(targets, residue).collect();
        let residue_has_targeted_group = targets.into_iter().any(|possible_target| {
            let possible_target = possible_target.into();
            residue
                .functional_groups
                .keys()
                .any(|&fg| Target::from_residue_and_group(residue, fg).matches(&possible_target))
        });

        if !non_free_groups.is_empty() {
            PolymerizerError::all_groups_occupied(residue, &non_free_groups)
        } else if residue_has_targeted_group {
            PolymerizerError::residue_not_in_polymer(residue)
        } else {
            PolymerizerError::no_matching_groups(residue, targets)
        }
    }

    fn diagnose_missing_target_group<T: Into<Target<&'p str>>>(
        &self,
        targets: &(impl IntoIterator<Item = T> + Copy),
        residue: &Residue<'a, 'p>,
        group: FunctionalGroup<'p>,
    ) -> PolymerizerError {
        let current_target = Target::from_residue_and_group(residue, group);

        let theoretically_possible = targets
            .into_iter()
            .any(|possible_target| current_target.matches(&possible_target.into()));
        let group_in_index = self
            .residue_groups(&[current_target], residue)
            .next()
            .is_some();
        let residue_has_group = residue.functional_groups.contains_key(&group);

        if !theoretically_possible {
            PolymerizerError::invalid_target(targets, &current_target)
        } else if group_in_index {
            PolymerizerError::group_occupied(group, residue)
        } else if residue_has_group {
            PolymerizerError::residue_not_in_polymer(residue)
        } else {
            PolymerizerError::nonexistent_group(group, residue)
        }
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

// FIXME: Should probably break this error handling out into a sub-module...
#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
pub(crate) enum PolymerizerError {
    #[error("no functional groups on residue {0} matching the target were free: {1}")]
    AllGroupsOccupied(Id, String),

    #[error("the functional group {0} of residue {1} was already {2}, but must be free")]
    GroupOccupied(String, Id, String),

    #[error("no functional groups on residue {0} matched the target {1}")]
    NoMatchingGroups(Id, String),

    #[error("the functional group {0} does not exist on residue {1}")]
    NonexistentGroup(String, Id),

    #[error("expected a target matching {0}, got {1}")]
    InvalidTarget(String, Target),

    #[error("residue {0} does not belong to the current polymer")]
    #[diagnostic(help(
        "the referenced residue was likely created by a different Polymerizer instance"
    ))]
    ResidueNotInPolymer(Id),

    #[error("residue {0} contains more than one free target group: {1}")]
    #[diagnostic(help("to resolve this ambiguity, specify an exact functional group to modify"))]
    AmbiguousGroup(Id, String),
}

impl PolymerizerError {
    fn all_groups_occupied(residue: &Residue, groups: &[(FunctionalGroup, bool)]) -> Self {
        let groups_with_states = Self::comma_list(
            groups.iter().map(|(fg, _)| {
                let fg_state = residue.group_state(fg).unwrap();
                format!("{fg} is {fg_state}")
            }),
            "and",
        );
        Self::AllGroupsOccupied(residue.id(), groups_with_states)
    }

    fn group_occupied(group: FunctionalGroup, residue: &Residue) -> Self {
        Self::GroupOccupied(
            group.to_string(),
            residue.id(),
            residue.group_state(&group).unwrap().to_string(),
        )
    }

    fn no_matching_groups<'a, T: Into<Target<&'a str>>>(
        residue: &Residue,
        valid_targets: &(impl IntoIterator<Item = T> + Copy),
    ) -> Self {
        let valid_targets = Self::comma_list(valid_targets.into_iter().map(Into::into), "or");
        Self::NoMatchingGroups(residue.id(), valid_targets)
    }

    fn nonexistent_group(group: FunctionalGroup, residue: &Residue) -> Self {
        Self::NonexistentGroup(group.to_string(), residue.id())
    }

    fn invalid_target<'a, T: Into<Target<&'a str>>>(
        valid_targets: &(impl IntoIterator<Item = T> + Copy),
        target: impl Into<Target>,
    ) -> Self {
        let valid_targets = Self::comma_list(valid_targets.into_iter().map(Into::into), "or");
        Self::InvalidTarget(valid_targets, target.into())
    }

    const fn residue_not_in_polymer(residue: &Residue) -> Self {
        Self::ResidueNotInPolymer(residue.id())
    }

    fn ambiguous_group(residue: &Residue, groups: &[FunctionalGroup]) -> Self {
        let groups = Self::comma_list(groups, "and");
        Self::AmbiguousGroup(residue.id(), groups)
    }

    // FIXME: No clue where this belongs...
    fn comma_list<I: Display>(items: impl IntoIterator<Item = I>, final_sep: &str) -> String {
        let mut items: Vec<_> = items.into_iter().map(|i| i.to_string()).collect();
        let len = items.len();
        if len > 1 {
            items.sort_unstable();
            let last = format!("{final_sep} {}", items.last().unwrap());
            *items.last_mut().unwrap() = last;
        }
        let sep = if len > 2 { ", " } else { " " };
        items.join(sep)
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_ron_snapshot;
    use itertools::Itertools;
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::{
        atoms::atomic_database::AtomicDatabase,
        polymers::{polymer_database::PolymerDatabase, target::Target},
        testing_tools::assert_miette_snapshot,
        FunctionalGroup, GroupState, Massive,
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
        let residues = STEM_RESIDUES.map(|abbr| polymerizer.residue(abbr).unwrap());
        assert_ron_snapshot!(residues, {
            ".**.isotopes, .**.functional_groups" => insta::sorted_redaction()
        });

        let residues = STEM_RESIDUES.map(|abbr| polymerizer.residue(abbr).unwrap().id());
        assert_eq!(residues, [5, 6, 7, 8]);

        let mut polymerizer = polymerizer.reset();
        let residues = STEM_RESIDUES.map(|abbr| polymerizer.residue(abbr).unwrap().id());
        assert_eq!(residues, [1, 2, 3, 4]);

        let nonexistent_residue = polymerizer.residue("?");
        assert_miette_snapshot!(nonexistent_residue);
    }

    #[test]
    fn find_free_group() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut murnac = polymerizer.residue("m").unwrap();

        let carboxyl = Target::new("Carboxyl", None, None);
        let carboxyl_group = polymerizer
            .find_free_group(&[carboxyl], &murnac, None)
            .unwrap();
        assert_eq!(
            carboxyl_group,
            FunctionalGroup::new("Carboxyl", "Lactyl Ether")
        );

        let hydroxyl = Target::new("Hydroxyl", None, None);
        let ambiguous_group = polymerizer.find_free_group(&[hydroxyl], &murnac, None);
        assert_miette_snapshot!(ambiguous_group);

        let murnac_groups = murnac.functional_groups.keys().copied().collect_vec();
        for group in murnac_groups {
            polymerizer.update_group(&mut murnac, group, GroupState::Acceptor);
        }
        let all_groups_occupied = polymerizer.find_free_group(&[hydroxyl], &murnac, None);
        assert_miette_snapshot!(all_groups_occupied);

        // Start a new polymer by resetting the polymerizer
        let mut polymerizer = polymerizer.reset();
        let residue_not_in_polymer = polymerizer.find_free_group(&[hydroxyl], &murnac, None);
        assert_miette_snapshot!(residue_not_in_polymer);

        let mut murnac = polymerizer.residue("m").unwrap();
        let amino = Target::new("Amino", None, None);
        let crazy = Target::new("Crazy", None, None);
        let no_matching_groups = polymerizer.find_free_group(&[amino, crazy], &murnac, None);
        assert_miette_snapshot!(no_matching_groups);

        let nonreducing_end = FunctionalGroup::new("Hydroxyl", "Nonreducing End");
        let invalid_target = polymerizer.find_free_group(&[crazy], &murnac, Some(nonreducing_end));
        assert_miette_snapshot!(invalid_target);

        polymerizer.update_group(&mut murnac, nonreducing_end, GroupState::Acceptor);
        let group_occupied =
            polymerizer.find_free_group(&[hydroxyl], &murnac, Some(nonreducing_end));
        assert_miette_snapshot!(group_occupied);

        polymerizer.update_group(&mut murnac, nonreducing_end, GroupState::Free);
        let hydroxyl_group = polymerizer
            .find_free_group(&[hydroxyl], &murnac, Some(nonreducing_end))
            .unwrap();
        assert_eq!(
            hydroxyl_group,
            FunctionalGroup::new("Hydroxyl", "Nonreducing End")
        );

        // Start a new polymer by resetting the polymerizer
        let mut polymerizer = polymerizer.reset();
        let residue_not_in_polymer =
            polymerizer.find_free_group(&[hydroxyl], &murnac, Some(nonreducing_end));
        assert_miette_snapshot!(residue_not_in_polymer);

        let alanine = polymerizer.residue("A").unwrap();
        let nonexistent_group =
            polymerizer.find_free_group(&[hydroxyl], &alanine, Some(nonreducing_end));
        assert_miette_snapshot!(nonexistent_group);
    }

    #[test]
    fn modify_group() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut murnac = polymerizer.residue("m").unwrap();
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

        let mut alanine = polymerizer.residue("A").unwrap();
        let nonexistent_group = polymerizer.modify_group("Red", &mut alanine, reducing_end);
        assert_miette_snapshot!(nonexistent_group);

        let nonexistent_modification = polymerizer.modify_group("Arg", &mut murnac, reducing_end);
        assert_miette_snapshot!(nonexistent_modification);
    }

    #[test]
    fn bond_groups() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut murnac = polymerizer.residue("m").unwrap();
        let mut alanine = polymerizer.residue("A").unwrap();
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
        let mut glcnac = polymerizer.residue("g").unwrap();
        let invalid_bond =
            polymerizer.bond_groups("Peptide", &mut alanine, c_terminal, &mut glcnac, n_terminal);
        assert_miette_snapshot!(invalid_bond);
        // When bonding fails due to the acceptor, make sure that the donor remains untouched
        assert!(alanine.group_state(&c_terminal).unwrap().is_free());

        let nonexistent_bond =
            polymerizer.bond_groups("Super", &mut murnac, lactyl, &mut alanine, n_terminal);
        assert_miette_snapshot!(nonexistent_bond);
    }

    #[test]
    fn modify() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut murnac = polymerizer.residue("m").unwrap();
        assert_eq!(murnac.monoisotopic_mass(), dec!(293.11106657336));

        polymerizer.modify("Anh", &mut murnac).unwrap();
        assert_eq!(murnac.monoisotopic_mass(), dec!(275.10050188933));
        let reducing_end = FunctionalGroup::new("Hydroxyl", "Reducing End");
        assert!(matches!(
            murnac.group_state(&reducing_end).unwrap(),
            GroupState::Modified(_)
        ));

        let all_groups_occupied = polymerizer.modify("Anh", &mut murnac);
        assert_miette_snapshot!(all_groups_occupied);

        // Start a new polymer by resetting the polymerizer
        let mut polymerizer = polymerizer.reset();
        let residue_not_in_polymer = polymerizer.modify("Anh", &mut murnac);
        assert_miette_snapshot!(residue_not_in_polymer);

        let mut murnac = polymerizer.residue("m").unwrap();
        let no_matching_groups = polymerizer.modify("Am", &mut murnac);
        assert_miette_snapshot!(no_matching_groups);

        let nonexistent_modification = polymerizer.modify("Arg", &mut murnac);
        assert_miette_snapshot!(nonexistent_modification);
    }

    // FIXME: Needs finishing!
    #[test]
    fn bond() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut murnac = polymerizer.residue("m").unwrap();
        let mut alanine = polymerizer.residue("A").unwrap();
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

        let all_groups_occupied = polymerizer.bond("Stem", &mut murnac, &mut alanine);
        assert_miette_snapshot!(all_groups_occupied);

        // Start a new polymer by resetting the polymerizer
        let mut polymerizer = polymerizer.reset();
        let residue_not_in_polymer = polymerizer.bond("Stem", &mut murnac, &mut alanine);
        assert_miette_snapshot!(residue_not_in_polymer);

        let mut alanine = polymerizer.residue("A").unwrap();
        let mut glcnac = polymerizer.residue("g").unwrap();
        let no_matching_groups = polymerizer.bond("Peptide", &mut alanine, &mut glcnac);
        assert_miette_snapshot!(no_matching_groups);
        // When bonding fails due to the acceptor, make sure that the donor remains untouched
        let c_terminal = FunctionalGroup::new("Carboxyl", "C-Terminal");
        assert!(alanine.group_state(&c_terminal).unwrap().is_free());

        let nonexistent_bond = polymerizer.bond("Super", &mut murnac, &mut alanine);
        assert_miette_snapshot!(nonexistent_bond);
    }
}
