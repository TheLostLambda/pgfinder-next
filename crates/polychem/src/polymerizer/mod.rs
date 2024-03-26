// FIXME: Good god, this file needs a lot of work...
mod errors;
pub(crate) use errors::{Error, PolymerizerError};

use std::{cmp::Ordering, slice};

use ahash::{HashMap, HashSet};

use crate::{
    atoms::atomic_database::AtomicDatabase,
    polymers::{
        polymer_database::{BondDescription, ModificationDescription, PolymerDatabase},
        target::{Index, Target},
    },
    AnyMod, AnyModification, Bond, BondTarget, Count, FunctionalGroup, GroupState, Id, NamedMod,
    PolychemError, Residue, Result,
};

#[derive(Clone)]
pub struct Polymerizer<'a, 'p> {
    atomic_db: &'a AtomicDatabase,
    polymer_db: &'p PolymerDatabase<'a>,
    residue_counter: Id,
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

    #[must_use]
    pub const fn atomic_db(&self) -> &'a AtomicDatabase {
        self.atomic_db
    }

    #[must_use]
    pub const fn polymer_db(&self) -> &'p PolymerDatabase<'a> {
        self.polymer_db
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

    pub fn chain(
        &mut self,
        abbrs: &[impl AsRef<str>],
        bond_kind: impl AsRef<str>,
    ) -> Result<Vec<Residue<'a, 'p>>> {
        let mut residues: Vec<_> = abbrs
            .iter()
            .map(|abbr| self.residue(abbr))
            .collect::<Result<_, _>>()?;
        self.bond_chain(&mut residues, bond_kind)?;
        Ok(residues)
    }

    pub fn modify(
        &mut self,
        modification: impl Into<AnyModification<'a, 'p>>,
        target: &mut Residue<'a, 'p>,
    ) -> Result<()> {
        let modification = modification.into();
        match modification.kind {
            AnyMod::Named(m) => {
                self.modify_with_optional_groups(m.abbr(), target, modification.multiplier)
            }
            AnyMod::Offset(_) => todo!(),
        }
    }

    pub fn modify_only_group(
        &mut self,
        abbr: impl AsRef<str>,
        target: &mut Residue<'a, 'p>,
    ) -> Result<()> {
        self.modify_with_optional_groups(abbr, target, 1)
    }

    pub fn modify_only_groups(
        &mut self,
        abbr: impl AsRef<str>,
        target: &mut Residue<'a, 'p>,
        number: Count,
    ) -> Result<()> {
        self.modify_with_optional_groups(abbr, target, number)
    }

    // PERF: Could create an `_unchecked` version for when you've already called `self.free_*_groups()` — skip straight
    // to `self.update_group()`!
    pub fn modify_group(
        &mut self,
        abbr: impl AsRef<str>,
        target: &mut Residue<'a, 'p>,
        target_group: FunctionalGroup<'p>,
    ) -> Result<()> {
        self.modify_with_optional_groups(abbr, target, [target_group])
    }

    // TODO: Add modify_groups() which takes a slice of target groups

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

    pub fn bond_chain(
        &mut self,
        residues: &mut [Residue<'a, 'p>],
        bond_kind: impl AsRef<str>,
    ) -> Result<()> {
        let bond_kind = bond_kind.as_ref();

        // NOTE: Doing this properly requires a `windows_mut()` method, which is blocked on lending iterators, but GATs
        // have now been stabilized, so the way is clear for those. Keep an eye out for standard library updates! For
        // now, this manual indexing and pattern-matching is a work-around!
        for i in 0..residues.len() - 1 {
            // SAFETY: The `unreachable!()` is safe, since `residues[i..=i + 1]` will always have two items in it
            let [donor, acceptor] = &mut residues[i..=i + 1] else {
                unreachable!()
            };

            self.bond(bond_kind, donor, acceptor)?;
        }

        Ok(())
    }
}

// FIXME: Add header for private section!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// FIXME: The groups stuff could be it's own submodule!
enum Groups<'p> {
    Any(usize),
    These(HashSet<FunctionalGroup<'p>>),
}

impl Groups<'_> {
    fn validate(self) -> Result<Self, PolymerizerError> {
        match self {
            Groups::Any(number) if number < 1 => Err(PolymerizerError::zero_group_number()),
            Groups::These(groups) if groups.is_empty() => Err(PolymerizerError::empty_group_set()),
            valid => Ok(valid),
        }
    }
}

// NOTE: Since this will only panic on 16-bit platforms, I couldn't test a `TryFrom` impl if I wanted to — I don't want
// to add an entire error-handling code path that could never possibly be run. I'm intentionally using `u32` instead of
// the more semantically correct `Count` here so that I get a compilation error if the type of `Count` is ever changed,
// so I can re-evaluate if this panic is still impossible on the hardware this program can actually be run on.
#[allow(clippy::fallible_impl_from)]
impl From<u32> for Groups<'_> {
    fn from(value: u32) -> Self {
        // SAFETY: This could only possibly panic on platforms with pointers smaller than 32 bits. If someone gets this
        // program running on a 16-bit platform, they 1) deserve a medal and 2) will have much bigger problems than
        // this panic...
        Self::Any(usize::try_from(value).unwrap())
    }
}

// FIXME: Update comment!
// NOTE: I can't unify these two to take anything that impls IntoIterator, because Count (a u32) might someday implement
// IntoIterator, and then that definition would overlap with the one above — this is, yet again, an ugly bit of code
// that's waiting on the stabilization of specialization...
// FIXME: This will need to be more flexible about the type it takes... I'll want to take Vec, slices, and HashSets?
impl<'p, const N: usize> From<[FunctionalGroup<'p>; N]> for Groups<'p> {
    fn from(value: [FunctionalGroup<'p>; N]) -> Self {
        Self::These(HashSet::from_iter(value))
    }
}

impl<'p> From<Option<FunctionalGroup<'p>>> for Groups<'p> {
    fn from(value: Option<FunctionalGroup<'p>>) -> Self {
        // NOTE: Annoyingly, `HashSet::from` only works for the std `RandomState`? If I wanted to change `from_iter`
        // into just `from` here, I would need to use `AHashSet` instead...
        value.map_or_else(|| Self::Any(1), |_| Self::These(HashSet::from_iter(value)))
    }
}

impl<'a, 'p> Polymerizer<'a, 'p> {
    // FIXME: Needs a good look-over...
    fn find_free_groups<T: Into<Target<&'p str>>>(
        &self,
        targets: &(impl IntoIterator<Item = T> + Copy),
        residue: &Residue<'a, 'p>,
        groups: impl Into<Groups<'p>>,
    ) -> Result<HashSet<FunctionalGroup<'p>>, PolymerizerError> {
        // PERF: There is likely some small performance gains to be made by moving to `Vec`s instead of `HashSet`s — as
        // long all returned results come from `free_groups` which already returns unique results
        let free_groups: HashSet<_> = self.free_residue_groups(targets, residue).collect();
        let groups = groups.into().validate()?;
        match groups {
            Groups::Any(number) => match free_groups.len().cmp(&number) {
                Ordering::Less => Err(self.diagnose_missing_groups_any(targets, residue, number)),
                Ordering::Equal => Ok(free_groups),
                Ordering::Greater => Err(PolymerizerError::ambiguous_groups(
                    residue,
                    number,
                    free_groups,
                )),
            },
            Groups::These(groups) => {
                if groups.is_subset(&free_groups) {
                    Ok(groups)
                } else {
                    let not_free: Vec<_> = groups.difference(&free_groups).copied().collect();
                    Err(self.diagnose_missing_groups_these(targets, residue, &not_free))
                }
            }
        }
    }

    fn update_groups(
        &mut self,
        target: &mut Residue<'a, 'p>,
        target_groups: &HashSet<FunctionalGroup<'p>>,
        group_state: GroupState<'a, 'p>,
    ) {
        for target_group in target_groups {
            let current_target = Target::from_residue_and_group(target, *target_group);

            // SAFETY: These `.unwrap()`s might panic if the target hasn't first been validated by `self.find_free_group()`!
            self.free_group_index
                .get_mut(current_target)
                .unwrap()
                .insert(target.id(), group_state.is_free());

            let target_state = target.group_state_mut(target_group).unwrap();
            *target_state = group_state;
        }
    }

    fn modify_with_optional_groups(
        &mut self,
        abbr: impl AsRef<str>,
        target: &mut Residue<'a, 'p>,
        groups: impl Into<Groups<'p>>,
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

        let target_groups = self
            .find_free_groups(&targets, target, groups)
            .map_err(|e| PolychemError::modification(name, abbr, target, e))?;

        let modified_state = GroupState::Modified(NamedMod {
            abbr,
            name,
            lost,
            gained,
        });
        self.update_groups(target, &target_groups, modified_state);

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
        let donor_groups = self
            .find_free_groups(&slice::from_ref(from), donor, donor_group)
            .map_err(|e| PolychemError::bond(kind, donor, acceptor, "donor", e))?;
        let acceptor_groups = self
            .find_free_groups(&slice::from_ref(to), acceptor, acceptor_group)
            .map_err(|e| PolychemError::bond(kind, donor, acceptor, "acceptor", e))?;

        let donor_state = GroupState::Donor(Bond {
            kind,
            lost,
            acceptor: BondTarget {
                residue: acceptor.id(),
                // SAFETY: The `count` of 1 provided to `find_free_groups` ensures that `acceptor_groups` will contain
                // at least one element
                // FIXME: Awful! Maybe a `find_free_group` that sets `find_free_groups`s `count` to 1?
                group: *acceptor_groups.iter().next().unwrap(),
            },
        });
        self.update_groups(donor, &donor_groups, donor_state);
        self.update_groups(acceptor, &acceptor_groups, GroupState::Acceptor);

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

    fn diagnose_missing_groups_any<T: Into<Target<&'p str>>>(
        &self,
        targets: &(impl IntoIterator<Item = T> + Copy),
        residue: &Residue<'a, 'p>,
        number: usize,
    ) -> PolymerizerError {
        let non_free_groups: Vec<_> = self.residue_groups(targets, residue).collect();
        let targeted_residue_groups: Vec<_> = targets
            .into_iter()
            .flat_map(|possible_target| {
                let possible_target = possible_target.into();
                residue
                    .functional_groups
                    .keys()
                    .copied()
                    .filter(move |&fg| {
                        Target::from_residue_and_group(residue, fg).matches(&possible_target)
                    })
            })
            .collect();

        if non_free_groups.len() >= number {
            PolymerizerError::groups_occupied(residue, &non_free_groups, number)
        } else if targeted_residue_groups.len() != non_free_groups.len() {
            PolymerizerError::residue_not_in_polymer(residue)
        } else {
            PolymerizerError::too_few_matching_groups(
                residue,
                targets,
                &targeted_residue_groups,
                number,
            )
        }
    }

    fn diagnose_missing_groups_these<T: Into<Target<&'p str>>>(
        &self,
        targets: &(impl IntoIterator<Item = T> + Copy),
        residue: &Residue<'a, 'p>,
        groups: &[FunctionalGroup<'p>],
    ) -> PolymerizerError {
        let mut errors: Vec<_> = groups
            .iter()
            .map(|&group| {
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
            })
            .collect();
        if errors.len() == 1 {
            // SAFETY: We've just checked that there is one item is present in the Vec, so this will never panic
            errors.remove(0)
        } else {
            PolymerizerError::multiple_missing_free_groups(errors)
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

#[cfg(test)]
mod tests {
    use insta::assert_ron_snapshot;
    use once_cell::sync::Lazy;
    use rust_decimal::Decimal;
    use rust_decimal_macros::dec;

    use crate::{
        testing_tools::assert_miette_snapshot, Charged, Massive, Modification, OffsetKind,
        OffsetMod,
    };

    use super::*;

    const STEM_RESIDUES: [&str; 4] = ["A", "E", "J", "A"];

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);

    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "polymer_database.kdl",
            include_str!("../../tests/data/polymer_database.kdl"),
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
    fn chain() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut residues = STEM_RESIDUES.map(|abbr| polymerizer.residue(abbr).unwrap());
        polymerizer.bond_chain(&mut residues, "Peptide").unwrap();
        assert_ron_snapshot!(residues, {
            ".**.composition, .**.lost" => "<FORMULA>",
            ".**.functional_groups" => insta::sorted_redaction()
        });

        let all_in_one_residues = polymerizer.chain(&STEM_RESIDUES, "Peptide").unwrap();
        assert_eq!(
            residues
                .iter()
                .map(Massive::monoisotopic_mass)
                .sum::<Decimal>(),
            all_in_one_residues
                .iter()
                .map(Massive::monoisotopic_mass)
                .sum::<Decimal>(),
        );

        assert_eq!(
            residues
                .iter()
                .map(Massive::monoisotopic_mass)
                .sum::<Decimal>(),
            STEM_RESIDUES
                .iter()
                .map(|abbr| polymerizer.residue(abbr).unwrap().monoisotopic_mass())
                .sum::<Decimal>()
                + Modification::new(
                    u32::try_from(residues.len() - 1).unwrap(),
                    OffsetMod::new(&ATOMIC_DB, OffsetKind::Remove, "H2O").unwrap()
                )
                .monoisotopic_mass()
        );

        let nonexistent_bond = polymerizer.bond_chain(&mut residues, "?");
        assert_eq!(
            nonexistent_bond.clone().unwrap_err(),
            polymerizer.chain(&STEM_RESIDUES, "?").unwrap_err()
        );
        assert_miette_snapshot!(nonexistent_bond);

        let nonexistent_residue = polymerizer.chain(&["?"], "Peptide");
        assert_miette_snapshot!(nonexistent_residue);
    }

    #[test]
    fn find_single_unambiguous_free_groups() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut murnac = polymerizer.residue("m").unwrap();

        let carboxyl = Target::new("Carboxyl", None, None);
        let carboxyl_group = polymerizer
            .find_free_groups(&[carboxyl], &murnac, 1)
            .unwrap();
        assert_eq!(carboxyl_group.len(), 1);
        assert!(carboxyl_group.contains(&FunctionalGroup::new("Carboxyl", "Lactyl Ether")));

        let zero_groups = polymerizer.find_free_groups(&[carboxyl], &murnac, 0);
        assert_miette_snapshot!(zero_groups);

        let hydroxyl = Target::new("Hydroxyl", None, None);
        let ambiguous_group = polymerizer.find_free_groups(&[hydroxyl], &murnac, 1);
        assert_miette_snapshot!(ambiguous_group);

        let murnac_groups = murnac.functional_groups.keys().copied().collect();
        polymerizer.update_groups(&mut murnac, &murnac_groups, GroupState::Acceptor);
        let all_groups_occupied = polymerizer.find_free_groups(&[hydroxyl], &murnac, 1);
        assert_miette_snapshot!(all_groups_occupied);

        // Start a new polymer by resetting the polymerizer
        let mut polymerizer = polymerizer.reset();
        let residue_not_in_polymer = polymerizer.find_free_groups(&[hydroxyl], &murnac, 1);
        assert_miette_snapshot!(residue_not_in_polymer);

        let murnac = polymerizer.residue("m").unwrap();
        let amino = Target::new("Amino", None, None);
        let crazy = Target::new("Crazy", None, None);
        let no_matching_groups = polymerizer.find_free_groups(&[amino, crazy], &murnac, 1);
        assert_miette_snapshot!(no_matching_groups);
    }

    #[test]
    fn find_single_specific_free_groups() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut murnac = polymerizer.residue("m").unwrap();

        let nonreducing_end = FunctionalGroup::new("Hydroxyl", "Nonreducing End");
        let hydroxyl = Target::new("Hydroxyl", None, None);
        polymerizer.update_groups(
            &mut murnac,
            &HashSet::from_iter([nonreducing_end]),
            GroupState::Acceptor,
        );
        let group_occupied = polymerizer.find_free_groups(&[hydroxyl], &murnac, [nonreducing_end]);
        assert_miette_snapshot!(group_occupied);

        let empty_groups = polymerizer.find_free_groups(&[hydroxyl], &murnac, []);
        assert_miette_snapshot!(empty_groups);

        polymerizer.update_groups(
            &mut murnac,
            &HashSet::from_iter([nonreducing_end]),
            GroupState::Free,
        );
        let hydroxyl_group = polymerizer
            .find_free_groups(&[hydroxyl], &murnac, [nonreducing_end])
            .unwrap();
        assert_eq!(hydroxyl_group.len(), 1);
        assert!(hydroxyl_group.contains(&FunctionalGroup::new("Hydroxyl", "Nonreducing End")));

        // Start a new polymer by resetting the polymerizer
        let mut polymerizer = polymerizer.reset();
        let residue_not_in_polymer =
            polymerizer.find_free_groups(&[hydroxyl], &murnac, [nonreducing_end]);
        assert_miette_snapshot!(residue_not_in_polymer);

        let alanine = polymerizer.residue("A").unwrap();
        let nonexistent_group =
            polymerizer.find_free_groups(&[hydroxyl], &alanine, [nonreducing_end]);
        assert_miette_snapshot!(nonexistent_group);

        let crazy = Target::new("Crazy", None, None);
        let invalid_target = polymerizer.find_free_groups(&[crazy], &murnac, [nonreducing_end]);
        assert_miette_snapshot!(invalid_target);
    }

    #[test]
    fn find_multiple_unambiguous_free_groups() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut murnac = polymerizer.residue("m").unwrap();

        let hydroxyl = Target::new("Hydroxyl", None, None);
        let hydroxyl_groups = polymerizer
            .find_free_groups(&[hydroxyl], &murnac, 3)
            .unwrap();
        assert_eq!(hydroxyl_groups.len(), 3);
        for location in ["Reducing End", "Nonreducing End", "6-Position"] {
            assert!(hydroxyl_groups.contains(&FunctionalGroup::new("Hydroxyl", location)));
        }

        let too_few_hydroxyl_groups = polymerizer.find_free_groups(&[hydroxyl], &murnac, 4);
        assert_miette_snapshot!(too_few_hydroxyl_groups);

        let ambiguous_group = polymerizer.find_free_groups(&[hydroxyl], &murnac, 2);
        assert_miette_snapshot!(ambiguous_group);

        let nonreducing_end = FunctionalGroup::new("Hydroxyl", "Nonreducing End");
        polymerizer.update_groups(
            &mut murnac,
            &HashSet::from_iter([nonreducing_end]),
            GroupState::Acceptor,
        );
        let groups_occupied = polymerizer.find_free_groups(&[hydroxyl], &murnac, 3);
        assert_miette_snapshot!(groups_occupied);

        let remaining_hydroxyl_groups = polymerizer
            .find_free_groups(&[hydroxyl], &murnac, 2)
            .unwrap();
        assert_eq!(remaining_hydroxyl_groups.len(), 2);
        for location in ["Reducing End", "6-Position"] {
            assert!(remaining_hydroxyl_groups.contains(&FunctionalGroup::new("Hydroxyl", location)));
        }

        let lysine = polymerizer.residue("K").unwrap();
        let n_terminal = Target::new("Amino", Some("N-Terminal"), None);
        let carboxyl = Target::new("Carboxyl", None, None);
        let terminals = polymerizer
            .find_free_groups(&[n_terminal, carboxyl], &lysine, 2)
            .unwrap();
        assert_eq!(terminals.len(), 2);
        assert!(terminals.contains(&FunctionalGroup::new("Amino", "N-Terminal")));
        assert!(terminals.contains(&FunctionalGroup::new("Carboxyl", "C-Terminal")));
    }

    #[test]
    fn find_multiple_specific_free_groups() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut murnac = polymerizer.residue("m").unwrap();

        let hydroxyl = Target::new("Hydroxyl", None, None);
        let reducing_end = FunctionalGroup::new("Hydroxyl", "Reducing End");
        let nonreducing_end = FunctionalGroup::new("Hydroxyl", "Nonreducing End");
        let six_position = FunctionalGroup::new("Hydroxyl", "6-Position");

        let one_hydroxyl_group = polymerizer
            .find_free_groups(&[hydroxyl], &murnac, [reducing_end])
            .unwrap();
        assert_eq!(one_hydroxyl_group.len(), 1);
        assert!(one_hydroxyl_group.contains(&reducing_end));

        let two_hydroxyl_groups = polymerizer
            .find_free_groups(&[hydroxyl], &murnac, [reducing_end, nonreducing_end])
            .unwrap();
        assert_eq!(two_hydroxyl_groups.len(), 2);
        assert!(two_hydroxyl_groups.contains(&reducing_end));
        assert!(two_hydroxyl_groups.contains(&nonreducing_end));

        let all_hydroxyl_groups = polymerizer
            .find_free_groups(
                &[hydroxyl],
                &murnac,
                [reducing_end, nonreducing_end, six_position],
            )
            .unwrap();
        assert_eq!(all_hydroxyl_groups.len(), 3);
        assert!(all_hydroxyl_groups.contains(&reducing_end));
        assert!(all_hydroxyl_groups.contains(&nonreducing_end));
        assert!(all_hydroxyl_groups.contains(&six_position));

        let n_terminal = FunctionalGroup::new("Amino", "N-Terminal");
        let invalid_target = polymerizer.find_free_groups(
            &[hydroxyl],
            &murnac,
            [reducing_end, nonreducing_end, six_position, n_terminal],
        );
        assert_miette_snapshot!(invalid_target);

        let crazy = Target::new("Crazy", None, None);
        let still_invalid_target = polymerizer.find_free_groups(
            &[hydroxyl, crazy],
            &murnac,
            [reducing_end, nonreducing_end, six_position, n_terminal],
        );
        assert_miette_snapshot!(still_invalid_target);

        let amino = Target::new("Amino", None, None);
        let nonexistent_group = polymerizer.find_free_groups(
            &[hydroxyl, crazy, amino],
            &murnac,
            [reducing_end, nonreducing_end, six_position, n_terminal],
        );
        assert_miette_snapshot!(nonexistent_group);

        polymerizer.update_groups(
            &mut murnac,
            &HashSet::from_iter([nonreducing_end]),
            GroupState::Acceptor,
        );
        let multiple_errors = polymerizer.find_free_groups(
            &[hydroxyl, crazy, amino],
            &murnac,
            [reducing_end, nonreducing_end, six_position, n_terminal],
        );
        assert_miette_snapshot!(multiple_errors);
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
    fn modify_only_group() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut murnac = polymerizer.residue("m").unwrap();
        assert_eq!(murnac.monoisotopic_mass(), dec!(293.11106657336));

        polymerizer.modify_only_group("Anh", &mut murnac).unwrap();
        assert_eq!(murnac.monoisotopic_mass(), dec!(275.10050188933));
        let reducing_end = FunctionalGroup::new("Hydroxyl", "Reducing End");
        assert!(matches!(
            murnac.group_state(&reducing_end).unwrap(),
            GroupState::Modified(_)
        ));

        let all_groups_occupied = polymerizer.modify_only_group("Anh", &mut murnac);
        assert_miette_snapshot!(all_groups_occupied);

        // Start a new polymer by resetting the polymerizer
        let mut polymerizer = polymerizer.reset();
        let residue_not_in_polymer = polymerizer.modify_only_group("Anh", &mut murnac);
        assert_miette_snapshot!(residue_not_in_polymer);

        let mut murnac = polymerizer.residue("m").unwrap();
        let no_matching_groups = polymerizer.modify_only_group("Am", &mut murnac);
        assert_miette_snapshot!(no_matching_groups);

        let nonexistent_modification = polymerizer.modify_only_group("Arg", &mut murnac);
        assert_miette_snapshot!(nonexistent_modification);
    }

    #[test]
    fn modify_only_groups() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut murnac = polymerizer.residue("m").unwrap();
        assert_eq!(murnac.monoisotopic_mass(), dec!(293.11106657336));

        polymerizer
            .modify_only_groups("Met", &mut murnac, 3)
            .unwrap();
        assert_eq!(murnac.monoisotopic_mass(), dec!(335.15801676674));
        let reducing_end = FunctionalGroup::new("Hydroxyl", "Reducing End");
        let nonreducing_end = FunctionalGroup::new("Hydroxyl", "Nonreducing End");
        let six_position = FunctionalGroup::new("Hydroxyl", "6-Position");
        for hydroxyl_group in [reducing_end, nonreducing_end, six_position] {
            assert!(matches!(
                murnac.group_state(&hydroxyl_group).unwrap(),
                GroupState::Modified(_)
            ));
        }

        let mut murnac = polymerizer.residue("m").unwrap();
        assert_eq!(murnac.monoisotopic_mass(), dec!(293.11106657336));
        polymerizer
            .modify_only_groups("Ca", &mut murnac, 4)
            .unwrap();
        assert_eq!(murnac.charge(), 4);
        assert_eq!(murnac.monoisotopic_mass(), dec!(448.927935519603480));

        let mut alanine = polymerizer.residue("A").unwrap();
        assert_eq!(alanine.monoisotopic_mass(), dec!(89.04767846918));
        polymerizer
            .modify_only_groups("Ca", &mut alanine, 2)
            .unwrap();
        assert_eq!(alanine.charge(), 2);
        assert_eq!(alanine.monoisotopic_mass(), dec!(166.956112942301740));

        let mut murnac = polymerizer.residue("m").unwrap();
        polymerizer
            .modify_only_groups("Met", &mut murnac, 3)
            .unwrap();
        let all_groups_occupied = polymerizer.modify_only_groups("Ca", &mut murnac, 4);
        assert_miette_snapshot!(all_groups_occupied);
        let still_all_groups_occupied = polymerizer.modify_only_groups("Ca", &mut murnac, 2);
        assert_miette_snapshot!(still_all_groups_occupied);

        assert_eq!(murnac.monoisotopic_mass(), dec!(335.15801676674));
        polymerizer
            .modify_only_groups("Ca", &mut murnac, 1)
            .unwrap();
        assert_eq!(murnac.monoisotopic_mass(), dec!(374.112234003300870));

        // Start a new polymer by resetting the polymerizer
        let mut polymerizer = polymerizer.reset();
        let residue_not_in_polymer = polymerizer.modify_only_groups("Ca", &mut alanine, 2);
        assert_miette_snapshot!(residue_not_in_polymer);

        let mut murnac = polymerizer.residue("m").unwrap();
        let no_matching_groups = polymerizer.modify_only_groups("Anh", &mut murnac, 2);
        assert_miette_snapshot!(no_matching_groups);
    }

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
