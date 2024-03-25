use std::fmt::Display;

use ahash::HashSet;
use miette::Diagnostic;
use thiserror::Error;

use crate::{polymers::target::Target, FunctionalGroup, Id, Residue};

#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
#[error(transparent)]
#[diagnostic(transparent)]
pub struct Error(#[from] PolymerizerError);

// FIXME: All of these `derive`s need reordering...
#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Ord, PartialOrd, Error)]
pub enum PolymerizerError {
    #[error(
        "{} functional groups on residue {0} matching the target were free{}: {3}",
        if .2 > &1 { "too few" } else { "no" },
        if .2 > &1 { format!(" ({}/{})", .1, .2) } else { String::new() }
    )]
    GroupsOccupied(Id, usize, usize, String),

    #[error("the functional group {0} of residue {1} was already {2}, but must be free")]
    GroupOccupied(String, Id, String),

    #[error("{} functional groups on residue {0} matched the target {1}{}",
        if .3 > &1 { "too few" } else { "no" },
        if .3 > &1 { format!(" ({}/{}): {}", .2, .3, .4) } else { String::new() }
    )]
    TooFewMatchingGroups(Id, String, usize, usize, String),

    #[error("the functional group {0} does not exist on residue {1}")]
    NonexistentGroup(String, Id),

    #[error("expected a target matching {0}, got {1}")]
    InvalidTarget(String, Target),

    #[error("residue {0} does not belong to the current polymer")]
    #[diagnostic(help(
        "the referenced residue was likely created by a different Polymerizer instance"
    ))]
    ResidueNotInPolymer(Id),

    #[error("residue {0} contains more than {1} free target group{}: {2}", if .1 == &1 { "" } else { "s" })]
    #[diagnostic(help(
        "to resolve this ambiguity, specify the exact functional group{} to modify",
        if .1 == &1 { "" } else { "s" }
    ))]
    AmbiguousGroups(Id, usize, String),

    #[error("attemped to find zero free groups, but you must look for at least one")]
    ZeroGroupNumber,

    #[error("received an empty group set, but you must provide at least one group to look for")]
    EmptyGroupSet,

    #[error("failed to find multiple of the requested free groups")]
    MultipleMissingFreeGroups(#[related] Vec<PolymerizerError>),
}

// FIXME: If this were made an impl on `Error`, then I could make the `PolymerizerError` private and these methods
// would be the only way to construct these errors (you couldn't build them manually accidentially)
impl PolymerizerError {
    pub(super) fn groups_occupied(
        residue: &Residue,
        groups: &[(FunctionalGroup, bool)],
        number: usize,
    ) -> Self {
        let free_groups = groups.iter().filter(|(_, free)| *free).count();
        let groups_with_states = Self::comma_list(
            groups.iter().map(|(fg, _)| {
                let fg_state = residue.group_state(fg).unwrap();
                format!("{fg} is {fg_state}")
            }),
            "and",
        );
        Self::GroupsOccupied(residue.id(), free_groups, number, groups_with_states)
    }

    pub(super) fn group_occupied(group: FunctionalGroup, residue: &Residue) -> Self {
        Self::GroupOccupied(
            group.to_string(),
            residue.id(),
            residue.group_state(&group).unwrap().to_string(),
        )
    }

    pub(super) fn too_few_matching_groups<'a, T: Into<Target<&'a str>>>(
        residue: &Residue,
        valid_targets: &(impl IntoIterator<Item = T> + Copy),
        found_groups: &[FunctionalGroup],
        number: usize,
    ) -> Self {
        let valid_targets = Self::comma_list(valid_targets.into_iter().map(Into::into), "or");
        let number_found = found_groups.len();
        let found_groups = Self::comma_list(found_groups, "and");
        Self::TooFewMatchingGroups(
            residue.id(),
            valid_targets,
            number_found,
            number,
            found_groups,
        )
    }

    pub(super) fn nonexistent_group(group: FunctionalGroup, residue: &Residue) -> Self {
        Self::NonexistentGroup(group.to_string(), residue.id())
    }

    pub(super) fn invalid_target<'a, T: Into<Target<&'a str>>>(
        valid_targets: &(impl IntoIterator<Item = T> + Copy),
        target: impl Into<Target>,
    ) -> Self {
        let valid_targets = Self::comma_list(valid_targets.into_iter().map(Into::into), "or");
        Self::InvalidTarget(valid_targets, target.into())
    }

    pub(super) const fn residue_not_in_polymer(residue: &Residue) -> Self {
        Self::ResidueNotInPolymer(residue.id())
    }

    pub(super) fn ambiguous_groups(
        residue: &Residue,
        number: usize,
        groups: HashSet<FunctionalGroup>,
    ) -> Self {
        let groups = Self::comma_list(groups, "and");
        Self::AmbiguousGroups(residue.id(), number, groups)
    }

    pub(super) const fn zero_group_number() -> Self {
        Self::ZeroGroupNumber
    }

    pub(super) const fn empty_group_set() -> Self {
        Self::EmptyGroupSet
    }

    pub(super) fn multiple_missing_free_groups(mut errors: Vec<Self>) -> Self {
        errors.sort_unstable();
        Self::MultipleMissingFreeGroups(errors)
    }

    // FIXME: No clue where this belongs...
    pub(super) fn comma_list<I: Display>(
        items: impl IntoIterator<Item = I>,
        final_sep: &str,
    ) -> String {
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
