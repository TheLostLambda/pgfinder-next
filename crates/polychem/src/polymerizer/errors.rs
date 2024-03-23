use std::fmt::Display;

use ahash::HashSet;
use miette::Diagnostic;
use thiserror::Error;

use crate::{polymers::target::Target, FunctionalGroup, Id, Residue};

#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
#[error(transparent)]
#[diagnostic(transparent)]
pub struct Error(#[from] PolymerizerError);

// FIXME: Should probably break this error handling out into a sub-module...
#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
pub enum PolymerizerError {
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
}

impl PolymerizerError {
    pub(super) fn all_groups_occupied(
        residue: &Residue,
        groups: &[(FunctionalGroup, bool)],
    ) -> Self {
        let groups_with_states = Self::comma_list(
            groups.iter().map(|(fg, _)| {
                let fg_state = residue.group_state(fg).unwrap();
                format!("{fg} is {fg_state}")
            }),
            "and",
        );
        Self::AllGroupsOccupied(residue.id(), groups_with_states)
    }

    pub(super) fn group_occupied(group: FunctionalGroup, residue: &Residue) -> Self {
        Self::GroupOccupied(
            group.to_string(),
            residue.id(),
            residue.group_state(&group).unwrap().to_string(),
        )
    }

    pub(super) fn no_matching_groups<'a, T: Into<Target<&'a str>>>(
        residue: &Residue,
        valid_targets: &(impl IntoIterator<Item = T> + Copy),
    ) -> Self {
        let valid_targets = Self::comma_list(valid_targets.into_iter().map(Into::into), "or");
        Self::NoMatchingGroups(residue.id(), valid_targets)
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
