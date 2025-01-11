use ahash::HashSet;
use miette::Diagnostic;
use std::fmt::Display;
use thiserror::Error;

use crate::{FunctionalGroup, GroupState, ResidueId, moieties::target::Target};

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Diagnostic, Error)]
pub enum FindFreeGroupsError {
    #[error("residue {residue_id} contains more than {needed_groups} free target group{}: {}",
        if .needed_groups == &1 { "" } else { "s" },
        comma_list(.free_group_names, "and")
    )]
    #[diagnostic(help(
        "to resolve this ambiguity, specify the exact functional group{} to modify",
        if .needed_groups == &1 { "" } else { "s" }
    ))]
    AmbiguousGroups {
        residue_id: ResidueId,
        needed_groups: usize,
        free_group_names: Vec<String>,
    },

    #[error("encountered a duplicate target group: {group_name}")]
    #[diagnostic(help("check for typos, or try removing the duplicate target group"))]
    DuplicateTargetGroup { group_name: String },

    #[error(
        "the functional group {group_name} of residue {residue_id} was already {group_state}, but must be free"
    )]
    GroupOccupied {
        group_name: String,
        residue_id: ResidueId,
        group_state: GroupState,
    },

    #[error(
        "{} functional groups on residue {residue_id} matching the target were free{}: {}",
        if .needed_groups > &1 { "too few" } else { "no" },
        if .needed_groups > &1 { format!(" ({}/{})", .free_groups, .needed_groups) } else { String::new() },
        comma_list(.group_names_and_states.iter().map(|(fg, gs)| format!("{fg} is {gs}")), "and")
    )]
    GroupsOccupied {
        residue_id: ResidueId,
        free_groups: usize,
        needed_groups: usize,
        group_names_and_states: Vec<(String, GroupState)>,
    },

    #[error(
        "expected a target group matching {}, but got {target_group}",
        comma_list(.valid_targets, "or")
    )]
    InvalidTargetGroup {
        valid_targets: Vec<Target>,
        target_group: Target,
    },

    #[error("failed to find multiple of the requested free groups")]
    MultipleMissingFreeGroups {
        #[related]
        errors: Vec<Self>,
    },

    #[error("no target groups were provided, but you must provide at least one group to look for")]
    NoTargetGroups,

    #[error("the functional group {group_name} does not exist on residue {residue_id}")]
    NonexistentGroup {
        group_name: String,
        residue_id: ResidueId,
    },

    #[error("{} functional groups on residue {residue_id} matched the target {}{}",
        if .needed_groups > &1 { "too few" } else { "no" },
        comma_list(.valid_targets, "or"),
        if .needed_groups > &1 {
            format!(" ({}/{}): {}", .free_groups, .needed_groups, comma_list(.free_group_names, "and"))
        } else {
            String::new()
        }
    )]
    TooFewMatchingGroups {
        residue_id: ResidueId,
        valid_targets: Vec<Target>,
        free_groups: usize,
        needed_groups: usize,
        free_group_names: Vec<String>,
    },

    #[error("attempted to find zero free groups, but you must look for at least one")]
    ZeroGroupNumber,
}

impl FindFreeGroupsError {
    pub fn ambiguous_groups(
        residue_id: ResidueId,
        needed_groups: usize,
        free_groups: &HashSet<FunctionalGroup>,
    ) -> Self {
        let mut free_group_names: Vec<_> = free_groups.iter().map(ToString::to_string).collect();
        free_group_names.sort_unstable();

        Self::AmbiguousGroups {
            residue_id,
            needed_groups,
            free_group_names,
        }
    }

    pub fn duplicate_target_group(group: &FunctionalGroup) -> Self {
        let group_name = group.to_string();

        Self::DuplicateTargetGroup { group_name }
    }

    pub fn group_occupied(
        group: &FunctionalGroup,
        residue_id: ResidueId,
        group_state: GroupState,
    ) -> Self {
        let group_name = group.to_string();

        Self::GroupOccupied {
            group_name,
            residue_id,
            group_state,
        }
    }

    pub fn groups_occupied(
        residue_id: ResidueId,
        needed_groups: usize,
        groups_and_states: Vec<(FunctionalGroup, GroupState)>,
    ) -> Self {
        let mut group_names_and_states: Vec<_> = groups_and_states
            .into_iter()
            .map(|(fg, gs)| (fg.to_string(), gs))
            .collect();
        group_names_and_states.sort_unstable();

        let free_groups = group_names_and_states
            .iter()
            .filter(|(_, gs)| gs.is_free())
            .count();

        Self::GroupsOccupied {
            residue_id,
            free_groups,
            needed_groups,
            group_names_and_states,
        }
    }

    pub fn invalid_target_group<'t, 'p, T: 'p>(
        valid_targets: &'t [T],
        target_group: Target<&'p str>,
    ) -> Self
    where
        &'t T: Into<Target<&'p str>>,
    {
        let mut valid_targets: Vec<_> = valid_targets
            .iter()
            .map(|t| Target::from(t.into()))
            .collect();
        valid_targets.sort_unstable();

        let target_group = Target::from(target_group);

        Self::InvalidTargetGroup {
            valid_targets,
            target_group,
        }
    }

    pub fn multiple_missing_free_groups(mut errors: Vec<Self>) -> Self {
        errors.sort_unstable();

        Self::MultipleMissingFreeGroups { errors }
    }

    pub fn nonexistent_group(group: &FunctionalGroup, residue_id: ResidueId) -> Self {
        let group_name = group.to_string();

        Self::NonexistentGroup {
            group_name,
            residue_id,
        }
    }

    pub fn too_few_matching_groups<'t, 'p, T: 'p>(
        residue_id: ResidueId,
        valid_targets: &'t [T],
        needed_groups: usize,
        free_groups: &HashSet<FunctionalGroup>,
    ) -> Self
    where
        &'t T: Into<Target<&'p str>>,
    {
        let mut valid_targets: Vec<_> = valid_targets
            .iter()
            .map(|t| Target::from(t.into()))
            .collect();
        valid_targets.sort_unstable();

        let mut free_group_names: Vec<_> = free_groups.iter().map(ToString::to_string).collect();
        free_group_names.sort_unstable();

        let free_groups = free_group_names.len();

        Self::TooFewMatchingGroups {
            residue_id,
            valid_targets,
            free_groups,
            needed_groups,
            free_group_names,
        }
    }
}

// FIXME: I'm not certain where this code belongs...
pub fn comma_list(items: impl IntoIterator<Item: Display>, final_sep: &str) -> String {
    let mut items: Vec<_> = items.into_iter().map(|i| i.to_string()).collect();
    let len = items.len();
    if len > 1 {
        let last_item = items.last_mut().unwrap();
        let last_item_with_sep = format!("{final_sep} {last_item}");
        *last_item = last_item_with_sep;
    }
    let sep = if len > 2 { ", " } else { " " };
    items.join(sep)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_comma_list() {
        let empty: [&str; 0] = [];
        assert_eq!(comma_list(empty, "whatever"), "");
        assert_eq!(comma_list(["A"], "whatever"), "A");
        assert_eq!(comma_list(["A", "B"], "or"), "A or B");
        assert_eq!(comma_list(["B", "C", "A"], "and"), "B, C, and A");
    }
}
