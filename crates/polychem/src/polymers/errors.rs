use ahash::HashSet;
use miette::Diagnostic;
use std::fmt::Display;
use thiserror::Error;

use crate::{moieties::target::Target, FunctionalGroup, GroupState, ResidueId};

#[derive(Debug, Clone, Diagnostic, Error)]
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

    #[error("attemped to find zero free groups, but you must look for at least one")]
    ZeroGroupNumber,
}

impl FindFreeGroupsError {
    pub(crate) fn ambiguous_groups(
        residue_id: ResidueId,
        needed_groups: usize,
        free_groups: &HashSet<FunctionalGroup>,
    ) -> Self {
        let free_group_names = free_groups.iter().map(ToString::to_string).collect();
        Self::AmbiguousGroups {
            residue_id,
            needed_groups,
            free_group_names,
        }
    }

    pub(crate) fn groups_occupied(
        residue_id: ResidueId,
        needed_groups: usize,
        groups_and_states: Vec<(FunctionalGroup, GroupState)>,
    ) -> Self {
        let group_names_and_states: Vec<_> = groups_and_states
            .into_iter()
            .map(|(fg, gs)| (fg.to_string(), gs))
            .collect();
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

    pub(crate) fn too_few_matching_groups<'p, T: Into<Target<&'p str>>>(
        residue_id: ResidueId,
        valid_targets: impl IntoIterator<Item = T>,
        needed_groups: usize,
        free_groups: &HashSet<FunctionalGroup>,
    ) -> Self {
        let valid_targets = valid_targets
            .into_iter()
            .map(|t| Target::from(t.into()))
            .collect();
        let free_group_names: Vec<_> = free_groups.iter().map(ToString::to_string).collect();
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
pub fn comma_list<I: Display>(items: impl IntoIterator<Item = I>, final_sep: &str) -> String {
    let mut items: Vec<_> = items.into_iter().map(|i| i.to_string()).collect();
    let len = items.len();
    if len > 1 {
        items.sort_unstable();
        let last_item = items.last_mut().unwrap();
        let last_item_with_sep = format!("{final_sep} {last_item}");
        *last_item = last_item_with_sep;
    }
    let sep = if len > 2 { ", " } else { " " };
    items.join(sep)
}
