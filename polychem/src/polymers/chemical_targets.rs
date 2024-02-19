// Standard Library Imports
use std::{
    collections::{HashMap, HashSet},
    fmt::Display,
    iter,
    ops::Deref,
};

// Public API ==========================================================================================================

#[derive(Copy, Clone, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
pub(super) struct Target<S: Deref<Target = str> = String> {
    group: S,
    location: Option<S>,
    residue: Option<S>,
}

impl<S: Deref<Target = str>> Target<S> {
    pub(super) const fn new(group: S, location: Option<S>, residue: Option<S>) -> Self {
        Self {
            group,
            location,
            residue,
        }
    }
}

#[derive(Clone, Eq, PartialEq, Debug)]
pub(super) struct TargetIndex<'a, V = ()> {
    index: GroupMap<'a, V>,
}

impl<'a, V> TargetIndex<'a, V> {
    pub(super) fn new() -> Self {
        Self {
            index: HashMap::new(),
        }
    }

    pub(super) fn insert(&mut self, target: impl Into<Target<&'a str>>, value: V) -> Option<V> {
        let Target {
            group,
            location,
            residue,
        } = target.into();
        self.index
            .entry(group)
            .or_default()
            .entry(location)
            .or_default()
            .insert(residue, value)
    }

    pub(super) fn get(&self, target: impl Into<Target<&'a str>>) -> Vec<&V> {
        self.get_entries(target)
            .into_iter()
            .map(|(_, v)| v)
            .collect()
    }

    pub(super) fn get_residues(&'a self, target: impl Into<Target<&'a str>>) -> HashSet<&'a str> {
        self.get_entries(target)
            .into_iter()
            .filter_map(|(&r, _)| r)
            .collect()
    }
}

impl<'a> TargetIndex<'a> {
    pub(super) fn insert_residue(&mut self, target: impl Into<Target<&'a str>>) -> bool {
        self.insert(target, ()).is_some()
    }
}

// Target Printing and Borrowing =======================================================================================

impl<S: Deref<Target = str>> Display for Target<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", &*self.group)?;
        if let Some(location) = self.location.as_deref() {
            write!(f, " at={location:?}")?;
        }
        if let Some(residue) = self.residue.as_deref() {
            write!(f, " of={residue:?}")?;
        }
        Ok(())
    }
}

impl<'a> From<&'a Target> for Target<&'a str> {
    fn from(value: &'a Target) -> Self {
        Self {
            group: &value.group,
            location: value.location.as_deref(),
            residue: value.residue.as_deref(),
        }
    }
}

// Collecting an Iterator of Targets Into a TargetIndex ================================================================

impl<'a, V, K: Into<Target<&'a str>>> FromIterator<(K, V)> for TargetIndex<'a, V> {
    fn from_iter<T: IntoIterator<Item = (K, V)>>(iter: T) -> Self {
        let mut index = Self::new();
        for (k, v) in iter {
            index.insert(k, v);
        }
        index
    }
}

impl<'a, K: Into<Target<&'a str>>> FromIterator<K> for TargetIndex<'a> {
    fn from_iter<T: IntoIterator<Item = K>>(iter: T) -> Self {
        iter.into_iter().zip(iter::repeat(())).collect()
    }
}

// Private Types and Methods ===========================================================================================

type GroupMap<'a, T> = HashMap<&'a str, LocationMap<'a, T>>;
type LocationMap<'a, T> = HashMap<Option<&'a str>, ResidueMap<'a, T>>;
type ResidueMap<'a, T> = HashMap<Option<&'a str>, T>;

impl<'a, V> TargetIndex<'a, V> {
    // NOTE: Is this is ever made public, I should actually apply this lint — ignoring it for now keeps the code for
    // this function simpler and just adds one `&` to `.get_residues()`
    #[allow(clippy::ref_option_ref)]
    fn get_entries(&self, target: impl Into<Target<&'a str>>) -> Vec<(&Option<&str>, &V)> {
        let Target {
            group,
            location,
            residue,
        } = target.into();

        // NOTE: Although this approach repeats Vec::new, it avoids the nesting of the equivalent if-let approach
        let Some(locations) = self.index.get(group) else {
            return Vec::new();
        };
        if location.is_none() {
            return locations.values().flatten().collect();
        }
        let Some(residues) = locations.get(&location) else {
            return Vec::new();
        };
        if residue.is_none() {
            return residues.iter().collect();
        }
        residues
            .get_key_value(&residue)
            .map_or_else(Vec::new, |entry| vec![entry])
    }
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;

    use super::{Target, TargetIndex};

    static TARGET_LIST: Lazy<[(Target<&str>, &str); 3]> = Lazy::new(|| {
        [
            (Target::new("Amino", None, None), "group"),
            (
                Target::new("Amino", Some("N-Terminal"), None),
                "group-location",
            ),
            (
                Target::new("Amino", Some("N-Terminal"), Some("Alanine")),
                "group-location-residue",
            ),
        ]
    });

    #[test]
    fn print_targets() {
        assert_eq!(TARGET_LIST[0].0.to_string(), r#""Amino""#);
        assert_eq!(TARGET_LIST[1].0.to_string(), r#""Amino" at="N-Terminal""#);
        assert_eq!(
            TARGET_LIST[2].0.to_string(),
            r#""Amino" at="N-Terminal" of="Alanine""#
        );
    }

    #[test]
    fn construct_target_index() {
        let mut for_index = TargetIndex::new();
        for &(target, value) in TARGET_LIST.iter() {
            for_index.insert(target, value);
        }
        let iter_index: TargetIndex<_> = TARGET_LIST.iter().copied().collect();
        assert_eq!(for_index, iter_index);
    }

    #[test]
    fn get_nonexistent_group() {
        let index: TargetIndex<_> = TARGET_LIST.iter().copied().collect();
        let amino = Target::new("Carboxyl", None, None);
        assert!(index.get(amino).is_empty());
    }

    #[test]
    fn get_group() {
        let index: TargetIndex<_> = TARGET_LIST.iter().copied().collect();
        let amino = Target::new("Amino", None, None);
        let mut values = index.get(amino);
        values.sort_unstable();
        assert_eq!(
            values,
            vec![&"group", &"group-location", &"group-location-residue"]
        );
    }

    #[test]
    fn get_group_location() {
        let index: TargetIndex<_> = TARGET_LIST.iter().copied().collect();
        let amino = Target::new("Amino", Some("N-Terminal"), None);
        let mut values = index.get(amino);
        values.sort_unstable();
        assert_eq!(values, vec![&"group-location", &"group-location-residue"]);
    }

    #[test]
    fn get_group_location_residue() {
        let index: TargetIndex<_> = TARGET_LIST.iter().copied().collect();
        let amino = Target::new("Amino", Some("N-Terminal"), Some("Alanine"));
        let mut values = index.get(amino);
        values.sort_unstable();
        assert_eq!(values, vec![&"group-location-residue"]);
    }

    static RESIDUE_LIST: Lazy<[Target<&str>; 6]> = Lazy::new(|| {
        [
            Target::new("Amino", Some("N-Terminal"), Some("Lysine")),
            Target::new("Amino", Some("Sidechain"), Some("Lysine")),
            Target::new("Carboxyl", Some("C-Terminal"), Some("Lysine")),
            Target::new("Amino", Some("N-Terminal"), Some("Aspartic Acid")),
            Target::new("Carboxyl", Some("Sidechain"), Some("Aspartic Acid")),
            Target::new("Carboxyl", Some("C-Terminal"), Some("Aspartic Acid")),
        ]
    });

    #[test]
    fn construct_residue_index() {
        let mut for_index = TargetIndex::new();
        for &residue in RESIDUE_LIST.iter() {
            for_index.insert_residue(residue);
        }
        let iter_index: TargetIndex = RESIDUE_LIST.iter().copied().collect();
        assert_eq!(for_index, iter_index);
    }

    #[test]
    fn get_nonexistant_residues() {
        let index: TargetIndex = RESIDUE_LIST.iter().copied().collect();
        let amino = Target::new("Hydroxyl", None, None);
        let values = index.get_residues(amino);
        assert!(values.is_empty());
    }

    #[test]
    fn get_amino_residues() {
        let index: TargetIndex = RESIDUE_LIST.iter().copied().collect();
        let amino = Target::new("Amino", None, None);
        let mut values = Vec::from_iter(index.get_residues(amino));
        values.sort_unstable();
        assert_eq!(values, vec!["Aspartic Acid", "Lysine"]);
    }

    #[test]
    fn get_carboxyl_residues() {
        let index: TargetIndex = RESIDUE_LIST.iter().copied().collect();
        let amino = Target::new("Carboxyl", None, None);
        let mut values = Vec::from_iter(index.get_residues(amino));
        values.sort_unstable();
        assert_eq!(values, vec!["Aspartic Acid", "Lysine"]);
    }

    #[test]
    fn get_amino_sidechain_residues() {
        let index: TargetIndex = RESIDUE_LIST.iter().copied().collect();
        let amino = Target::new("Amino", Some("Sidechain"), None);
        let mut values = Vec::from_iter(index.get_residues(amino));
        values.sort_unstable();
        assert_eq!(values, vec!["Lysine"]);
    }

    #[test]
    fn get_carboxyl_sidechain_residues() {
        let index: TargetIndex = RESIDUE_LIST.iter().copied().collect();
        let amino = Target::new("Carboxyl", Some("Sidechain"), None);
        let mut values = Vec::from_iter(index.get_residues(amino));
        values.sort_unstable();
        assert_eq!(values, vec!["Aspartic Acid"]);
    }

    #[test]
    fn get_aspartic_acid_n_terminal() {
        let index: TargetIndex = RESIDUE_LIST.iter().copied().collect();
        let amino = Target::new("Amino", Some("N-Terminal"), Some("Aspartic Acid"));
        let mut values = Vec::from_iter(index.get_residues(amino));
        values.sort_unstable();
        assert_eq!(values, vec!["Aspartic Acid"]);
    }

    #[test]
    fn get_aspartic_acid_nonexistant_terminal() {
        let index: TargetIndex = RESIDUE_LIST.iter().copied().collect();
        let amino = Target::new("Carboxyl", Some("N-Terminal"), Some("Aspartic Acid"));
        let values = index.get_residues(amino);
        assert!(values.is_empty());
    }

    #[test]
    fn get_nonexistant_amino_n_terminal() {
        let index: TargetIndex = RESIDUE_LIST.iter().copied().collect();
        let amino = Target::new("Amino", Some("N-Terminal"), Some("Glycine"));
        let values = index.get_residues(amino);
        assert!(values.is_empty());
    }
}
