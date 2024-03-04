// Standard Library Imports
use std::{
    collections::{hash_map::Entry, HashMap, HashSet},
    fmt::{Display, Formatter},
    iter,
    ops::Deref,
};

// Public API ==========================================================================================================

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
pub(crate) struct Target<S: Deref<Target = str> = String> {
    pub(crate) group: S,
    pub(crate) location: Option<S>,
    pub(crate) residue: Option<S>,
}

impl<S: Deref<Target = str>> Target<S> {
    pub(crate) const fn new(group: S, location: Option<S>, residue: Option<S>) -> Self {
        Self {
            group,
            location,
            residue,
        }
    }
}

impl<S: Deref<Target = str> + Eq> Target<S> {
    pub(crate) fn matches(&self, other: &Self) -> bool {
        self.group == other.group
            && (other.location.is_none() || self.location == other.location)
            && (other.residue.is_none() || self.residue == other.residue)
    }
}

#[derive(Clone, Eq, PartialEq, Debug)]
pub(crate) struct TargetIndex<'a, V = ()> {
    index: GroupMap<'a, V>,
}

impl<'a, V> TargetIndex<'a, V> {
    pub(crate) fn new() -> Self {
        Self {
            index: HashMap::new(),
        }
    }

    pub(crate) fn insert(&mut self, target: impl Into<Target<&'a str>>, value: V) -> Option<V> {
        match self.entry(target) {
            Entry::Occupied(mut e) => Some(e.insert(value)),
            Entry::Vacant(e) => {
                e.insert(value);
                None
            }
        }
    }

    pub(crate) fn entry(
        &mut self,
        target: impl Into<Target<&'a str>>,
    ) -> Entry<Option<&'a str>, V> {
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
            .entry(residue)
    }

    pub(crate) fn get_mut(&mut self, target: impl Into<Target<&'a str>>) -> Option<&mut V> {
        let Target {
            group,
            location,
            residue,
        } = target.into();
        self.index
            .get_mut(group)
            .and_then(|ls| ls.get_mut(&location))
            .and_then(|rs| rs.get_mut(&residue))
    }

    pub(crate) fn contains_target(&self, target: impl Into<Target<&'a str>>) -> bool {
        // PERF: Could be further optimized, see commit 9044078
        !self.matches_with_targets(target).is_empty()
    }

    pub(crate) fn matches(&self, target: impl Into<Target<&'a str>>) -> Vec<&V> {
        // PERF: Could also be faster if I never constructed the Targets I'm throwing away here
        self.matches_with_targets(target)
            .into_iter()
            .map(|(_, v)| v)
            .collect()
    }

    pub(crate) fn matches_with_targets(
        &self,
        target: impl Into<Target<&'a str>>,
    ) -> Vec<(Target<&'a str>, &V)> {
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
            return locations
                .iter()
                .flat_map(|(&l, rs)| rs.iter().map(move |(&r, v)| (Target::new(group, l, r), v)))
                .collect();
        }
        let Some(residues) = locations.get(&location) else {
            return Vec::new();
        };

        if residue.is_none() {
            return residues
                .iter()
                .map(|(&r, v)| (Target::new(group, location, r), v))
                .collect();
        }
        residues.get(&residue).map_or_else(Vec::new, |entry| {
            vec![(Target::new(group, location, residue), entry)]
        })
    }
}

// Target Printing and Borrowing =======================================================================================

impl<S: Deref<Target = str>> Display for Target<S> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
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

// Private Types =======================================================================================================

type GroupMap<'a, T> = HashMap<&'a str, LocationMap<'a, T>>;
type LocationMap<'a, T> = HashMap<Option<&'a str>, ResidueMap<'a, T>>;
type ResidueMap<'a, T> = HashMap<Option<&'a str>, T>;

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

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
    fn match_targets() {
        // More general targets match more specific ones
        assert!(TARGET_LIST[2].0.matches(&TARGET_LIST[1].0));
        assert!(TARGET_LIST[1].0.matches(&TARGET_LIST[0].0));

        // But not the other way around
        assert!(!TARGET_LIST[0].0.matches(&TARGET_LIST[1].0));
        assert!(!TARGET_LIST[1].0.matches(&TARGET_LIST[2].0));

        // Mismatched groups never match
        let matchless = Target::new("Hydroxyl", None, None);
        assert!(!TARGET_LIST[0].0.matches(&matchless));
        assert!(!TARGET_LIST[1].0.matches(&matchless));
        assert!(!TARGET_LIST[2].0.matches(&matchless));

        assert!(!matchless.matches(&TARGET_LIST[0].0));
        assert!(!matchless.matches(&TARGET_LIST[1].0));
        assert!(!matchless.matches(&TARGET_LIST[2].0));
    }

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
    fn insert_duplicate_target() {
        let mut iter_index: TargetIndex<_> = TARGET_LIST.iter().take(2).copied().collect();
        assert_eq!(
            iter_index.insert(TARGET_LIST[0].0, "new group"),
            Some("group")
        );
        assert_eq!(
            iter_index.insert(TARGET_LIST[1].0, "new group-location"),
            Some("group-location")
        );
        assert_eq!(
            iter_index.insert(TARGET_LIST[2].0, "first group-location-residue"),
            None
        );
        assert_eq!(
            iter_index.insert(TARGET_LIST[0].0, "final group"),
            Some("new group")
        );
        assert_eq!(
            iter_index.insert(TARGET_LIST[1].0, "final group-location"),
            Some("new group-location")
        );
        assert_eq!(
            iter_index.insert(TARGET_LIST[2].0, "new group-location-residue"),
            Some("first group-location-residue")
        );
    }

    #[test]
    fn get_nonexistent_group() {
        let index: TargetIndex<_> = TARGET_LIST.iter().copied().collect();
        let amino = Target::new("Carboxyl", None, None);
        assert!(index.matches(amino).is_empty());
    }

    #[test]
    fn get_group() {
        let index: TargetIndex<_> = TARGET_LIST.iter().copied().collect();
        let amino = Target::new("Amino", None, None);
        let mut values = index.matches(amino);
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
        let mut values = index.matches(amino);
        values.sort_unstable();
        assert_eq!(values, vec![&"group-location", &"group-location-residue"]);
    }

    #[test]
    fn get_group_location_residue() {
        let index: TargetIndex<_> = TARGET_LIST.iter().copied().collect();
        let amino = Target::new("Amino", Some("N-Terminal"), Some("Alanine"));
        let mut values = index.matches(amino);
        values.sort_unstable();
        assert_eq!(values, vec![&"group-location-residue"]);
    }

    #[test]
    fn update_values() {
        let mut index: TargetIndex<_> = TARGET_LIST.iter().copied().collect();
        let amino = Target::new("Amino", Some("N-Terminal"), Some("Alanine"));
        let mut values = index.matches(amino);
        values.sort_unstable();
        assert_eq!(values, vec![&"group-location-residue"]);

        *index.get_mut(amino).unwrap() = "residue-group-location?";
        let mut values = index.matches(amino);
        values.sort_unstable();
        assert_eq!(values, vec![&"residue-group-location?"]);
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

    impl<'a, V> TargetIndex<'a, V> {
        // NOTE: Probably not super useful public API, but it does make the logic of `TargetIndex::get_key_value` way
        // easier to test — it it does end up being useful outside of these tests, it can be made public again
        fn get_residues(&'a self, target: impl Into<Target<&'a str>>) -> HashSet<&'a str> {
            self.matches_with_targets(target)
                .into_iter()
                .filter_map(|(t, _)| t.residue)
                .collect()
        }
    }

    #[test]
    fn construct_residue_index() {
        let mut for_index = TargetIndex::new();
        for &residue in RESIDUE_LIST.iter() {
            for_index.insert(residue, ());
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
