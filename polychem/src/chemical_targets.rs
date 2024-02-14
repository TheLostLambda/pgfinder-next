use std::collections::HashMap;

#[derive(Clone, PartialEq, Eq, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
pub(super) struct Target {
    group: String,
    location: Option<String>,
    residue: Option<String>,
}

impl Target {
    pub fn new<S: Into<String>>(group: S, location: Option<S>, residue: Option<S>) -> Self {
        Self {
            group: group.into(),
            location: location.map(|s| s.into()),
            residue: residue.map(|s| s.into()),
        }
    }
}

type GroupMap<T> = HashMap<String, LocationMap<T>>;
type LocationMap<T> = HashMap<Option<String>, ResidueMap<T>>;
type ResidueMap<T> = HashMap<Option<String>, T>;

#[derive(Clone, PartialEq, Eq, Debug)]
pub(super) struct TargetIndex<T = ()> {
    // FIXME: From a pure time-complexity standpoint, this should be ideal, but I honestly don't expect more than a
    // couple dozen entries to ever end up here — it'll likely be even fewer, so maybe just searching a (sorted?)
    // Vec could be better?
    index: GroupMap<T>,
}

// NOTE: This is manually implemented because #[derive(Default)] is wrongly convinced that T needs to implement Default
impl<T> Default for TargetIndex<T> {
    fn default() -> Self {
        Self {
            index: HashMap::new(),
        }
    }
}

impl<K: Into<Target>, V> FromIterator<(K, V)> for TargetIndex<V> {
    fn from_iter<T: IntoIterator<Item = (K, V)>>(iter: T) -> Self {
        let mut index = Self::new();
        for (k, v) in iter {
            index.insert_value(k, v);
        }
        index
    }
}

impl<T> TargetIndex<T> {
    fn new() -> Self {
        Self::default()
    }

    fn insert_value(&mut self, target: impl Into<Target>, value: T) -> Option<T> {
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

    fn get_entries(&self, target: impl Into<Target>) -> Vec<(&Option<String>, &T)> {
        let Target {
            group,
            location,
            residue,
        } = target.into();
        let Some(index) = self.index.get(&group) else {
            return Vec::new();
        };
        let Some(index) = location.and_then(|l| index.get(&Some(l))) else {
            return index.values().flatten().collect();
        };
        if let Some(entry) = residue.and_then(|r| index.get_key_value(&Some(r))) {
            vec![entry]
        } else {
            index.iter().collect()
        }
    }

    fn get_residues(&self, target: impl Into<Target>) -> Vec<&String> {
        self.get_entries(target)
            .into_iter()
            .filter_map(|(r, _)| r.as_ref())
            .collect()
    }

    fn get(&self, target: impl Into<Target>) -> Vec<&T> {
        self.get_entries(target)
            .into_iter()
            .map(|(_, v)| v)
            .collect()
    }
    // FIXME: for get_residues, just drop None from the Vec
}

impl TargetIndex {
    fn insert(&mut self, target: impl Into<Target>) -> bool {
        self.insert_value(target, ()).is_some()
    }
}

#[cfg(test)]
mod tests {
    use once_cell::sync::Lazy;

    use super::{Target, TargetIndex};

    static TARGET_LIST: Lazy<[(Target, &str); 3]> = Lazy::new(|| {
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
    fn construct_target_index() {
        let mut for_index = TargetIndex::new();
        for (target, value) in TARGET_LIST.clone() {
            for_index.insert_value(target.clone(), value);
        }
        let iter_index: TargetIndex<_> = TARGET_LIST.iter().cloned().collect();
        assert_eq!(for_index, iter_index);
    }

    #[test]
    fn get_nonexistent_group() {
        let index: TargetIndex<_> = TARGET_LIST.iter().cloned().collect();
        let amino = Target::new("Carboxyl", None, None);
        assert!(index.get(amino).is_empty());
    }

    #[test]
    fn get_group() {
        let index: TargetIndex<_> = TARGET_LIST.iter().cloned().collect();
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
        let index: TargetIndex<_> = TARGET_LIST.iter().cloned().collect();
        let amino = Target::new("Amino", Some("N-Terminal"), None);
        let mut values = index.get(amino);
        values.sort_unstable();
        assert_eq!(values, vec![&"group-location", &"group-location-residue"]);
    }

    #[test]
    fn get_group_location_residue() {
        let index: TargetIndex<_> = TARGET_LIST.iter().cloned().collect();
        let amino = Target::new("Amino", Some("N-Terminal"), Some("Alanine"));
        let mut values = index.get(amino);
        values.sort_unstable();
        assert_eq!(values, vec![&"group-location-residue"]);
    }
}
