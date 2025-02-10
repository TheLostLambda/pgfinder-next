use std::fmt::{self, Display, Formatter};

use crate::GroupState;

impl Display for GroupState {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Free => "free",
                Self::Modified(..) => "modified",
                Self::Donor(..) => "a donor",
                Self::Acceptor(..) => "an acceptor",
            }
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::{BondId, ModificationId};

    use super::*;

    #[test]
    fn display_group_state() {
        let states = [
            GroupState::Free,
            GroupState::Modified(ModificationId(0)),
            GroupState::Donor(BondId(0)),
            GroupState::Acceptor(BondId(0)),
        ];
        assert_eq!(
            states.map(|gs| gs.to_string()),
            [
                "free".to_owned(),
                "modified".to_owned(),
                "a donor".to_owned(),
                "an acceptor".to_owned()
            ]
        );
    }

    #[test]
    fn default_and_is_variant_derives() {
        let free = GroupState::default();
        assert!(free.is_free());
        assert!(!free.is_modified());
        assert!(!free.is_donor());
        assert!(!free.is_acceptor());
    }
}
