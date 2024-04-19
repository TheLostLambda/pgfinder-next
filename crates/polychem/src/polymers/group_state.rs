use std::fmt::{self, Display, Formatter};

use crate::GroupState;

impl Display for GroupState {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                GroupState::Free => "free",
                GroupState::Modified(_) => "modified",
                GroupState::Donor(_) => "a donor",
                GroupState::Acceptor(_) => "an acceptor",
            }
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_and_is_variant_derives() {
        let free = GroupState::default();
        assert!(free.is_free());
        assert!(!free.is_modified());
        assert!(!free.is_donor());
        assert!(!free.is_acceptor());
    }
}
