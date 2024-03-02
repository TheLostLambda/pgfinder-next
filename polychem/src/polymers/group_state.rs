use std::fmt::{Display, Formatter};

use crate::GroupState;

impl GroupState<'_, '_> {
    pub fn is_free(&self) -> bool {
        matches!(self, Self::Free)
    }
}

impl Display for GroupState<'_, '_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                GroupState::Free => "Free",
                GroupState::Modified(_) => "Modified",
                GroupState::Donor(_) => "Donor",
                GroupState::Acceptor => "Acceptor",
            }
        )
    }
}
