use std::fmt::{self, Display, Formatter};

use crate::GroupState;

impl GroupState<'_, '_> {
    #[must_use]
    pub const fn is_free(&self) -> bool {
        matches!(self, Self::Free)
    }
}

impl Display for GroupState<'_, '_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                GroupState::Free => "free",
                GroupState::Modified(_) => "modified",
                GroupState::Donor(_) => "a donor",
                GroupState::Acceptor => "an acceptor",
            }
        )
    }
}
