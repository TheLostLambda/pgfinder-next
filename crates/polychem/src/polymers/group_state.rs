use std::fmt::{self, Display, Formatter};

use crate::GroupState;

impl GroupState<'_, '_> {
    #[must_use]
    pub const fn is_free(&self) -> bool {
        matches!(self, Self::Free)
    }

    #[must_use]
    pub const fn is_modified(&self) -> bool {
        matches!(self, Self::Modified(_))
    }

    #[must_use]
    pub const fn is_donor(&self) -> bool {
        matches!(self, Self::Donor(_))
    }

    #[must_use]
    pub const fn is_acceptor(&self) -> bool {
        matches!(self, Self::Acceptor)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_free() {
        let free = GroupState::Free;
        assert!(free.is_free());
        assert!(!free.is_modified());
        assert!(!free.is_donor());
        assert!(!free.is_acceptor());
    }
}
