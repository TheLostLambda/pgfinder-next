use std::fmt::{self, Display, Formatter};

use crate::MassNumber;

impl Display for MassNumber {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl From<u32> for MassNumber {
    fn from(value: u32) -> Self {
        Self(value)
    }
}
