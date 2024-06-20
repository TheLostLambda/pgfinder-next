use std::num::NonZero;

use crate::MassNumber;

impl MassNumber {
    pub(crate) fn new(n: u32) -> Option<Self> {
        NonZero::new(n).map(Self)
    }
}
