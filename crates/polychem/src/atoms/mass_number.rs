use std::num::NonZeroU32;

use crate::MassNumber;

impl MassNumber {
    pub(crate) fn new(n: u32) -> Option<Self> {
        NonZeroU32::new(n).map(Self)
    }
}
