use crate::{Abundance, Charge, Mass, Mz};
use std::ops::Mul;

impl Mass {
    pub(crate) fn checked_div(self, charge: Charge) -> Option<Mz> {
        (charge.0 != 0).then(|| self / charge)
    }
}

impl Mul<Abundance> for Mass {
    type Output = Mass;

    fn mul(self, rhs: Abundance) -> Self::Output {
        Mass(self.0 * rhs.0)
    }
}
