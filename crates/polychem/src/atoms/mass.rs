use rust_decimal::Decimal;

use crate::{Charge, Mass, Mz};

impl Mass {
    pub(crate) fn with_charge(self, charge: Charge) -> Option<Mz> {
        (!charge.is_zero()).then(|| Mz(self.0 / Decimal::from(charge.0)))
    }
}

impl From<Mass> for Decimal {
    fn from(value: Mass) -> Self {
        value.0
    }
}
