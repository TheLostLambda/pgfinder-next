use std::ops::Div;

use rust_decimal::Decimal;

use crate::{Charge, Mass, Mz};

impl Div<Charge> for Mass {
    type Output = Mz;

    fn div(self, rhs: Charge) -> Self::Output {
        Mz(self.0 / Decimal::from(rhs.0))
    }
}
