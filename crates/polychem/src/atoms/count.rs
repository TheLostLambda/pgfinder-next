use std::{
    fmt::{self, Display, Formatter},
    num::NonZeroU32,
    ops::Mul,
};

use rust_decimal::Decimal;

use crate::{Charge, Count, Mass};

impl Count {
    pub(crate) fn new(n: u32) -> Option<Self> {
        NonZeroU32::new(n).map(Self)
    }
}

impl Mul<Mass> for Count {
    type Output = Mass;

    fn mul(self, rhs: Mass) -> Self::Output {
        Mass(Decimal::from(self.0.get()) * rhs.0)
    }
}

impl Mul<Charge> for Count {
    type Output = Charge;

    fn mul(self, rhs: Charge) -> Self::Output {
        Charge(i64::from(self.0.get()) * rhs.0)
    }
}

impl Display for Count {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let count = self.0.get();
        if count > 1 {
            write!(f, "{}", count)?;
        }
        Ok(())
    }
}

impl Default for Count {
    fn default() -> Self {
        Self(NonZeroU32::new(1).unwrap())
    }
}
