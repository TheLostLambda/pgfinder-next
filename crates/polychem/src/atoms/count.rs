use std::{
    fmt::{self, Display, Formatter},
    num::NonZeroU32,
    ops::Mul,
};

use rust_decimal::Decimal;

use crate::{AverageMass, Charge, Count, Mass, MonoisotopicMass};

impl Count {
    pub(crate) fn new(n: u32) -> Option<Self> {
        NonZeroU32::new(n).map(Self)
    }
}

macro_rules! mass_mul_impls {
    // NOTE: `$mass_type` is a `tt` since it actually has to play the role of both a type (`ty`) and expression (`expr`)
    // in this impl, and `tt` appears to be the only way to pull off that sort of "metavariable polymorphism"
    ($($mass_type:tt),+ $(,)?) => {
        $(
            impl Mul<$mass_type> for Count {
                type Output = $mass_type;

                fn mul(self, rhs: $mass_type) -> Self::Output {
                    $mass_type(Decimal::from(self.0.get()) * rhs.0)
                }
            }
        )+
    };
}

mass_mul_impls!(Mass, MonoisotopicMass, AverageMass);

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
            write!(f, "{count}")?;
        }
        Ok(())
    }
}

impl Default for Count {
    fn default() -> Self {
        Self(NonZeroU32::new(1).unwrap())
    }
}
