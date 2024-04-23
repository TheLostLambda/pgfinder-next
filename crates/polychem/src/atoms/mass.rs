use crate::{Abundance, AverageMass, Mass, MonoisotopicMass};
use std::ops::Mul;

macro_rules! mass_conversion_impls {
    ($($mass_type:ty),+ $(,)?) => {
        $(
            impl From<$mass_type> for Mass {
                fn from(value: $mass_type) -> Self {
                    Self(value.0)
                }
            }
            impl From<Mass> for $mass_type {
                fn from(value: Mass) -> Self {
                    Self(value.0)
                }
            }
        )+
    };
}

mass_conversion_impls!(MonoisotopicMass, AverageMass);

impl Mul<Abundance> for Mass {
    type Output = Mass;

    fn mul(self, rhs: Abundance) -> Self::Output {
        Mass(self.0 * rhs.0)
    }
}
