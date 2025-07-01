use std::{
    fmt::{self, Display, Formatter},
    num::NonZero,
    ops::Mul,
};

use rust_decimal::Decimal;
use static_assertions::const_assert;
use thiserror::Error;

use crate::{AverageMass, Charge, Count, Mass, MonoisotopicMass};

impl Count {
    const BITS: u32 = NonZero::<u32>::BITS;

    pub fn new(n: u32) -> Option<Self> {
        NonZero::new(n).map(Self)
    }
}

impl From<Count> for u32 {
    fn from(value: Count) -> Self {
        value.0.get()
    }
}

#[derive(Clone, Debug, Error)]
#[error("`Count` values must be non-zero")]
pub struct Zero;

impl TryFrom<u32> for Count {
    type Error = Zero;

    fn try_from(value: u32) -> Result<Self, Self::Error> {
        Self::new(value).ok_or(Zero)
    }
}

// NOTE: Since this will only panic on 16-bit platforms, I couldn't test a `TryFrom` impl if I wanted to — I don't want
// to add an entire error-handling code path that could never possibly be run.
#[expect(clippy::fallible_impl_from)]
impl From<Count> for usize {
    fn from(value: Count) -> Self {
        const_assert!(usize::BITS >= Count::BITS);
        // SAFETY: The above assertion prevents compilation on platforms with fewer bits than the type used to
        // represent `Count`. If this code compiles, then this `.unwrap()` will never panic.
        Self::try_from(value.0.get()).unwrap()
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
        Self::new(1).unwrap()
    }
}
