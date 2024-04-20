use std::ops::Div;

use rust_decimal::Decimal;

use crate::{AverageMass, AverageMz, Charge, MonoisotopicMass, MonoisotopicMz};

macro_rules! mz_div_impls {
    // NOTE: `$mz_type` is a `tt` since it actually has to play the role of both a type (`ty`) and expression (`expr`)
    // in this impl, and `tt` appears to be the only way to pull off that sort of "metavariable polymorphism"
    ($($mass_type:ty => $mz_type:tt),+ $(,)?) => {
        $(
            impl Div<Charge> for $mass_type {
                type Output = $mz_type;

                fn div(self, rhs: Charge) -> Self::Output {
                    $mz_type(self.0 / Decimal::from(rhs.0))
                }
            }
        )+
    };
}

mz_div_impls!(MonoisotopicMass => MonoisotopicMz, AverageMass => AverageMz);
