// Standard Library Imports
use std::{cmp::Ordering, ops::RangeInclusive};

// External Crate Imports
use derive_more::{From, Into};
use mzdata::mzpeaks::Tolerance;

// Public API ==========================================================================================================

// PERF: Try out `f32` and see if that saves enough space to speed things up!
#[derive(Copy, Clone, Debug, From, Into)]
pub struct TotalFloat(f64);
pub type Mz = TotalFloat;
pub type Minutes = TotalFloat;

impl Mz {
    pub fn ppm_window(mz: f64, ppm: f64) -> RangeInclusive<Self> {
        let (min_mz, max_mz) = Tolerance::PPM(ppm).bounds(mz);
        Self(min_mz)..=Self(max_mz)
    }
}

// Implementing Ord for TotalFloat  ====================================================================================

impl Ord for TotalFloat {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

impl PartialOrd for TotalFloat {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for TotalFloat {
    fn eq(&self, other: &Self) -> bool {
        // NOTE: This is *not* equivalent to `self.0 == other.0`! The `.total_cmp()` used in the `Ord` implementation
        // differs from the default `PartialEq` implementation for floats, and it's vital that the `Ord` and
        // `PartialEq` implementations agree on which values are equal.
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for TotalFloat {}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use std::ops::{Bound, RangeBounds};

    use assert_float_eq::assert_float_absolute_eq;

    use super::*;

    #[test]
    fn total_float_traits() {
        // Test `From` and `Into` impls
        let a: TotalFloat = 941.407_703.into();
        let b = TotalFloat::from(471.711_128);
        assert_float_absolute_eq!(a.into(), 941.407_703);
        assert_float_absolute_eq!(b.into(), 471.711_128);

        // Test `Copy` and `Ord` impl
        assert_eq!(a.cmp(&b), Ordering::Greater);
        assert_eq!(b.cmp(&a), Ordering::Less);
        assert_eq!(a.cmp(&a), Ordering::Equal);
        assert_eq!(b.cmp(&b), Ordering::Equal);
    }

    #[test]
    fn mz_ppm_window() {
        let window = Mz::ppm_window(941.407_703, 7.5);
        assert!(window.contains(&Mz::from(941.407_703)));

        let Bound::Included(&TotalFloat(start)) = window.start_bound() else {
            panic!("expected inclusive start bound");
        };
        let Bound::Included(&TotalFloat(end)) = window.end_bound() else {
            panic!("expected inclusive end bound")
        };
        assert_float_absolute_eq!(start, 941.400_642);
        assert_float_absolute_eq!(end, 941.414_764);
    }
}
