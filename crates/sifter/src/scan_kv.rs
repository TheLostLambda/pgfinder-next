// Standard Library Imports
use std::ops::RangeInclusive;

// External Crate Imports
use mzdata::mzpeaks::Tolerance;

// Local Crate Imports
use crate::{
    peaks::Peaks,
    total_float::{Minutes, Mz},
};

// Public API ==========================================================================================================

// PERF: After getting some benchmarks in place, try cutting the size of this struct in half (f64 → f32) — the 64 bits
// definitely isn't needed precision-wise, so if it helps performance, it's probably worth the small amount of type
// casting! Don't forget about padding! All fields must shrink to the same size.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct ScanKey {
    precursor: Mz,
    scan_number: usize,
}

// PERF: Shrinking this `start_time` to a `f32` won't shrink the size of this struct. If I want to avoid any packing
// anywhere, then I should move `start_time` back to `ScanInfo`, then shrink all of those fields to 4 bytes. But I'll
// need to benchmark if *increasing* the key size, whilst *decreasing* the overall K + V size actually helps
// performance. It's possible that smaller keys are faster anyways!
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct ScanValue {
    start_time: Minutes,
    peaks: Peaks,
}

impl ScanKey {
    pub fn new(precursor: f64, scan_number: usize) -> Self {
        Self {
            precursor: Mz::from(precursor),
            scan_number,
        }
    }

    pub fn ppm_window(mz: f64, ppm: f64) -> RangeInclusive<Self> {
        let (min_mz, max_mz) = Tolerance::PPM(ppm).bounds(mz);
        Self::new(min_mz, 0)..=Self::new(max_mz, usize::MAX)
    }
}

impl ScanValue {
    pub fn new(start_time: f64, peaks: Peaks) -> Self {
        Self {
            start_time: Minutes::from(start_time),
            peaks,
        }
    }
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use std::ops::{Bound, RangeBounds};

    use assert_float_eq::assert_float_absolute_eq;

    use super::*;

    #[test]
    fn scan_key_ppm_window() {
        let window = ScanKey::ppm_window(471.711_128, 10.0);
        assert!(window.contains(&ScanKey::new(471.711_128, 42)));

        let Bound::Included(&ScanKey {
            precursor: start_mz,
            scan_number: start_scan,
        }) = window.start_bound()
        else {
            panic!("expected inclusive start bound");
        };
        assert_float_absolute_eq!(start_mz.into(), 471.706_411);
        assert_eq!(start_scan, 0);

        let Bound::Included(&ScanKey {
            precursor: end_mz,
            scan_number: end_scan,
        }) = window.end_bound()
        else {
            panic!("expected inclusive end bound");
        };
        assert_float_absolute_eq!(end_mz.into(), 471.715_845);
        assert_eq!(end_scan, usize::MAX);
    }
}
