// Standard Library Imports
use std::ops::RangeInclusive;

// Public API ==========================================================================================================

pub trait PpmWindow: Sized {
    fn ppm_window(mz: f64, ppm: f64) -> RangeInclusive<Self>;
}
