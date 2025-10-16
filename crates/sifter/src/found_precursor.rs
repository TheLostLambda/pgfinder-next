// Local Crate Imports
use crate::FoundPrecursor;

// Public API ==========================================================================================================

impl<'p> FoundPrecursor<'p, '_> {
    #[must_use]
    pub fn theoretical_name(&self) -> &'p str {
        self.theoretical.name()
    }

    #[must_use]
    pub fn theoretical_mz(&self) -> f64 {
        self.theoretical.mz()
    }

    // MISSING: I don't want to promise that this method is `const` in my API...
    #[expect(clippy::missing_const_for_fn)]
    #[must_use]
    pub fn observed_mz(&self) -> f64 {
        self.observed_mz
    }

    // MISSING: I don't want to promise that this method is `const` in my API...
    #[expect(clippy::missing_const_for_fn)]
    #[must_use]
    pub fn scan_number(&self) -> usize {
        self.scan_number
    }

    // MISSING: I don't want to promise that this method is `const` in my API...
    #[expect(clippy::missing_const_for_fn)]
    #[must_use]
    pub fn start_time(&self) -> f64 {
        self.start_time
    }
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use assert_float_eq::assert_float_absolute_eq;

    use crate::NamedIon;

    use super::*;

    #[test]
    fn found_precursor_getters() {
        let theoretical = NamedIon::new("gm-AEJA", 942.41);
        let ff = FoundPrecursor::new(&theoretical, 942.4121, 44, 12.345);

        assert_eq!(ff.theoretical_name(), "gm-AEJA");
        assert_float_absolute_eq!(ff.theoretical_mz(), 942.41);
        assert_float_absolute_eq!(ff.observed_mz(), 942.4121);
        assert_eq!(ff.scan_number(), 44);
        assert_float_absolute_eq!(ff.start_time(), 12.345);
    }
}
