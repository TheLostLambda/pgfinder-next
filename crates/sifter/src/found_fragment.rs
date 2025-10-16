// Local Crate Imports
use crate::FoundFragment;

// Public API ==========================================================================================================

impl FoundFragment<'_> {
    #[must_use]
    pub fn theoretical_precursor_name(&self) -> &str {
        self.theoretical_precursor.name()
    }

    #[must_use]
    pub fn theoretical_fragment_name(&self) -> &str {
        self.theoretical_fragment.name()
    }

    #[must_use]
    pub fn theoretical_precursor_mz(&self) -> f64 {
        self.theoretical_precursor.mz()
    }

    #[must_use]
    pub fn theoretical_fragment_mz(&self) -> f64 {
        self.theoretical_fragment.mz()
    }

    // MISSING: I don't want to promise that this method is `const` in my API...
    #[expect(clippy::missing_const_for_fn)]
    #[must_use]
    pub fn observed_precursor_mz(&self) -> f64 {
        self.observed_precursor_mz
    }

    // MISSING: I don't want to promise that this method is `const` in my API...
    #[expect(clippy::missing_const_for_fn)]
    #[must_use]
    pub fn observed_fragment_mz(&self) -> f64 {
        self.observed_fragment_mz
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
    fn found_fragment_getters() {
        let ff = FoundFragment::new(
            NamedIon::new("gm-AEJA", 942.41),
            NamedIon::new("J", 173.09),
            942.4121,
            173.0923,
            44,
            12.345,
        );

        assert_eq!(ff.theoretical_precursor_name(), "gm-AEJA");
        assert_eq!(ff.theoretical_fragment_name(), "J");
        assert_float_absolute_eq!(ff.theoretical_precursor_mz(), 942.41);
        assert_float_absolute_eq!(ff.theoretical_fragment_mz(), 173.09);
        assert_float_absolute_eq!(ff.observed_precursor_mz(), 942.4121);
        assert_float_absolute_eq!(ff.observed_fragment_mz(), 173.0923);
        assert_eq!(ff.scan_number(), 44);
        assert_float_absolute_eq!(ff.start_time(), 12.345);
    }
}
