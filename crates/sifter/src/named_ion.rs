// Standard Library Imports
use std::borrow::Cow;

// Local Crate Imports
use crate::NamedIon;

// Public API ==========================================================================================================

impl<'n> NamedIon<'n> {
    pub fn new(name: impl Into<Cow<'n, str>>, mz: impl Into<f64>) -> Self {
        Self {
            name: name.into(),
            mz: mz.into(),
        }
    }

    #[must_use]
    pub fn name(&self) -> &str {
        self.name.as_ref()
    }

    // MISSING: I don't want to promise that this method is `const` in my API — especially when no other part of
    // `NamedIon` is...
    #[expect(clippy::missing_const_for_fn)]
    #[must_use]
    pub fn mz(&self) -> f64 {
        self.mz
    }
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use assert_float_eq::assert_float_absolute_eq;

    use super::*;

    #[test]
    fn named_ion_new() {
        let from_str = NamedIon::new("magic", 42);
        let from_string = NamedIon::new(String::from("magic"), 42.0);
        assert_eq!(from_str, from_string);
    }

    #[test]
    fn named_ion_getters() {
        let name: String = "mewo".into();
        let mz: u8 = 255;
        let named_ion = NamedIon::new(name, mz);

        assert_eq!(named_ion.name(), "mewo");
        assert_float_absolute_eq!(named_ion.mz(), 2.55e2);
    }
}
