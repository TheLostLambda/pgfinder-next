use std::fmt::{Display, Formatter};

use crate::FunctionalGroup;

impl FunctionalGroup {
    pub fn new(name: impl Into<String>, location: impl Into<String>) -> Self {
        let name = name.into();
        let location = location.into();
        Self { name, location }
    }
}

impl Display for FunctionalGroup {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?} at={:?}", self.name, self.location)
    }
}

#[cfg(test)]
mod tests {
    use crate::FunctionalGroup;

    #[test]
    fn display() {
        let n_terminal = FunctionalGroup::new("Amino", "N-Terminal");
        assert_eq!(n_terminal.to_string(), r#""Amino" at="N-Terminal""#);
        let c_terminal = FunctionalGroup::new("Carboxyl", "C-Terminal");
        assert_eq!(c_terminal.to_string(), r#""Carboxyl" at="C-Terminal""#);
    }
}
