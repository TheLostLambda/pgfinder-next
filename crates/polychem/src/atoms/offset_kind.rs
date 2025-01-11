use std::fmt::{self, Display, Formatter};
use std::ops::{Neg, Not};

use crate::OffsetKind;

impl OffsetKind {
    pub fn offset<T: Neg<Output = T>>(self, value: T) -> T {
        match self {
            Self::Add => value,
            Self::Remove => -value,
        }
    }
}

impl Not for OffsetKind {
    type Output = Self;

    fn not(self) -> Self::Output {
        match self {
            Self::Add => Self::Remove,
            Self::Remove => Self::Add,
        }
    }
}

impl Display for OffsetKind {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}", match self {
            Self::Add => "+",
            Self::Remove => "-",
        })
    }
}

#[cfg(test)]
mod tests {
    use rust_decimal_macros::dec;

    use super::*;

    #[test]
    fn offset() {
        let add = OffsetKind::Add;
        let remove = OffsetKind::Remove;
        assert_eq!(add.offset(1), 1);
        assert_eq!(remove.offset(1), -1);
        assert_eq!(add.offset(dec!(3.14)), dec!(3.14));
        assert_eq!(remove.offset(dec!(3.14)), dec!(-3.14));
    }

    #[test]
    fn not() {
        let add = OffsetKind::Add;
        let remove = OffsetKind::Remove;
        assert_eq!(!add, remove);
        assert_eq!(!remove, add);
    }

    #[test]
    fn offset_kind_display() {
        let add = OffsetKind::Add;
        assert_eq!(add.to_string(), "+");
        let remove = OffsetKind::Remove;
        assert_eq!(remove.to_string(), "-");
    }
}
