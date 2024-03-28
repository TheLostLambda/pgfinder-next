use std::fmt::{self, Display, Formatter};

use rust_decimal::Decimal;

use crate::{Charge, OffsetKind};

impl Display for OffsetKind {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Add => "+",
                Self::Remove => "-",
            }
        )
    }
}

impl From<OffsetKind> for Decimal {
    fn from(value: OffsetKind) -> Self {
        Charge::from(value).into()
    }
}

impl From<OffsetKind> for Charge {
    fn from(value: OffsetKind) -> Self {
        match value {
            OffsetKind::Add => 1,
            OffsetKind::Remove => -1,
        }
    }
}

#[cfg(test)]
mod tests {
    use rust_decimal_macros::dec;

    use super::*;

    #[test]
    fn offset_kind_display() {
        let add = OffsetKind::Add;
        assert_eq!(add.to_string(), "+");
        let remove = OffsetKind::Remove;
        assert_eq!(remove.to_string(), "-");
    }

    #[test]
    fn into_charge() {
        let add = OffsetKind::Add;
        assert_eq!(Charge::from(add), 1);
        let remove = OffsetKind::Remove;
        assert_eq!(Charge::from(remove), -1);
    }

    #[test]
    fn into_decimal() {
        let add = OffsetKind::Add;
        assert_eq!(Decimal::from(add), dec!(1));
        let remove = OffsetKind::Remove;
        assert_eq!(Decimal::from(remove), dec!(-1));
    }
}
