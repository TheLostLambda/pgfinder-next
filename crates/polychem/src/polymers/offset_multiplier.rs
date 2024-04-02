use std::cmp::Ordering;

use crate::{Count, OffsetKind, OffsetMultiplier, SignedCount};

use super::errors::OffsetMultiplierError;

impl TryFrom<SignedCount> for OffsetMultiplier {
    type Error = OffsetMultiplierError;

    fn try_from(value: SignedCount) -> Result<Self, Self::Error> {
        let abs_value = value.unsigned_abs();
        let count =
            Count::try_from(abs_value).map_err(|_| OffsetMultiplierError::TooLarge(abs_value))?;
        match value.cmp(&0) {
            Ordering::Less => Ok(Self(OffsetKind::Remove, count)),
            Ordering::Equal => Err(OffsetMultiplierError::Zero),
            Ordering::Greater => Ok(Self(OffsetKind::Add, count)),
        }
    }
}

impl From<OffsetMultiplier> for SignedCount {
    fn from(value: OffsetMultiplier) -> Self {
        let OffsetMultiplier(kind, count) = value;
        Self::from(kind) * Self::from(count)
    }
}

#[cfg(test)]
mod tests {
    use crate::testing_tools::assert_miette_snapshot;

    use super::*;

    #[test]
    fn signed_count_to_offset_multiplier() {
        let negative: SignedCount = -3;
        let zero: SignedCount = 0;
        let positive: SignedCount = 2;
        assert_eq!(
            negative.try_into(),
            Ok(OffsetMultiplier(OffsetKind::Remove, 3))
        );
        assert_miette_snapshot!(OffsetMultiplier::try_from(zero));
        assert_eq!(
            positive.try_into(),
            Ok(OffsetMultiplier(OffsetKind::Add, 2))
        );

        let too_big = SignedCount::MAX;
        let too_small = SignedCount::MIN;
        assert_miette_snapshot!(OffsetMultiplier::try_from(too_big));
        assert_miette_snapshot!(OffsetMultiplier::try_from(too_small));
    }

    #[test]
    fn offset_multiplier_to_signed_count() {
        let max = OffsetMultiplier(OffsetKind::Add, Count::MAX);
        let min = OffsetMultiplier(OffsetKind::Remove, Count::MAX);
        assert_eq!(SignedCount::from(max), 4_294_967_295);
        assert_eq!(SignedCount::from(min), -4_294_967_295);
    }
}
