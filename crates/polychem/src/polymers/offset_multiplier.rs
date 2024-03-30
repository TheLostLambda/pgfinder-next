use std::cmp::Ordering;

use crate::{Count, OffsetKind, OffsetMultiplier, PolychemError, SignedCount};

impl TryFrom<SignedCount> for OffsetMultiplier {
    // FIXME: Actually add a PolychemError variant!
    type Error = PolychemError;

    fn try_from(value: SignedCount) -> Result<Self, Self::Error> {
        // FIXME: Replace unwraps with ?, and fill in the todo!
        Ok(match value.cmp(&0) {
            Ordering::Less => Self(
                OffsetKind::Remove,
                Count::try_from(value.unsigned_abs()).unwrap(),
            ),
            Ordering::Equal => todo!(),
            Ordering::Greater => Self(OffsetKind::Add, Count::try_from(value).unwrap()),
        })
    }
}

impl From<OffsetMultiplier> for SignedCount {
    fn from(value: OffsetMultiplier) -> Self {
        let OffsetMultiplier(kind, count) = value;
        Self::from(kind) * Self::from(count)
    }
}
