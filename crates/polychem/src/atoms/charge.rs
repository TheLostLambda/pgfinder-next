use crate::Charge;

impl Charge {
    pub(crate) fn abs(self) -> Self {
        Self(self.0.abs())
    }
}
