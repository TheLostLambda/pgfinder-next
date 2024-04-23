use crate::Charge;

impl Charge {
    pub(crate) const fn abs(self) -> Self {
        Self(self.0.abs())
    }
}
