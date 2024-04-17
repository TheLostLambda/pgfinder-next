use crate::Charge;

impl Charge {
    pub(crate) fn abs(&self) -> Self {
        Self(self.0.abs())
    }

    pub(crate) fn is_zero(&self) -> bool {
        self.0 == 0
    }
}
