use std::fmt::{self, Display, Formatter};

use crate::{Charge, Charged, Mass, Particle};

use super::{
    atomic_database::{AtomicDatabase, ParticleDescription},
    errors::AtomicLookupError,
};

impl<'a> Particle<'a> {
    pub(crate) fn new(
        db: &'a AtomicDatabase,
        symbol: impl AsRef<str>,
    ) -> Result<Self, AtomicLookupError> {
        let symbol = symbol.as_ref();
        let (symbol, ParticleDescription { name, mass, charge }) = db
            .particles
            .get_key_value(symbol)
            .ok_or_else(|| AtomicLookupError::particle(symbol))?;
        Ok(Self {
            symbol,
            name,
            mass,
            charge,
        })
    }

    pub(crate) const fn mass(&self) -> Mass {
        *self.mass
    }
}

impl Display for Particle<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let symbol = self.symbol;
        write!(f, "{symbol}")
    }
}

impl Charged for Particle<'_> {
    fn charge(&self) -> Charge {
        *self.charge
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use rust_decimal_macros::dec;
    use std::sync::LazyLock;

    use crate::testing_tools::assert_miette_snapshot;

    use super::*;

    static DB: LazyLock<AtomicDatabase> = LazyLock::new(AtomicDatabase::default);

    #[test]
    fn new_particle() {
        // Sucessfully lookup particles that exist
        assert_debug_snapshot!(Particle::new(&DB, "p").unwrap());
        assert_debug_snapshot!(Particle::new(&DB, "e").unwrap());
        // Fail to lookup particles that don't exist
        assert_miette_snapshot!(Particle::new(&DB, "m"));
    }

    #[test]
    fn particle_display() {
        let p = Particle::new(&DB, "p").unwrap();
        assert_eq!(p.to_string(), "p");
        let e = Particle::new(&DB, "e").unwrap();
        assert_eq!(e.to_string(), "e");
    }

    #[test]
    fn particle_masses() {
        let proton = Particle::new(&DB, "p").unwrap();
        assert_eq!(proton.mass(), Mass(dec!(1.007276466621)));
        let electron = Particle::new(&DB, "e").unwrap();
        assert_eq!(electron.mass(), Mass(dec!(0.000548579909065)));
    }

    #[test]
    fn particle_charges() {
        let proton = Particle::new(&DB, "p").unwrap();
        assert_eq!(proton.charge(), Charge(1));
        let electron = Particle::new(&DB, "e").unwrap();
        assert_eq!(electron.mass(), Mass(dec!(0.000548579909065)));
        assert_eq!(electron.charge(), Charge(-1));

        assert_eq!(proton.charge(), -electron.charge());
    }
}
