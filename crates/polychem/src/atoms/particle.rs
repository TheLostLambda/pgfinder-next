use std::fmt::{self, Display, Formatter};

use rust_decimal::Decimal;

use crate::{Charge, Charged, Massive, Particle};

use super::{
    atomic_database::{AtomicDatabase, ParticleDescription},
    errors::AtomicLookupError,
};

impl<'a> Particle<'a> {
    pub(super) fn new(
        db: &'a AtomicDatabase,
        symbol: impl AsRef<str>,
    ) -> Result<Self, AtomicLookupError> {
        let symbol = symbol.as_ref();
        let (symbol, ParticleDescription { name, mass, charge }) = db
            .particles
            .get_key_value(symbol)
            .ok_or_else(|| AtomicLookupError::Particle(symbol.to_owned()))?;
        Ok(Self {
            symbol,
            name,
            mass,
            charge,
        })
    }
}

impl Display for Particle<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let symbol = self.symbol;
        write!(f, "{symbol}")
    }
}

impl Massive for Particle<'_> {
    fn monoisotopic_mass(&self) -> Decimal {
        *self.mass
    }

    fn average_mass(&self) -> Decimal {
        *self.mass
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
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use crate::testing_tools::assert_miette_snapshot;

    use super::*;

    static DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);

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
        assert_eq!(*proton.mass, dec!(1.007276466621));
        assert_eq!(proton.monoisotopic_mass(), *proton.mass);
        assert_eq!(proton.average_mass(), *proton.mass);
        let electron = Particle::new(&DB, "e").unwrap();
        assert_eq!(*electron.mass, dec!(0.000548579909065));
        assert_eq!(electron.monoisotopic_mass(), *electron.mass);
        assert_eq!(electron.average_mass(), *electron.mass);
    }

    #[test]
    fn particle_charges() {
        let proton = Particle::new(&DB, "p").unwrap();
        assert_eq!(proton.charge(), 1);
        let electron = Particle::new(&DB, "e").unwrap();
        assert_eq!(*electron.mass, dec!(0.000548579909065));
        assert_eq!(electron.charge(), -1);

        assert_eq!(proton.charge(), -electron.charge());
    }
}
