use rust_decimal::Decimal;

use crate::{Charge, Charged, Massive, Mz, Particle};

use super::{
    atomic_database::{AtomicDatabase, ParticleDescription},
    AtomicLookupError,
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

impl Mz for Particle<'_> {}

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use miette::Result;
    use once_cell::sync::Lazy;

    use crate::testing_tools::assert_miette_snapshot;

    use super::{AtomicDatabase, Particle};

    static DB: Lazy<AtomicDatabase> = Lazy::new(|| {
        AtomicDatabase::from_kdl(
            "atomic_database.kdl",
            include_str!("../../atomic_database.kdl"),
        )
        .unwrap()
    });

    #[test]
    fn new_particle() -> Result<()> {
        // Sucessfully lookup particles that exist
        assert_debug_snapshot!(Particle::new(&DB, "p")?);
        assert_debug_snapshot!(Particle::new(&DB, "e")?);
        // Fail to lookup particles that don't exist
        assert_miette_snapshot!(Particle::new(&DB, "m"));
        Ok(())
    }
}
