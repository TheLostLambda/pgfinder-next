use crate::Particle;

use super::{atomic_database::AtomicDatabase, AtomicLookupError};

impl Particle {
    pub(super) fn new(
        db: &AtomicDatabase,
        symbol: impl AsRef<str>,
    ) -> std::result::Result<Self, AtomicLookupError> {
        let symbol = symbol.as_ref();
        db.particles
            .get(symbol)
            .cloned()
            .ok_or_else(|| AtomicLookupError::Particle(symbol.to_owned()))
    }
}

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
