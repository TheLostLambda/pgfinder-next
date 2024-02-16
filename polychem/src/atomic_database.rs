// Standard Library Imports
use std::{collections::HashMap, ops::Deref, str::FromStr};

// External Crate Imports
use knuffel::{
    ast::{self, Integer, Literal, Radix, TypeName},
    decode::{Context, Kind},
    errors::{DecodeError, ExpectedType},
    span::Spanned,
    traits::ErrorSpan,
    Decode, DecodeScalar,
};
use miette::{Diagnostic, Result};
use rust_decimal::Decimal;
use thiserror::Error;

// Local Module Imports
use super::{Charge, Element, Isotope, MassNumber, Particle};

// Public API ==========================================================================================================

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct AtomicDatabase {
    pub(super) elements: HashMap<String, Element>,
    pub(super) particles: HashMap<String, Particle>,
}

impl AtomicDatabase {
    pub fn from_kdl(file_name: impl AsRef<str>, text: impl AsRef<str>) -> Result<Self> {
        let parsed_db: AtomicDatabaseKdl = knuffel::parse(file_name.as_ref(), text.as_ref())?;
        let elements = parsed_db
            .elements
            .into_iter()
            .map(ElementEntry::from)
            .collect();
        let particles = parsed_db
            .particles
            .into_iter()
            .map(ParticleEntry::from)
            .collect();
        Ok(Self {
            elements,
            particles,
        })
    }
}

// KDL File Schema =====================================================================================================

#[derive(Debug, Decode)]
struct AtomicDatabaseKdl {
    #[knuffel(child, unwrap(children))]
    elements: Vec<ElementKdl>,
    #[knuffel(child, unwrap(children))]
    particles: Vec<ParticleKdl>,
}

#[derive(Debug, Decode)]
struct ElementKdl {
    #[knuffel(node_name)]
    symbol: ElementSymbol,
    #[knuffel(argument)]
    name: String,
    #[knuffel(children(name = "isotope", non_empty))]
    isotopes: Vec<IsotopeKdl>,
}

#[derive(Debug, Decode)]
struct ParticleKdl {
    #[knuffel(node_name)]
    symbol: ParticleSymbol,
    #[knuffel(argument)]
    name: String,
    #[knuffel(child, unwrap(argument))]
    mass: DecimalKdl,
    #[knuffel(child, unwrap(argument))]
    charge: Charge,
}

#[derive(Debug, Decode)]
struct IsotopeKdl {
    #[knuffel(argument)]
    mass_number: MassNumber,
    #[knuffel(argument)]
    relative_mass: DecimalKdl,
    #[knuffel(argument)]
    abundance: Option<DecimalKdl>,
}

// Lossless Parsing of KDL Numbers to Decimal ==========================================================================

#[derive(Debug, Default)]
struct DecimalKdl(Decimal);

impl<S: ErrorSpan> DecodeScalar<S> for DecimalKdl {
    fn type_check(type_name: &Option<Spanned<TypeName, S>>, ctx: &mut Context<S>) {
        if let Some(t) = type_name {
            ctx.emit_error(DecodeError::TypeName {
                span: t.span().clone(),
                found: Some(t.deref().clone()),
                expected: ExpectedType::no_type(),
                rust_type: "Decimal",
            });
        }
    }

    fn raw_decode(
        value: &Spanned<Literal, S>,
        ctx: &mut Context<S>,
    ) -> Result<Self, DecodeError<S>> {
        match &**value {
            Literal::Decimal(ast::Decimal(s)) | Literal::Int(Integer(Radix::Dec, s)) => {
                let res = if s.contains(['e', 'E']) {
                    Decimal::from_scientific(s)
                } else {
                    Decimal::from_str_exact(s)
                };
                match res {
                    Ok(d) => Ok(Self(d)),
                    Err(e) => {
                        ctx.emit_error(DecodeError::conversion(value, Box::new(e)));
                        Ok(Self::default())
                    }
                }
            }
            unsupported => {
                ctx.emit_error(DecodeError::unsupported(
                    value,
                    format!(
                        "expected a decimal number, found {}",
                        Kind::from(unsupported)
                    ),
                ));
                Ok(Self::default())
            }
        }
    }
}

// Element and Particle Symbol Validation ==============================================================================

#[derive(Debug)]
struct ElementSymbol(String);

impl FromStr for ElementSymbol {
    type Err = InvalidAtomicSymbolError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.to_owned();
        let chrs: Vec<_> = s.chars().collect();
        match chrs[..] {
            [f] if f.is_ascii_uppercase() => Ok(Self(s)),
            [f, l] if f.is_ascii_uppercase() && l.is_ascii_lowercase() => Ok(Self(s)),
            _ => Err(InvalidAtomicSymbolError::Element(s)),
        }
    }
}

#[derive(Debug)]
struct ParticleSymbol(String);

impl FromStr for ParticleSymbol {
    type Err = InvalidAtomicSymbolError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.to_owned();
        if s.len() == 1 && s.chars().next().unwrap().is_ascii_lowercase() {
            Ok(Self(s))
        } else {
            Err(InvalidAtomicSymbolError::Particle(s))
        }
    }
}

#[derive(Eq, PartialEq, Debug, Diagnostic, Error)]
enum InvalidAtomicSymbolError {
    #[error(
        "expected a single uppercase ASCII letter optionally followed by a lowercase ASCII letter, got {0:?}"
    )]
    Element(String),
    #[error("expected a single lowercase ASCII letter, got {0:?}")]
    Particle(String),
}

// Conversion From Parsed KDL to Internal Representation ===============================================================

type ElementEntry = (String, Element);

impl From<ElementKdl> for ElementEntry {
    fn from(
        ElementKdl {
            symbol,
            name,
            isotopes,
        }: ElementKdl,
    ) -> Self {
        let isotopes = isotopes.into_iter().map(IsotopeEntry::from).collect();
        (
            symbol.0.clone(),
            Element {
                symbol: symbol.0,
                name,
                mass_number: None,
                isotopes,
            },
        )
    }
}

type IsotopeEntry = (MassNumber, Isotope);

impl From<IsotopeKdl> for IsotopeEntry {
    fn from(
        IsotopeKdl {
            mass_number,
            relative_mass,
            abundance,
        }: IsotopeKdl,
    ) -> Self {
        (
            mass_number,
            Isotope {
                relative_mass: relative_mass.0,
                abundance: abundance.map(|a| a.0),
            },
        )
    }
}

type ParticleEntry = (String, Particle);

impl From<ParticleKdl> for ParticleEntry {
    fn from(
        ParticleKdl {
            symbol,
            name,
            mass,
            charge,
        }: ParticleKdl,
    ) -> Self {
        (
            symbol.0.clone(),
            Particle {
                symbol: symbol.0,
                name,
                mass: mass.0,
                charge,
            },
        )
    }
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use indoc::indoc;
    use insta::assert_debug_snapshot;
    use knuffel::{self, Decode};
    use miette::Result;
    use rust_decimal_macros::dec;

    use super::{AtomicDatabase, AtomicDatabaseKdl, DecimalKdl, ElementKdl, ParticleKdl};

    use crate::testing_tools::assert_miette_snapshot;

    const KDL: &str = include_str!("../atomic_database.kdl");

    #[test]
    fn parse_default_atomic_database() {
        let db: AtomicDatabaseKdl = knuffel::parse("atomic_database.kdl", KDL).unwrap();
        // Basic property checking
        assert_eq!(db.elements.len(), 120); // 118 + 2 for deuterium and tritium
        assert_eq!(db.particles.len(), 2);
        assert_eq!(db.elements.iter().flat_map(|e| &e.isotopes).count(), 356);
        // Full snapshot test
        assert_debug_snapshot!(db);
    }

    #[test]
    fn build_default_atomic_database() {
        let db = AtomicDatabase::from_kdl("atomic_database.kdl", KDL).unwrap();
        // Basic property checking
        assert_eq!(db.elements.len(), 120); // 118 + 2 for deuterium and tritium
        assert_eq!(db.particles.len(), 2);
        assert_eq!(
            db.elements.iter().flat_map(|(_, e)| &e.isotopes).count(),
            356
        );
    }

    #[test]
    fn uppercase_particle_symbol() {
        let kdl = indoc! {r#"
            P "Proton" {
              mass 1.007276466621
              charge +1
            }
        "#};
        let res = knuffel::parse::<Vec<ParticleKdl>>("test", kdl);
        assert!(res.is_err());
        assert_miette_snapshot!(res);
    }

    #[test]
    fn multiple_char_particle_symbol() {
        let kdl = indoc! {r#"
            pr "Proton" {
              mass 1.007276466621
              charge +1
            }
        "#};
        let res = knuffel::parse::<Vec<ParticleKdl>>("test", kdl);
        assert!(res.is_err());
        assert_miette_snapshot!(res);
    }

    #[test]
    fn lowercase_element_symbol() {
        let kdl = indoc! {r#"
            d "Deuterium" {
              isotope 2 2.01410177812 1
            }
        "#};
        let res = knuffel::parse::<Vec<ElementKdl>>("test", kdl);
        assert!(res.is_err());
        assert_miette_snapshot!(res);
    }

    #[test]
    fn double_uppercase_element_symbol() {
        let kdl = indoc! {r#"
            DT "Deuterium" {
              isotope 2 2.01410177812 1
            }
        "#};
        let res = knuffel::parse::<Vec<ElementKdl>>("test", kdl);
        assert!(res.is_err());
        assert_miette_snapshot!(res);
    }

    #[test]
    fn three_char_element_symbol() {
        let kdl = indoc! {r#"
            Deu "Deuterium" {
              isotope 2 2.01410177812 1
            }
        "#};
        let res = knuffel::parse::<Vec<ElementKdl>>("test", kdl);
        assert!(res.is_err());
        assert_miette_snapshot!(res);
    }

    #[test]
    fn element_without_isotopes() {
        let kdl = indoc! {r#"
            D "Deuterium" {
              // isotope 2 2.01410177812 1
            }
        "#};
        let res = knuffel::parse::<Vec<ElementKdl>>("test", kdl);
        assert!(res.is_err());
        assert_miette_snapshot!(res);
    }

    #[derive(Debug, Decode)]
    struct Lossless(#[knuffel(argument)] DecimalKdl);

    #[test]
    fn decimal_underflow() {
        let kdl = "lossless 0.00000_00000_00000_00000_00000_0001";
        let res = knuffel::parse::<Vec<Lossless>>("test", kdl);
        assert!(res.is_err());
        assert_miette_snapshot!(res);
    }

    #[test]
    fn decimal_scientific() {
        let kdl = "lossless 5.485_799_090_65e-4";
        let res = knuffel::parse::<Vec<Lossless>>("test", kdl);
        assert!(&res.is_ok());
        assert_eq!(res.unwrap()[0].0 .0, dec!(0.000548579909065));
    }

    #[test]
    fn decimal_from_integer() {
        let kdl = "lossless 1";
        let res = knuffel::parse::<Vec<Lossless>>("test", kdl);
        assert!(&res.is_ok());
        assert_eq!(res.unwrap()[0].0 .0, dec!(1));
    }

    #[test]
    fn decimal_lack_of_precision() {
        let kdl = "lossless 1e-42";
        let res = knuffel::parse::<Vec<Lossless>>("test", kdl);
        assert!(res.is_err());
        assert_miette_snapshot!(res);
    }

    #[test]
    fn decimal_illegal_type() {
        let kdl = "lossless (pi)3.14";
        let res = knuffel::parse::<Vec<Lossless>>("test", kdl);
        assert!(res.is_err());
        assert_miette_snapshot!(res);
    }

    #[test]
    fn decimal_from_bool() {
        let kdl = "lossless true";
        let res = knuffel::parse::<Vec<Lossless>>("test", kdl);
        assert!(res.is_err());
        assert_miette_snapshot!(res);
    }
}
