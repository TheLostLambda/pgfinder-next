use super::{Element, MassNumber, Particle};
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
use std::{borrow::BorrowMut, collections::HashMap, ops::Deref, str::FromStr};
use thiserror::Error;

// pub struct ChemicalDatabase {
//     elements: HashMap<String, Element>,
//     particles: HashMap<String, Particle>,
// }

// impl ChemicalDatabase {
//     pub fn from_kdl(kdl: impl AsRef<str>) -> miette::Result<Self> {
//         todo!()
//     }
// }

#[derive(Decode, Debug)]
struct ChemicalDatabaseKdl {
    #[knuffel(child, unwrap(children))]
    elements: Vec<ElementKdl>,
    #[knuffel(child, unwrap(children))]
    particles: Vec<ParticleKdl>,
}

#[derive(Decode, Debug)]
struct ElementKdl {
    #[knuffel(node_name)]
    symbol: ElementSymbol,
    #[knuffel(argument)]
    name: String,
    #[knuffel(children(name = "isotope"))]
    isotopes: Vec<IsotopeKdl>,
}

#[derive(Decode, Debug)]
struct ParticleKdl {
    #[knuffel(node_name)]
    symbol: ParticleSymbol,
    #[knuffel(argument)]
    name: String,
    #[knuffel(child, unwrap(argument))]
    mass: DecimalKdl,
    #[knuffel(child, unwrap(argument))]
    charge: i32,
}

#[derive(Decode, Debug)]
struct IsotopeKdl {
    #[knuffel(argument)]
    mass_number: MassNumber,
    #[knuffel(argument)]
    relative_mass: DecimalKdl,
    #[knuffel(argument)]
    abundance: Option<DecimalKdl>,
}

#[derive(Debug, Default, PartialEq, Eq)]
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
        match value.deref() {
            Literal::Decimal(ast::Decimal(s)) | Literal::Int(Integer(Radix::Dec, s)) => {
                let res = if s.contains(['e', 'E']) {
                    Decimal::from_scientific(s)
                } else {
                    Decimal::from_str_exact(s)
                };
                match res {
                    Ok(d) => Ok(Self(d)),
                    Err(e) => {
                        ctx.emit_error(DecodeError::Conversion {
                            span: value.span().clone(),
                            source: Box::new(e),
                        });
                        Ok(Default::default())
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
                Ok(Default::default())
            }
        }
    }
}

#[derive(Error, Diagnostic, PartialEq, Eq, Debug)]
enum InvalidChemicalSymbolError {
    #[error(
        "expected a single uppercase ASCII letter optionally followed by a lowercase ASCII letter, got {0:?}"
    )]
    Element(String),
    #[error("expected a single lowercase ASCII letter, got {0:?}")]
    Particle(String),
}

#[derive(Debug)]
struct ElementSymbol(String);

impl FromStr for ElementSymbol {
    type Err = InvalidChemicalSymbolError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.to_owned();
        let chrs: Vec<_> = s.chars().collect();
        match chrs[..] {
            [f] if f.is_ascii_uppercase() => Ok(Self(s)),
            [f, l] if f.is_ascii_uppercase() && l.is_ascii_lowercase() => Ok(Self(s)),
            _ => Err(InvalidChemicalSymbolError::Element(s)),
        }
    }
}

#[derive(Debug)]
struct ParticleSymbol(String);

impl FromStr for ParticleSymbol {
    type Err = InvalidChemicalSymbolError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.to_owned();
        if s.len() == 1 && s.chars().next().unwrap().is_ascii_lowercase() {
            Ok(Self(s))
        } else {
            Err(InvalidChemicalSymbolError::Particle(s))
        }
    }
}

#[cfg(test)]
mod tests {
    use std::error::Error;

    use indoc::indoc;
    use insta::assert_debug_snapshot;
    use knuffel::{self, Decode};
    use miette::{Diagnostic, Result};
    use rust_decimal_macros::dec;

    use super::{
        ChemicalDatabaseKdl, DecimalKdl, ElementKdl, InvalidChemicalSymbolError, ParticleKdl,
    };

    const KDL: &str = include_str!("chemistry.kdl");

    #[test]
    fn parse_default_chemical_database() -> Result<()> {
        let db: ChemicalDatabaseKdl = knuffel::parse("chemistry.kdl", KDL)?;
        // Basic property checking
        assert_eq!(db.elements.len(), 120); // 118 + 2 for deuterium and tritium
        assert_eq!(db.particles.len(), 2);
        assert_eq!(db.elements.iter().flat_map(|e| &e.isotopes).count(), 356);
        // Full snapshot test
        assert_debug_snapshot!(db);
        Ok(())
    }

    #[test]
    fn uppercase_particle_symbol() -> Result<()> {
        let kdl = indoc! {r#"
            P "Proton" {
              mass 1.007276466621
              charge +1
            }
        "#};
        let res = knuffel::parse::<Vec<ParticleKdl>>("test", kdl);
        assert!(res.is_err());
        assert_eq!(
            extract_err::<InvalidChemicalSymbolError>(&res.unwrap_err()),
            &InvalidChemicalSymbolError::Particle("P".to_string())
        );
        Ok(())
    }

    #[test]
    fn multiple_char_particle_symbol() -> Result<()> {
        let kdl = indoc! {r#"
            pr "Proton" {
              mass 1.007276466621
              charge +1
            }
        "#};
        let res = knuffel::parse::<Vec<ParticleKdl>>("test", kdl);
        assert!(res.is_err());
        assert_eq!(
            extract_err::<InvalidChemicalSymbolError>(&res.unwrap_err()),
            &InvalidChemicalSymbolError::Particle("pr".to_string())
        );
        Ok(())
    }

    #[test]
    fn lowercase_element_symbol() -> Result<()> {
        let kdl = indoc! {r#"
            d "Deuterium" {
              isotope 2 2.01410177812 1
            }
        "#};
        let res = knuffel::parse::<Vec<ElementKdl>>("test", kdl);
        assert!(res.is_err());
        assert_eq!(
            extract_err::<InvalidChemicalSymbolError>(&res.unwrap_err()),
            &InvalidChemicalSymbolError::Element("d".to_string())
        );
        Ok(())
    }

    #[test]
    fn double_uppercase_element_symbol() -> Result<()> {
        let kdl = indoc! {r#"
            DT "Deuterium" {
              isotope 2 2.01410177812 1
            }
        "#};
        let res = knuffel::parse::<Vec<ElementKdl>>("test", kdl);
        assert!(res.is_err());
        assert_eq!(
            extract_err::<InvalidChemicalSymbolError>(&res.unwrap_err()),
            &InvalidChemicalSymbolError::Element("DT".to_string())
        );
        Ok(())
    }

    #[test]
    fn three_char_element_symbol() -> Result<()> {
        let kdl = indoc! {r#"
            Deu "Deuterium" {
              isotope 2 2.01410177812 1
            }
        "#};
        let res = knuffel::parse::<Vec<ElementKdl>>("test", kdl);
        assert!(res.is_err());
        assert_eq!(
            extract_err::<InvalidChemicalSymbolError>(&res.unwrap_err()),
            &InvalidChemicalSymbolError::Element("Deu".to_string())
        );
        Ok(())
    }

    #[derive(Decode, Debug, PartialEq, Eq)]
    struct Lossless(#[knuffel(argument)] DecimalKdl);

    #[test]
    fn decimal_underflow() -> Result<()> {
        let kdl = "lossless 0.00000_00000_00000_00000_00000_0001";
        let res = knuffel::parse::<Vec<Lossless>>("test", kdl);
        assert!(res.is_err());
        assert_eq!(
            extract_err::<rust_decimal::Error>(&res.unwrap_err()),
            &rust_decimal::Error::Underflow
        );
        Ok(())
    }

    #[test]
    fn decimal_scientific() -> Result<()> {
        let kdl = "lossless 5.485_799_090_65e-4";
        let res = knuffel::parse::<Vec<Lossless>>("test", kdl);
        assert!(&res.is_ok());
        assert_eq!(res.unwrap()[0].0 .0, dec!(0.000548579909065));
        Ok(())
    }

    #[test]
    fn decimal_from_integer() -> Result<()> {
        let kdl = "lossless 1";
        let res = knuffel::parse::<Vec<Lossless>>("test", kdl);
        assert!(&res.is_ok());
        assert_eq!(res.unwrap()[0].0 .0, dec!(1));
        Ok(())
    }

    #[test]
    fn decimal_lack_of_precision() -> Result<()> {
        let kdl = "lossless 1e-42";
        let res = knuffel::parse::<Vec<Lossless>>("test", kdl);
        assert!(res.is_err());
        assert_eq!(
            extract_err::<rust_decimal::Error>(&res.unwrap_err()),
            &rust_decimal::Error::ScaleExceedsMaximumPrecision(42)
        );
        Ok(())
    }

    #[test]
    fn decimal_illegal_type() -> Result<()> {
        let kdl = "lossless (pi)3.14";
        let res = knuffel::parse::<Vec<Lossless>>("test", kdl);
        assert!(res.is_err());
        assert_debug_snapshot!(&res.unwrap_err().related().unwrap().next().unwrap());
        Ok(())
    }

    #[test]
    fn decimal_from_bool() -> Result<()> {
        let kdl = "lossless true";
        let res = knuffel::parse::<Vec<Lossless>>("test", kdl);
        assert!(res.is_err());
        assert_debug_snapshot!(&res.unwrap_err().related().unwrap().next().unwrap());
        Ok(())
    }

    fn extract_err<T: Error + 'static>(e: &knuffel::Error) -> &T {
        e.related()
            .unwrap()
            .next()
            .unwrap()
            .source()
            .unwrap()
            .downcast_ref::<T>()
            .unwrap()
    }
}
