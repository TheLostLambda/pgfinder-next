use nom::{
    branch::alt,
    character::complete::char,
    combinator::{map, opt},
    sequence::{pair, terminated, tuple},
    Parser,
};
use nom_miette::{into, map_res, wrap_err, FromExternalError, LabeledErrorKind, LabeledParseError};

use crate::{
    polymerizer::Polymerizer, AnyModification, AtomicDatabase, Count, Modification, NamedMod,
    OffsetMod, PolychemError, PolymerDatabase,
};

use super::{
    chemical_composition::chemical_composition,
    count,
    errors::{ParseResult, PolychemErrorKind},
    offset_kind,
};

// FIXME: The errors for these parsers need to be tested and improved!
// FIXME: Ensure that users of these parsers *don't* need to use nom-miette!

/// Any Modification = Named Modification | Offset Modification
pub fn any<'a, 'p, 's, K>(
    polymerizer: &Polymerizer<'a, 'p>,
    identifier: impl Parser<&'s str, &'s str, LabeledParseError<'s, K>>,
) -> impl FnMut(&'s str) -> ParseResult<AnyModification<'a, 'p>, K>
where
    K: LabeledErrorKind + From<PolychemErrorKind> + FromExternalError<'s, Box<PolychemError>>,
{
    alt((
        into(named(polymerizer.polymer_db(), identifier)),
        into(offset::<K>(polymerizer.atomic_db())),
    ))
}

// FIXME: I probably need to add a lot of `wrap_err`s around these parsers!
/// Named Modification = [ Multiplier ] , Identifier
pub fn named<'a, 'p, 's, K>(
    db: &'p PolymerDatabase<'a>,
    identifier: impl Parser<&'s str, &'s str, LabeledParseError<'s, K>>,
) -> impl FnMut(&'s str) -> ParseResult<Modification<NamedMod<'a, 'p>>, K>
where
    K: LabeledErrorKind + From<PolychemErrorKind> + FromExternalError<'s, Box<PolychemError>>,
{
    let named_mod = map_res(identifier, |abbr| NamedMod::new(db, abbr));
    let parser = pair(opt(into(multiplier)), named_mod);
    map(parser, |(multiplier, named_mod)| {
        Modification::new(multiplier.unwrap_or(1), named_mod)
    })
}

/// Offset Modification = Offset Kind , [ Multiplier ] ,
///   Chemical Composition ;
pub fn offset<'a, 's, K: From<PolychemErrorKind> + LabeledErrorKind>(
    db: &'a AtomicDatabase,
) -> impl FnMut(&'s str) -> ParseResult<Modification<OffsetMod<'a>>, K> {
    let parser = tuple((offset_kind, opt(multiplier), chemical_composition(db)));
    // FIXME: Remove this `into` combinator after updating all sub-parsers to be generic! Well, maybe I only make the
    // public parsers generic...
    into(map(parser, |(kind, multiplier, composition)| {
        Modification::new(
            multiplier.unwrap_or(1),
            OffsetMod::new_with_composition(kind, composition),
        )
    }))
}

/// Multiplier = Count , "x" ;
fn multiplier(i: &str) -> ParseResult<Count> {
    let parser = terminated(count, char('x'));
    wrap_err(parser, PolychemErrorKind::ExpectedMultiplier)(i)
}

#[cfg(test)]
mod tests {
    use nom::{
        bytes::complete::tag,
        character::complete::{alpha1, alphanumeric1},
        combinator::recognize,
        multi::many0,
    };
    use once_cell::sync::Lazy;
    use rust_decimal_macros::dec;

    use super::*;
    use crate::{Charged, Massive};

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);

    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "polymer_database.kdl",
            include_str!("../../tests/data/polymer_database.kdl"),
        )
        .unwrap()
    });

    /// Identifier = letter , { letter | digit | "_" } ;
    fn identifier(i: &str) -> ParseResult<&str> {
        recognize(pair(alpha1, many0(alt((alphanumeric1, tag("_"))))))(i)
    }

    #[test]
    fn test_multiplier() {
        // Valid Multipliers
        assert_eq!(multiplier("1x"), Ok(("", 1)));
        assert_eq!(multiplier("10x"), Ok(("", 10)));
        assert_eq!(multiplier("422x"), Ok(("", 422)));
        assert_eq!(multiplier("9999x"), Ok(("", 9999)));
        // Invalid Multipliers
        assert!(multiplier("1").is_err());
        assert!(multiplier("10").is_err());
        assert!(multiplier("422").is_err());
        assert!(multiplier("9999").is_err());
        assert!(multiplier("0").is_err());
        assert!(multiplier("01").is_err());
        assert!(multiplier("00145").is_err());
        assert!(multiplier("H").is_err());
        assert!(multiplier("p").is_err());
        assert!(multiplier("+H").is_err());
        assert!(multiplier("[H]").is_err());
        // Multiple Multipliers
        assert_eq!(multiplier("1xOH"), Ok(("OH", 1)));
        assert_eq!(multiplier("42xHeH"), Ok(("HeH", 42)));
    }

    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_named_modification() {
        let mut named_modification = named(&POLYMER_DB, identifier);
        macro_rules! assert_offset_mass {
            ($input:literal, $output:literal, $mass:expr) => {
                let (rest, modification) = named_modification($input).unwrap();
                assert_eq!(rest, $output);
                assert_eq!(modification.monoisotopic_mass(), $mass);
            };
        }
        // Valid Named Modifications
        assert_offset_mass!("Am", "", dec!(-0.98401558291));
        assert_offset_mass!("Ac", "", dec!(42.01056468403));
        assert_offset_mass!("Poly", "", dec!(77.95068082490));
        assert_offset_mass!("DeAc", "", dec!(-42.01056468403));
        assert_offset_mass!("Red", "", dec!(2.01565006446));
        assert_offset_mass!("Anh", "", dec!(-18.01056468403));
        assert_offset_mass!("1xAm", "", dec!(-0.98401558291));
        assert_offset_mass!("2xRed", "", dec!(4.03130012892));
        assert_offset_mass!("3xAnh", "", dec!(-54.03169405209));
        // Invalid Named Modifications
        assert!(named_modification(" H2O").is_err());
        assert!(named_modification("1").is_err());
        assert!(named_modification("9999").is_err());
        assert!(named_modification("0").is_err());
        assert!(named_modification("00145").is_err());
        assert!(named_modification("+H").is_err());
        assert!(named_modification("[H]").is_err());
        assert!(named_modification("Øof").is_err());
        assert!(named_modification("-Ac").is_err());
        assert!(named_modification("_Ac").is_err());
        assert!(named_modification("+Am").is_err());
        assert!(named_modification("-2xAm").is_err());
        assert!(named_modification("(Am)").is_err());
        assert!(named_modification("-4xH2O").is_err());
        assert!(named_modification("-2p").is_err());
        assert!(named_modification("+C2H2O-2e").is_err());
        assert!(named_modification("-3xC2H2O-2e").is_err());
        assert!(named_modification("+NH3+p").is_err());
        assert!(named_modification("+2xD2O").is_err());
        assert!(named_modification("-2x[2H]2O").is_err());
        // Non-Existent Named Modifications
        assert!(named_modification("Blue").is_err());
        assert!(named_modification("Hydro").is_err());
        assert!(named_modification("1xAm2").is_err());
        assert!(named_modification("2xR_ed").is_err());
        // Multiple Named Modifications
        assert_offset_mass!("Anh, Am", ", Am", dec!(-18.01056468403));
        assert_offset_mass!("1xAm)JAA", ")JAA", dec!(-0.98401558291));
    }

    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_offset_modification() {
        let mut offset_modification = offset::<PolychemErrorKind>(&ATOMIC_DB);
        macro_rules! assert_offset_mz {
            ($input:literal, $output:literal, $mass:expr, $charge:literal) => {
                let (rest, modification) = offset_modification($input).unwrap();
                assert_eq!(rest, $output);
                assert_eq!(modification.monoisotopic_mass(), $mass);
                assert_eq!(modification.charge(), $charge);
            };
        }
        // Valid Offset Modifications
        assert_offset_mz!("+H2O", "", dec!(18.01056468403), 0);
        assert_offset_mz!("-H2O", "", dec!(-18.01056468403), 0);
        assert_offset_mz!("+2xH2O", "", dec!(36.02112936806), 0);
        assert_offset_mz!("-4xH2O", "", dec!(-72.04225873612), 0);
        assert_offset_mz!("-2p", "", dec!(-2.014552933242), -2);
        assert_offset_mz!("+H", "", dec!(1.00782503223), 0);
        assert_offset_mz!("+C2H2O-2e", "", dec!(42.009467524211870), 2);
        assert_offset_mz!("-3xC2H2O-2e", "", dec!(-126.02840257263561), -6);
        assert_offset_mz!("+NH3+p", "", dec!(18.033825567741), 1);
        assert_offset_mz!("+2xD2O", "", dec!(40.04623635162), 0);
        assert_offset_mz!("-2x[2H]2O", "", dec!(-40.04623635162), 0);
        assert_offset_mz!("+[37Cl]5-2p", "", dec!(182.814960076758), -2);
        assert_offset_mz!("-NH2[99Tc]", "", dec!(-114.92497486889), 0);
        // Invalid Offset Modifications
        assert!(offset_modification(" ").is_err());
        assert!(offset_modification("H2O").is_err());
        assert!(offset_modification("(-H2O)").is_err());
        assert!(offset_modification("+0xH2O").is_err());
        assert!(offset_modification("2xH2O").is_err());
        assert!(offset_modification("-2x3xH2O").is_err());
        assert!(offset_modification("-2x+H2O").is_err());
        assert!(offset_modification("+2[2H]").is_err());
        assert!(offset_modification("-[H+p]O").is_err());
        assert!(offset_modification("+NH2[100Tc]").is_err());
        // Multiple Offset Modifications
        assert_offset_mz!("+[37Cl]5-2p10", "10", dec!(182.814960076758), -2);
        assert_offset_mz!("+[2H]2O*H2O", "*H2O", dec!(20.02311817581), 0);
        assert_offset_mz!("+NH2{100Tc", "{100Tc", dec!(16.01872406889), 0);
        assert_offset_mz!("+C11H12N2O2 H2O", " H2O", dec!(204.08987763476), 0);
    }
}
