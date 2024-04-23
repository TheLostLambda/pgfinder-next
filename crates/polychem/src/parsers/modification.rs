use nom::{branch::alt, character::complete::char, sequence::terminated, Parser};
use nom_miette::{wrap_err, LabeledParseError};

use crate::{Count, ModificationId, Polymer};

use super::{
    count,
    errors::{ParseResult, PolychemErrorKind},
};

// FIXME: The errors for these parsers need to be tested and improved!
// FIXME: Ensure that users of these parsers *don't* need to use nom-miette!

/// Any Modification = Named Modification | Offset Modification
pub fn any<'a, 'p, 's, K>(
    polymer: &mut Polymer<'a, 'p>,
    identifier: impl Parser<&'s str, &'s str, LabeledParseError<'s, K>>,
) -> impl FnMut(&'s str) -> ParseResult<ModificationId, K>
where
    K: From<PolychemErrorKind> + From<nom::error::ErrorKind>,
{
    alt((named(polymer, identifier), offset::<K>(polymer)))
}

// FIXME: I probably need to add a lot of `wrap_err`s around these parsers!
/// Named Modification = [ Multiplier ] , Identifier
pub fn named<'a, 'p, 's, K>(
    _polymer: &mut Polymer<'a, 'p>,
    mut identifier: impl Parser<&'s str, &'s str, LabeledParseError<'s, K>>,
) -> impl FnMut(&'s str) -> ParseResult<ModificationId, K>
where
    K: From<PolychemErrorKind> + From<nom::error::ErrorKind>,
{
    // FIXME: Obviously get rid of this lint-silencing todo hack
    let _ = identifier.parse("");
    |_| todo!()
}

/// Offset Modification = Offset Kind , [ Multiplier ] ,
///   Chemical Composition ;
pub fn offset<'a, 's, K>(
    _polymer: &mut Polymer<'a, '_>,
) -> impl FnMut(&'s str) -> ParseResult<ModificationId, K>
where
    K: From<PolychemErrorKind>,
{
    |_| todo!()
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
        sequence::pair,
    };
    use once_cell::sync::Lazy;

    use super::*;
    use crate::{AtomicDatabase, PolymerDatabase};

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
        macro_rules! assert_multiplier {
            ($input:literal, $output:literal, $count:expr) => {
                assert_eq!(
                    multiplier($input),
                    Ok(($output, Count::new($count).unwrap()))
                );
            };
        }
        // Valid Multipliers
        assert_multiplier!("1x", "", 1);
        assert_multiplier!("10x", "", 10);
        assert_multiplier!("422x", "", 422);
        assert_multiplier!("9999x", "", 9999);
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
        assert_multiplier!("1xOH", "OH", 1);
        assert_multiplier!("42xHeH", "HeH", 42);
    }

    #[ignore]
    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_named_modification() {
        // TODO: Restore from git
        todo!()
    }

    #[ignore]
    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_offset_modification() {
        // TODO: Restore from git
        todo!()
    }
}
