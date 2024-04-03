use nom::{
    character::complete::{char, one_of, satisfy, u32},
    combinator::{cut, map, not},
    sequence::preceded,
};
use nom_miette::expect;

use crate::{Count, OffsetKind};

use self::errors::{CompositionErrorKind, ParseResult};

mod chemical_composition;
pub mod errors;

// Re-exports:
pub use chemical_composition::chemical_composition;

/// Count = digit - "0" , { digit } ;
pub fn count(i: &str) -> ParseResult<Count> {
    let not_zero = expect(
        cut(not(char('0'))),
        CompositionErrorKind::ExpectedNoLeadingZero,
    );
    let digits = expect(u32, CompositionErrorKind::ExpectedDigit);
    preceded(not_zero, digits)(i)
}

/// Offset Kind = "+" | "-" ;
pub fn offset_kind(i: &str) -> ParseResult<OffsetKind> {
    map(one_of("+-"), |c| match c {
        '+' => OffsetKind::Add,
        '-' => OffsetKind::Remove,
        _ => unreachable!(),
    })(i)
}

/// uppercase
///   = "A" | "B" | "C" | "D" | "E" | "F" | "G"
///   | "H" | "I" | "J" | "K" | "L" | "M" | "N"
///   | "O" | "P" | "Q" | "R" | "S" | "T" | "U"
///   | "V" | "W" | "X" | "Y" | "Z"
///   ;
pub fn uppercase(i: &str) -> ParseResult<char> {
    let parser = satisfy(|c| c.is_ascii_uppercase());
    expect(parser, CompositionErrorKind::ExpectedUppercase)(i)
}

/// lowercase
///   = "a" | "b" | "c" | "d" | "e" | "f" | "g"
///   | "h" | "i" | "j" | "k" | "l" | "m" | "n"
///   | "o" | "p" | "q" | "r" | "s" | "t" | "u"
///   | "v" | "w" | "x" | "y" | "z"
///   ;
pub fn lowercase(i: &str) -> ParseResult<char> {
    let parser = satisfy(|c| c.is_ascii_lowercase());
    expect(parser, CompositionErrorKind::ExpectedLowercase)(i)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_uppercase() {
        // Ensure the complete uppercase ASCII alphabet is present
        for c in 'A'..='Z' {
            assert_eq!(uppercase(&c.to_string()), Ok(("", c)));
        }
        // Ensure the complete lowercase ASCII alphabet is absent
        for c in 'a'..='z' {
            assert!(uppercase(&c.to_string()).is_err());
        }
        // Ensure only one character is parsed
        assert_eq!(uppercase("Hg"), Ok(("g", 'H')));
        assert_eq!(uppercase("HG"), Ok(("G", 'H')));
    }

    #[test]
    fn test_lowercase() {
        // Ensure the complete lowercase ASCII alphabet is present
        for c in 'a'..='z' {
            assert_eq!(lowercase(&c.to_string()), Ok(("", c)));
        }
        // Ensure the complete uppercase ASCII alphabet is absent
        for c in 'A'..='Z' {
            assert!(lowercase(&c.to_string()).is_err());
        }
        // Ensure only one character is parsed
        assert_eq!(lowercase("hg"), Ok(("g", 'h')));
        assert_eq!(lowercase("hG"), Ok(("G", 'h')));
    }

    #[test]
    fn test_count() {
        // Valid Counts
        assert_eq!(count("1"), Ok(("", 1)));
        assert_eq!(count("10"), Ok(("", 10)));
        assert_eq!(count("422"), Ok(("", 422)));
        assert_eq!(count("9999"), Ok(("", 9999)));
        // Invalid Counts
        assert!(count("0").is_err());
        assert!(count("01").is_err());
        assert!(count("00145").is_err());
        assert!(count("H").is_err());
        assert!(count("p").is_err());
        assert!(count("+H").is_err());
        assert!(count("[H]").is_err());
        // Multiple Counts
        assert_eq!(count("1OH"), Ok(("OH", 1)));
        assert_eq!(count("42HeH"), Ok(("HeH", 42)));
    }
}
