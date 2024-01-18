// FIXME: Starting with getting the "happy path" working, then adding in nice error reporting! Don't forget to test!
// FIXME: Make sure to update all of the `is_err()` tests to check that the error contains the correct, rich info...
// FIXME: Order functions in the same way as the EBNF

use nom::{
    branch::alt,
    character::complete::{char, satisfy, u32},
    combinator::{map, map_res, opt, recognize},
    sequence::{delimited, pair},
    IResult,
};

use super::{chemical_database::ChemicalDatabase, Element, MassNumber, OffsetKind, Particle};

/// Element = uppercase , [ lowercase ] ;
fn element_symbol(i: &str) -> IResult<&str, &str> {
    recognize(pair(uppercase, opt(lowercase)))(i)
}

/// Isotope = "[" , Integer , Element , "]" ;
fn isotope_expr(i: &str) -> IResult<&str, (MassNumber, &str)> {
    delimited(char('['), pair(u32, element_symbol), char(']'))(i)
}

/// Element = uppercase , [ lowercase ] ;
fn element<'a>(db: &'a ChemicalDatabase) -> impl FnMut(&'a str) -> IResult<&'a str, Element> {
    map_res(element_symbol, |symbol| Element::new(db, symbol))
}

/// Isotope = "[" , Integer , Element , "]" ;
fn isotope<'a>(db: &'a ChemicalDatabase) -> impl FnMut(&'a str) -> IResult<&'a str, Element> {
    map_res(isotope_expr, |(mass_number, symbol)| {
        Element::new_isotope(db, symbol, mass_number)
    })
}

/// Particle = lowercase ;
fn particle<'a>(db: &'a ChemicalDatabase) -> impl FnMut(&'a str) -> IResult<&'a str, Particle> {
    map_res(lowercase, |symbol| Particle::new(db, symbol.to_string()))
}

/// Offset Kind = "+" | "-" ;
fn offset_kind(i: &str) -> IResult<&str, OffsetKind> {
    map(alt((char('+'), char('-'))), |c| match c {
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
fn uppercase(i: &str) -> IResult<&str, char> {
    satisfy(|c| c.is_ascii_uppercase())(i)
}

/// lowercase
///   = "a" | "b" | "c" | "d" | "e" | "f" | "g"
///   | "h" | "i" | "j" | "k" | "l" | "m" | "n"
///   | "o" | "p" | "q" | "r" | "s" | "t" | "u"
///   | "v" | "w" | "x" | "y" | "z"
///   ;
fn lowercase(i: &str) -> IResult<&str, char> {
    satisfy(|c| c.is_ascii_lowercase())(i)
}

#[cfg(test)]
mod tests {
    use miette::Result;
    use once_cell::sync::Lazy;

    use super::*;

    static DB: Lazy<ChemicalDatabase> = Lazy::new(|| {
        ChemicalDatabase::from_kdl("chemistry.kdl", include_str!("chemistry.kdl")).unwrap()
    });

    #[test]
    fn test_uppercase() -> Result<()> {
        // Ensure the complete uppercase ASCII alphabet is present
        for c in 'A'..='Z' {
            assert_eq!(uppercase(&c.to_string()), Ok(("", c)))
        }
        // Ensure the complete lowercase ASCII alphabet is absent
        for c in 'a'..='z' {
            assert!(uppercase(&c.to_string()).is_err())
        }
        // Ensure only one character is parsed
        assert_eq!(uppercase("Hg"), Ok(("g", 'H')));
        assert_eq!(uppercase("HG"), Ok(("G", 'H')));
        // FIXME: Add a check for a richer error message when parsing fails
        Ok(())
    }

    #[test]
    fn test_lowercase() -> Result<()> {
        // Ensure the complete lowercase ASCII alphabet is present
        for c in 'a'..='z' {
            assert_eq!(lowercase(&c.to_string()), Ok(("", c)))
        }
        // Ensure the complete uppercase ASCII alphabet is absent
        for c in 'A'..='Z' {
            assert!(lowercase(&c.to_string()).is_err())
        }
        // Ensure only one character is parsed
        assert_eq!(lowercase("hg"), Ok(("g", 'h')));
        assert_eq!(lowercase("hG"), Ok(("G", 'h')));
        // FIXME: Add a check for a richer error message when parsing fails
        Ok(())
    }

    #[test]
    fn test_element_symbol() -> Result<()> {
        // Valid Element Symbols
        assert_eq!(element_symbol("H"), Ok(("", "H")));
        assert_eq!(element_symbol("He"), Ok(("", "He")));
        // Invalid Element Symbols
        assert!(element_symbol("p").is_err());
        assert!(element_symbol("ep").is_err());
        assert!(element_symbol("1H").is_err());
        assert!(element_symbol("+H").is_err());
        assert!(element_symbol("[H]").is_err());
        // Multiple Element Symbols
        assert_eq!(element_symbol("OH"), Ok(("H", "O")));
        assert_eq!(element_symbol("HeH"), Ok(("H", "He")));
        Ok(())
    }

    #[test]
    fn test_element() -> Result<()> {
        let mut element = element(&DB);
        macro_rules! assert_element_name {
            ($input:literal, $output:literal, $name:literal) => {
                assert_eq!(
                    element($input).map(|(r, e)| (r, e.name)),
                    Ok(($output, $name.to_string()))
                );
            };
        }
        // Valid Elements
        assert_element_name!("H", "", "Hydrogen");
        assert_element_name!("He", "", "Helium");
        // Invalid Elements
        assert!(element("p").is_err());
        assert!(element("ep").is_err());
        assert!(element("1H").is_err());
        assert!(element("+H").is_err());
        assert!(element("[H]").is_err());
        // Non-Existent Elements
        assert!(element("X").is_err());
        assert!(element("To").is_err());
        // Multiple Elements
        assert_element_name!("OH", "H", "Oxygen");
        assert_element_name!("HeH", "H", "Helium");
        Ok(())
    }

    #[test]
    fn test_isotope_expr() -> Result<()> {
        // Valid Isotope Expressions
        assert_eq!(isotope_expr("[1H]"), Ok(("", (1, "H"))));
        assert_eq!(isotope_expr("[18O]"), Ok(("", (18, "O"))));
        assert_eq!(isotope_expr("[37Cl]"), Ok(("", (37, "Cl"))));
        // Invalid Isotope Expressions
        assert!(isotope_expr("H").is_err());
        assert!(isotope_expr("[H]").is_err());
        assert!(isotope_expr("[H2]").is_err());
        assert!(isotope_expr("[18OH]").is_err());
        assert!(isotope_expr("[18]").is_err());
        assert!(isotope_expr("[[18O]]").is_err());
        assert!(isotope_expr("[-18O]").is_err());
        assert!(isotope_expr("[+18O]").is_err());
        // Multiple Isotope Expressions
        assert_eq!(isotope_expr("[13C]O2"), Ok(("O2", (13, "C"))));
        assert_eq!(isotope_expr("[3He]H"), Ok(("H", (3, "He"))));
        Ok(())
    }

    #[test]
    fn test_isotope() -> Result<()> {
        let mut isotope = isotope(&DB);
        macro_rules! assert_isotope_name {
            ($input:literal, $output:literal, $name:literal) => {
                assert_eq!(
                    isotope($input)
                        .map(|(r, e)| (r, format!("{}-{}", e.name, e.mass_number.unwrap()))),
                    Ok(($output, $name.to_string()))
                );
            };
        }
        // Valid Isotopes
        assert_isotope_name!("[1H]", "", "Hydrogen-1");
        assert_isotope_name!("[18O]", "", "Oxygen-18");
        assert_isotope_name!("[37Cl]", "", "Chlorine-37");
        // Invalid Isotopes
        assert!(isotope("H").is_err());
        assert!(isotope("[H]").is_err());
        assert!(isotope("[H2]").is_err());
        assert!(isotope("[18OH]").is_err());
        assert!(isotope("[18]").is_err());
        assert!(isotope("[[18O]]").is_err());
        assert!(isotope("[-18O]").is_err());
        assert!(isotope("[+18O]").is_err());
        // Non-Existent Isotopes
        assert!(isotope("[42X]").is_err());
        assert!(isotope("[99To]").is_err());
        assert!(isotope("[15C]").is_err());
        assert!(isotope("[100Tc]").is_err());
        // Multiple Isotopes
        assert_isotope_name!("[13C]O2", "O2", "Carbon-13");
        assert_isotope_name!("[3He]H", "H", "Helium-3");
        Ok(())
    }

    #[test]
    fn test_particle() -> Result<()> {
        let mut particle = particle(&DB);
        macro_rules! assert_particle_name {
            ($input:literal, $output:literal, $name:literal) => {
                assert_eq!(
                    particle($input).map(|(r, e)| (r, e.name)),
                    Ok(($output, $name.to_string()))
                );
            };
        }
        // Valid Particles
        assert_particle_name!("p", "", "Proton");
        assert_particle_name!("e", "", "Electron");
        // Invalid Particles
        assert!(particle("P").is_err());
        assert!(particle("Ep").is_err());
        assert!(particle("1p").is_err());
        assert!(particle("+e").is_err());
        assert!(particle("[p]").is_err());
        // Non-Existent Particles
        assert!(particle("m").is_err());
        assert!(particle("g").is_err());
        // Multiple Particles
        assert_particle_name!("ep", "p", "Electron");
        assert_particle_name!("pe", "e", "Proton");
        Ok(())
    }

    #[test]
    fn test_offset_kind() -> Result<()> {
        // Valid Offset Kinds
        assert_eq!(offset_kind("+"), Ok(("", OffsetKind::Add)));
        assert_eq!(offset_kind("-"), Ok(("", OffsetKind::Remove)));
        // Invalid Offset Kinds
        assert!(offset_kind("p").is_err());
        assert!(offset_kind("H").is_err());
        assert!(offset_kind("1H").is_err());
        assert!(offset_kind("1+H").is_err());
        assert!(offset_kind("[H]").is_err());
        // Multiple Offset Kinds
        assert_eq!(offset_kind("+-"), Ok(("-", OffsetKind::Add)));
        assert_eq!(offset_kind("--"), Ok(("-", OffsetKind::Remove)));
        Ok(())
    }
}
