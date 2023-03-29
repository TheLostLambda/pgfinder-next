use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::{alpha1, alphanumeric1, char, one_of, satisfy, space0, u8},
    combinator::{map, recognize},
    multi::{many0, separated_list1},
    sequence::{delimited, pair, terminated, tuple},
    IResult, Parser,
};
use serde::{Deserialize, Serialize};

type Multimer = ();
type Monomer = ();
type Glycan = ();
type Monosaccharide = char;
type Modifications = Vec<Modification>;
type Moiety = String;
type Peptide = ();
type AminoAcid = char;
type LateralChain = ();

// FIXME: Add missing derives!
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum Modification {
    Add(Moiety),
    Remove(Moiety),
}

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum Crosslink {
    Ambiguous,
    Between(u8, u8),
}

/// Multimer = Monomer , { ( "~" (* Glycosidic Bond *) | Crosslink ) , Monomer } ;
pub fn multimer(i: &str) -> IResult<&str, Multimer> {
    todo!()
}

/// Monomer = Glycan , [ "-" , Peptide ] ;
pub fn monomer(i: &str) -> IResult<&str, Monomer> {
    todo!()
}

/// Glycan = Monosaccharide , [ Modifications ] ,
///   { Monosaccharide , [ Modifications ] } ;
pub fn glycan(i: &str) -> IResult<&str, Glycan> {
    todo!()
}

/// Monosaccharide = lowercase ;
pub fn monosaccharide(i: &str) -> IResult<&str, Monosaccharide> {
    satisfy(|c| c.is_ascii_lowercase())(i)
}

/// Modifications = "(" , ( "+" | "-" ) , Moiety ,
///   { "," , { space } , ( "+" | "-" ) , Moiety } , ")" ;
pub fn modifications(i: &str) -> IResult<&str, Modifications> {
    let parser = delimited(
        char('('),
        separated_list1(terminated(char(','), space0), tuple((one_of("+-"), moiety))),
        char(')'),
    );
    map(parser, |mods| {
        mods.into_iter()
            .map(|(sign, moiety)| match sign {
                '+' => Modification::Add(moiety),
                '-' => Modification::Remove(moiety),
                _ => unreachable!(),
            })
            .collect()
    })(i)
}

/// Moiety = letter , { letter | digit | "_" } ;
pub fn moiety(i: &str) -> IResult<&str, Moiety> {
    map(
        recognize(pair(alpha1, many0(alt((alphanumeric1, tag("_")))))),
        str::to_string,
    )(i)
}

/// Peptide = AminoAcid , [ Modifications ] ,
///   [ LateralChain ] , { AminoAcid , [ Modifications ] ,
///   [ LateralChain ] } ;
pub fn peptide(i: &str) -> IResult<&str, Peptide> {
    todo!()
}

/// AminoAcid = uppercase ;
pub fn amino_acid(i: &str) -> IResult<&str, AminoAcid> {
    satisfy(|c| c.is_ascii_uppercase())(i)
}

/// LateralChain = "[" , AminoAcid , [ Modifications ] ,
///   { AminoAcid , [ Modifications ] } , "]" ;
pub fn lateral_chain(i: &str) -> IResult<&str, LateralChain> {
    // delimited(char('['), peptide, char(']'))(i)
    todo!()
}

/// Crosslink
///   = "=" (* Ambiguous Crosslink *)
///   | "=" , position , ( "<" | ">" ) , position , "=" ,
///     { ":" , "=" , position , ( "<" | ">" ) , position , "=" }
///   ;
pub fn crosslink(i: &str) -> IResult<&str, Crosslink> {
    let ambiguous = char('=').map(|_| Crosslink::Ambiguous);
    let between = delimited(
        char('='),
        tuple((one_of("12345"), one_of("<>"), one_of("12345"))),
        char('='),
    )
    .map(|(a, dir, b)| {
        let a = a.to_digit(10).unwrap() as u8;
        let b = b.to_digit(10).unwrap() as u8;
        match dir {
            '<' => Crosslink::Between(b, a),
            '>' => Crosslink::Between(a, b),
            _ => unreachable!(),
        }
    });
    alt((ambiguous, between))(i)
}

#[cfg(test)]
mod tests {
    use std::error::Error;

    use super::*;

    #[test]
    fn test_monosaccharide() -> Result<(), Box<dyn Error>> {
        // Ensure the complete lowercase ASCII alphabet is present
        for c in 'a'..='z' {
            assert_eq!(monosaccharide(&c.to_string()).unwrap(), ("", c));
        }
        Ok(())
    }

    #[test]
    fn test_modifications() -> Result<(), Box<dyn Error>> {
        assert_eq!(
            modifications("(+Ac)")?,
            ("", vec![Modification::Add("Ac".to_string())])
        );
        assert_eq!(
            modifications("(-Ac)")?,
            ("", vec![Modification::Remove("Ac".to_string())])
        );
        assert_eq!(
            modifications("(+Ac,+Am)")?,
            (
                "",
                vec![
                    Modification::Add("Ac".to_string()),
                    Modification::Add("Am".to_string())
                ]
            )
        );
        assert_eq!(
            modifications("(+Ac, +Am)")?,
            (
                "",
                vec![
                    Modification::Add("Ac".to_string()),
                    Modification::Add("Am".to_string())
                ]
            )
        );
        assert_eq!(
            modifications("(+Ac, -Am)")?,
            (
                "",
                vec![
                    Modification::Add("Ac".to_string()),
                    Modification::Remove("Am".to_string())
                ]
            )
        );
        Ok(())
    }

    // TODO: Finish filling in these tests!
    #[test]
    fn test_moiety() -> Result<(), Box<dyn Error>> {
        Ok(())
    }

    #[test]
    fn test_amino_acid() -> Result<(), Box<dyn Error>> {
        // Ensure the complete uppercase ASCII alphabet is present
        for c in 'A'..='Z' {
            assert_eq!(amino_acid(&c.to_string()).unwrap(), ("", c));
        }
        Ok(())
    }

    #[test]
    fn test_crosslink() -> Result<(), Box<dyn Error>> {
        Ok(())
    }
}
