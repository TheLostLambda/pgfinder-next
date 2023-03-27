use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::{alpha1, alphanumeric1, char, one_of, satisfy, space0},
    combinator::{map, recognize},
    multi::{many0, separated_list1},
    sequence::{delimited, pair, terminated, tuple},
    IResult,
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
type Crosslink = ();

// FIXME: Add missing derives!
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum Modification {
    Add(Moiety),
    Remove(Moiety),
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

/// LateralChain = "[" , Peptide , "]" ;
pub fn lateral_chain(i: &str) -> IResult<&str, LateralChain> {
    todo!()
}

/// Crosslink
///   = "=" (* Ambiguous Crosslink *)
///   | "=" , digit , ( "<" | ">" ) , digit , "=" ,
///     { ":" , "=" , digit , ( "<" | ">" ) , digit , "=" }
///   ;
pub fn crosslink(i: &str) -> IResult<&str, Crosslink> {
    todo!()
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
    fn test_amino_acid() -> Result<(), Box<dyn Error>> {
        // Ensure the complete uppercase ASCII alphabet is present
        for c in 'A'..='Z' {
            assert_eq!(amino_acid(&c.to_string()).unwrap(), ("", c));
        }
        Ok(())
    }

    // TODO: Write tests for `moiety` and `modifications`!
}
