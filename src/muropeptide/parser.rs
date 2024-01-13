use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::{alpha1, alphanumeric1, char, one_of, satisfy},
    combinator::{map, opt, recognize},
    multi::{many0, many1, separated_list1},
    sequence::{delimited, pair, preceded, tuple},
    IResult, Parser,
};
use serde::{Deserialize, Serialize};

// TODO: Write a "reverse-parser" Display impl that goes from AST to string
// TODO: Should all of these tuples be structs with named fields?
// FIXME: Make as much private as possible!
// FIXME: Check that all EBNF is up to date!
pub type Multimer = (Vec<Monomer>, Vec<Crosslink>);
pub type Monomer = (Glycan, Option<Peptide>);
pub type Glycan = Vec<(Monosaccharide, Option<Modifications>)>;
pub type Monosaccharide = char;
pub type Modifications = Vec<Modification>;
pub type Moiety = String;
pub type Peptide = Vec<(AminoAcid, Option<Modifications>, Option<LateralChain>)>;
pub type AminoAcid = char;
pub type LateralChain = Vec<(AminoAcid, Option<Modifications>)>;

// FIXME: Add missing derives!
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum Modification {
    Add(Moiety),
    Remove(Moiety),
}

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum Crosslink {
    Glycosidic,
    Ambiguous,
    // TODO: Would this be better as `Between(Vec<Bond>)`, where `Bond` contains the
    // `DonatesTo(u8, u8)` and `AcceptsFrom(u8, u8)` variants?
    Explicit(Vec<(u8, BondDirection, u8)>),
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum BondDirection {
    DonatesTo,
    AcceptsFrom,
}

/// Multimer = Monomer , { Crosslink , Monomer } ;
pub fn multimer(i: &str) -> IResult<&str, Multimer> {
    // NOTE: I don't love parsing the string twice, but I can't justify a more complex single-pass
    // solution without running some performance benchmarks...
    let (o, monomers) = separated_list1(crosslink, monomer)(i)?;
    let (_, crosslinks) = opt(delimited(
        monomer,
        separated_list1(monomer, crosslink),
        monomer,
    ))(i)?;
    Ok((o, (monomers, crosslinks.unwrap_or_default())))
}

/// Monomer = Glycan , [ "-" , Peptide ] ;
pub fn monomer(i: &str) -> IResult<&str, Monomer> {
    pair(glycan, opt(preceded(char('-'), peptide)))(i)
}

/// Glycan = Monosaccharide , [ Modifications ] ,
///   { Monosaccharide , [ Modifications ] } ;
pub fn glycan(i: &str) -> IResult<&str, Glycan> {
    many1(pair(monosaccharide, opt(modifications)))(i)
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
        separated_list1(char(','), tuple((one_of("+-"), moiety))),
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
    many1(tuple((amino_acid, opt(modifications), opt(lateral_chain))))(i)
}

/// AminoAcid = uppercase ;
pub fn amino_acid(i: &str) -> IResult<&str, AminoAcid> {
    satisfy(|c| c.is_ascii_uppercase())(i)
}

/// LateralChain = "[" , AminoAcid , [ Modifications ] ,
///   { AminoAcid , [ Modifications ] } , "]" ;
pub fn lateral_chain(i: &str) -> IResult<&str, LateralChain> {
    delimited(
        char('['),
        many1(pair(amino_acid, opt(modifications))),
        char(']'),
    )(i)
}

/// Crosslink
///   = "=" , position , ( "<" | ">" ) , position , "=" ,
///     { ":" , "=" , position , ( "<" | ">" ) , position , "=" }
///   | "=" (* Ambiguous Crosslink *)
///   | "~" (* Glycosidic Bond *)
///   ;
pub fn crosslink(i: &str) -> IResult<&str, Crosslink> {
    let glycosidic = char('~').map(|_| Crosslink::Glycosidic);
    let ambiguous = char('=').map(|_| Crosslink::Ambiguous);
    let between = delimited(
        char('='),
        separated_list1(
            char(':'),
            tuple((one_of("12345"), one_of("<>"), one_of("12345"))),
        ),
        char('='),
    )
    .map(|v| {
        Crosslink::Explicit(
            v.into_iter()
                .map(|(a, dir, b)| {
                    let a = a.to_digit(10).unwrap() as u8;
                    let b = b.to_digit(10).unwrap() as u8;
                    (
                        a,
                        match dir {
                            '<' => BondDirection::AcceptsFrom,
                            '>' => BondDirection::DonatesTo,
                            _ => unreachable!(),
                        },
                        b,
                    )
                })
                .collect(),
        )
    });
    alt((between, ambiguous, glycosidic))(i)
}

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use std::error::Error;

    use super::*;

    #[ignore]
    #[test]
    fn test_multimer() -> Result<(), Box<dyn Error>> {
        assert_debug_snapshot!(multimer("gm~gm-AEJA")?);
        assert_debug_snapshot!(multimer("gm~gm(-Ac)-AEJA")?);
        assert_debug_snapshot!(multimer("gm-AEJA=gm-AEJ")?);
        assert_debug_snapshot!(multimer("gm(+An,-Ac)~gm")?);
        assert_debug_snapshot!(multimer("gm-AEJA=4>3=gm-AEJ")?);
        assert_debug_snapshot!(multimer("gm-AEJA=3<3=gm-AEJ")?);
        assert_debug_snapshot!(multimer("gm-AEJ=3<3:3>3=gm-AQJ")?);
        assert_debug_snapshot!(multimer("gm-AEJA=4>3=gm-AEJ=3>3=gm-AQJA")?);
        assert_debug_snapshot!(multimer("gm-AEJA=3<3=gm-AEJ=3<4=gm-AQJA")?);
        Ok(())
    }

    #[ignore]
    #[test]
    fn test_monomer() -> Result<(), Box<dyn Error>> {
        assert_debug_snapshot!(monomer("gm")?);
        assert_debug_snapshot!(monomer("gm-AEJA")?);
        assert_debug_snapshot!(monomer("m-AEJA")?);
        assert_debug_snapshot!(monomer("gm-AE(+Am)JA")?);
        assert_debug_snapshot!(monomer("gm(+An,-Ac)-AEJ")?);
        assert_debug_snapshot!(monomer("gm-AQK[AA]AA")?);
        Ok(())
    }

    #[ignore]
    #[test]
    fn test_glycan() -> Result<(), Box<dyn Error>> {
        assert_debug_snapshot!(glycan("g")?);
        assert_debug_snapshot!(glycan("m")?);
        assert_debug_snapshot!(glycan("gm")?);
        assert_debug_snapshot!(glycan("gm(-Ac)")?);
        assert_debug_snapshot!(glycan("g(+Am)m(-Ac)")?);
        assert_debug_snapshot!(glycan("gm(-Ac,+Am)")?);
        Ok(())
    }

    #[ignore]
    #[test]
    fn test_monosaccharide() {
        // Ensure the complete lowercase ASCII alphabet is present
        for c in 'a'..='z' {
            assert_eq!(monosaccharide(&c.to_string()).unwrap(), ("", c));
        }
    }

    #[ignore]
    #[test]
    fn test_modifications() -> Result<(), Box<dyn Error>> {
        assert!(modifications("()").is_err());
        assert_debug_snapshot!(modifications("(+Ac)")?);
        assert_debug_snapshot!(modifications("(-Ac)")?);
        assert_debug_snapshot!(modifications("(+Ac,+Am)")?);
        assert_debug_snapshot!(modifications("(+Ac,-Am,+OH)")?);
        Ok(())
    }

    #[ignore]
    #[test]
    fn test_moiety() -> Result<(), Box<dyn Error>> {
        assert_eq!(moiety("Ac")?, ("", "Ac".to_string()));
        assert_eq!(moiety("h2o")?, ("", "h2o".to_string()));
        assert_eq!(moiety("h_2o")?, ("", "h_2o".to_string()));
        assert!(moiety("+Ac").is_err());
        assert!(moiety("2Am").is_err());
        assert!(moiety("_Am").is_err());
        assert_eq!(moiety("Ac+H2O")?, ("+H2O", "Ac".to_string()));
        Ok(())
    }

    #[ignore]
    #[test]
    fn test_peptide() -> Result<(), Box<dyn Error>> {
        assert_debug_snapshot!(peptide("A")?);
        assert_debug_snapshot!(peptide("E")?);
        assert_debug_snapshot!(peptide("AEJA")?);
        assert_debug_snapshot!(peptide("AEJ(+Am)")?);
        assert_debug_snapshot!(peptide("A(+Am)E(-Ac)")?);
        assert_debug_snapshot!(peptide("AE(-Ac,+Am)")?);
        assert_debug_snapshot!(peptide("AQK[AA]AA")?);
        assert_debug_snapshot!(peptide("AQK(-H)[AA(+Am)]AA")?);
        Ok(())
    }

    #[ignore]
    #[test]
    fn test_amino_acid() {
        // Ensure the complete uppercase ASCII alphabet is present
        for c in 'A'..='Z' {
            assert_eq!(amino_acid(&c.to_string()).unwrap(), ("", c));
        }
    }

    #[ignore]
    #[test]
    fn test_lateral_chain() -> Result<(), Box<dyn Error>> {
        assert_debug_snapshot!(lateral_chain("[A]")?);
        assert_debug_snapshot!(lateral_chain("[G]")?);
        assert_debug_snapshot!(lateral_chain("[AG]")?);
        assert_debug_snapshot!(lateral_chain("[AG(-Ac)]")?);
        assert_debug_snapshot!(lateral_chain("[A(+Am)G(-Ac)]")?);
        assert_debug_snapshot!(lateral_chain("[AG(-Ac,+Am)]")?);
        Ok(())
    }

    #[ignore]
    #[test]
    fn test_crosslink() -> Result<(), Box<dyn Error>> {
        assert_debug_snapshot!(crosslink("~")?);
        assert_debug_snapshot!(crosslink("=")?);
        assert_debug_snapshot!(crosslink("=4>3=")?);
        assert_debug_snapshot!(crosslink("=3<4=")?);
        assert_debug_snapshot!(crosslink("=3<3:3>3=")?);
        Ok(())
    }
}
