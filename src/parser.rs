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
    use std::error::Error;

    use super::*;

    #[test]
    fn test_multimer() -> Result<(), Box<dyn Error>> {
        assert_eq!(
            multimer("gm~gm-AEJA")?,
            (
                "",
                (
                    vec![
                        (vec![('g', None), ('m', None)], None),
                        (
                            vec![('g', None), ('m', None)],
                            Some(vec![
                                ('A', None, None),
                                ('E', None, None),
                                ('J', None, None),
                                ('A', None, None)
                            ])
                        )
                    ],
                    vec![Crosslink::Glycosidic]
                )
            )
        );
        assert_eq!(
            multimer("gm~gm(-Ac)-AEJA")?,
            (
                "",
                (
                    vec![
                        (vec![('g', None), ('m', None)], None),
                        (
                            vec![
                                ('g', None),
                                ('m', Some(vec![Modification::Remove("Ac".to_string())]))
                            ],
                            Some(vec![
                                ('A', None, None),
                                ('E', None, None),
                                ('J', None, None),
                                ('A', None, None)
                            ])
                        )
                    ],
                    vec![Crosslink::Glycosidic]
                )
            )
        );
        assert_eq!(
            multimer("gm-AEJA=gm-AEJ")?,
            (
                "",
                (
                    vec![
                        (
                            vec![('g', None), ('m', None)],
                            Some(vec![
                                ('A', None, None),
                                ('E', None, None),
                                ('J', None, None),
                                ('A', None, None)
                            ])
                        ),
                        (
                            vec![('g', None), ('m', None)],
                            Some(vec![
                                ('A', None, None),
                                ('E', None, None),
                                ('J', None, None),
                            ])
                        )
                    ],
                    vec![Crosslink::Ambiguous]
                )
            )
        );
        assert_eq!(
            multimer("gm(+An,-Ac)~gm")?,
            (
                "",
                (
                    vec![
                        (
                            vec![
                                ('g', None),
                                (
                                    'm',
                                    Some(vec![
                                        Modification::Add("An".to_string()),
                                        Modification::Remove("Ac".to_string())
                                    ])
                                )
                            ],
                            None
                        ),
                        (vec![('g', None), ('m', None)], None)
                    ],
                    vec![Crosslink::Glycosidic]
                )
            )
        );
        assert_eq!(
            multimer("gm-AEJA=4>3=gm-AEJ")?,
            (
                "",
                (
                    vec![
                        (
                            vec![('g', None), ('m', None)],
                            Some(vec![
                                ('A', None, None),
                                ('E', None, None),
                                ('J', None, None),
                                ('A', None, None)
                            ])
                        ),
                        (
                            vec![('g', None), ('m', None)],
                            Some(vec![
                                ('A', None, None),
                                ('E', None, None),
                                ('J', None, None)
                            ])
                        )
                    ],
                    vec![Crosslink::Explicit(vec![(4, BondDirection::DonatesTo, 3)])]
                )
            )
        );
        assert_eq!(
            multimer("gm-AEJA=3<3=gm-AEJ")?,
            (
                "",
                (
                    vec![
                        (
                            vec![('g', None), ('m', None)],
                            Some(vec![
                                ('A', None, None),
                                ('E', None, None),
                                ('J', None, None),
                                ('A', None, None)
                            ])
                        ),
                        (
                            vec![('g', None), ('m', None)],
                            Some(vec![
                                ('A', None, None),
                                ('E', None, None),
                                ('J', None, None)
                            ])
                        )
                    ],
                    vec![Crosslink::Explicit(vec![(
                        3,
                        BondDirection::AcceptsFrom,
                        3
                    )])]
                )
            )
        );
        assert_eq!(
            multimer("gm-AEJ=3<3:3>3=gm-AQJ")?,
            (
                "",
                (
                    vec![
                        (
                            vec![('g', None), ('m', None)],
                            Some(vec![
                                ('A', None, None),
                                ('E', None, None),
                                ('J', None, None),
                            ])
                        ),
                        (
                            vec![('g', None), ('m', None)],
                            Some(vec![
                                ('A', None, None),
                                ('Q', None, None),
                                ('J', None, None)
                            ])
                        )
                    ],
                    vec![Crosslink::Explicit(vec![
                        (3, BondDirection::AcceptsFrom, 3),
                        (3, BondDirection::DonatesTo, 3)
                    ])]
                )
            )
        );
        assert_eq!(
            multimer("gm-AEJA=4>3=gm-AEJ=3>3=gm-AQJA")?,
            (
                "",
                (
                    vec![
                        (
                            vec![('g', None), ('m', None)],
                            Some(vec![
                                ('A', None, None),
                                ('E', None, None),
                                ('J', None, None),
                                ('A', None, None)
                            ])
                        ),
                        (
                            vec![('g', None), ('m', None)],
                            Some(vec![
                                ('A', None, None),
                                ('E', None, None),
                                ('J', None, None)
                            ])
                        ),
                        (
                            vec![('g', None), ('m', None)],
                            Some(vec![
                                ('A', None, None),
                                ('Q', None, None),
                                ('J', None, None),
                                ('A', None, None)
                            ])
                        )
                    ],
                    vec![
                        Crosslink::Explicit(vec![(4, BondDirection::DonatesTo, 3)]),
                        Crosslink::Explicit(vec![(3, BondDirection::DonatesTo, 3)])
                    ]
                )
            )
        );
        assert_eq!(
            multimer("gm-AEJA=3<3=gm-AEJ=3<4=gm-AQJA")?,
            (
                "",
                (
                    vec![
                        (
                            vec![('g', None), ('m', None)],
                            Some(vec![
                                ('A', None, None),
                                ('E', None, None),
                                ('J', None, None),
                                ('A', None, None)
                            ])
                        ),
                        (
                            vec![('g', None), ('m', None)],
                            Some(vec![
                                ('A', None, None),
                                ('E', None, None),
                                ('J', None, None)
                            ])
                        ),
                        (
                            vec![('g', None), ('m', None)],
                            Some(vec![
                                ('A', None, None),
                                ('Q', None, None),
                                ('J', None, None),
                                ('A', None, None)
                            ])
                        )
                    ],
                    vec![
                        Crosslink::Explicit(vec![(3, BondDirection::AcceptsFrom, 3)]),
                        Crosslink::Explicit(vec![(3, BondDirection::AcceptsFrom, 4)])
                    ]
                )
            )
        );
        Ok(())
    }

    #[test]
    fn test_monomer() -> Result<(), Box<dyn Error>> {
        assert_eq!(monomer("gm")?, ("", (vec![('g', None), ('m', None)], None)));
        assert_eq!(
            monomer("gm-AEJA")?,
            (
                "",
                (
                    vec![('g', None), ('m', None)],
                    Some(vec![
                        ('A', None, None),
                        ('E', None, None),
                        ('J', None, None),
                        ('A', None, None)
                    ])
                )
            )
        );
        assert_eq!(
            monomer("m-AEJA")?,
            (
                "",
                (
                    vec![('m', None)],
                    Some(vec![
                        ('A', None, None),
                        ('E', None, None),
                        ('J', None, None),
                        ('A', None, None)
                    ])
                )
            )
        );
        assert_eq!(
            monomer("gm-AE(+Am)JA")?,
            (
                "",
                (
                    vec![('g', None), ('m', None)],
                    Some(vec![
                        ('A', None, None),
                        ('E', Some(vec![Modification::Add("Am".to_string())]), None),
                        ('J', None, None),
                        ('A', None, None)
                    ])
                )
            )
        );
        assert_eq!(
            monomer("gm(+An,-Ac)-AEJ")?,
            (
                "",
                (
                    vec![
                        ('g', None),
                        (
                            'm',
                            Some(vec![
                                Modification::Add("An".to_string()),
                                Modification::Remove("Ac".to_string())
                            ])
                        )
                    ],
                    Some(vec![
                        ('A', None, None),
                        ('E', None, None),
                        ('J', None, None),
                    ])
                )
            )
        );
        assert_eq!(
            monomer("gm-AQK[AA]AA")?,
            (
                "",
                (
                    vec![('g', None), ('m', None)],
                    Some(vec![
                        ('A', None, None),
                        ('Q', None, None),
                        ('K', None, Some(vec![('A', None), ('A', None)])),
                        ('A', None, None),
                        ('A', None, None),
                    ])
                )
            )
        );
        Ok(())
    }

    #[test]
    fn test_glycan() -> Result<(), Box<dyn Error>> {
        assert_eq!(glycan("g")?, ("", vec![('g', None)]));
        assert_eq!(glycan("m")?, ("", vec![('m', None)]));
        assert_eq!(glycan("gm")?, ("", vec![('g', None), ('m', None)]));
        assert_eq!(
            glycan("gm(-Ac)")?,
            (
                "",
                vec![
                    ('g', None),
                    ('m', Some(vec![Modification::Remove("Ac".to_string())]))
                ]
            )
        );
        assert_eq!(
            glycan("g(+Am)m(-Ac)")?,
            (
                "",
                vec![
                    ('g', Some(vec![Modification::Add("Am".to_string())])),
                    ('m', Some(vec![Modification::Remove("Ac".to_string())]))
                ]
            )
        );
        assert_eq!(
            glycan("gm(-Ac,+Am)")?,
            (
                "",
                vec![
                    ('g', None),
                    (
                        'm',
                        Some(vec![
                            Modification::Remove("Ac".to_string()),
                            Modification::Add("Am".to_string())
                        ])
                    )
                ]
            )
        );
        Ok(())
    }

    #[test]
    fn test_monosaccharide() {
        // Ensure the complete lowercase ASCII alphabet is present
        for c in 'a'..='z' {
            assert_eq!(monosaccharide(&c.to_string()).unwrap(), ("", c));
        }
    }

    #[test]
    fn test_modifications() -> Result<(), Box<dyn Error>> {
        assert!(modifications("()").is_err());
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
            modifications("(+Ac,-Am,+OH)")?,
            (
                "",
                vec![
                    Modification::Add("Ac".to_string()),
                    Modification::Remove("Am".to_string()),
                    Modification::Add("OH".to_string())
                ]
            )
        );
        Ok(())
    }

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

    #[test]
    fn test_peptide() -> Result<(), Box<dyn Error>> {
        assert_eq!(peptide("A")?, ("", vec![('A', None, None)]));
        assert_eq!(peptide("E")?, ("", vec![('E', None, None)]));
        assert_eq!(
            peptide("AEJA")?,
            (
                "",
                vec![
                    ('A', None, None),
                    ('E', None, None),
                    ('J', None, None),
                    ('A', None, None)
                ]
            )
        );
        assert_eq!(
            peptide("AEJ(+Am)")?,
            (
                "",
                vec![
                    ('A', None, None),
                    ('E', None, None),
                    ('J', Some(vec![Modification::Add("Am".to_string())]), None)
                ]
            )
        );
        assert_eq!(
            peptide("A(+Am)E(-Ac)")?,
            (
                "",
                vec![
                    ('A', Some(vec![Modification::Add("Am".to_string())]), None),
                    (
                        'E',
                        Some(vec![Modification::Remove("Ac".to_string())]),
                        None
                    )
                ]
            )
        );
        assert_eq!(
            peptide("AE(-Ac,+Am)")?,
            (
                "",
                vec![
                    ('A', None, None),
                    (
                        'E',
                        Some(vec![
                            Modification::Remove("Ac".to_string()),
                            Modification::Add("Am".to_string())
                        ]),
                        None
                    )
                ]
            )
        );
        assert_eq!(
            peptide("AQK[AA]AA")?,
            (
                "",
                vec![
                    ('A', None, None),
                    ('Q', None, None),
                    ('K', None, Some(vec![('A', None), ('A', None)])),
                    ('A', None, None),
                    ('A', None, None),
                ]
            )
        );
        assert_eq!(
            peptide("AQK(-H)[AA(+Am)]AA")?,
            (
                "",
                vec![
                    ('A', None, None),
                    ('Q', None, None),
                    (
                        'K',
                        Some(vec![Modification::Remove("H".to_string())]),
                        Some(vec![
                            ('A', None),
                            ('A', Some(vec![Modification::Add("Am".to_string())]))
                        ])
                    ),
                    ('A', None, None),
                    ('A', None, None),
                ]
            )
        );
        Ok(())
    }

    #[test]
    fn test_amino_acid() {
        // Ensure the complete uppercase ASCII alphabet is present
        for c in 'A'..='Z' {
            assert_eq!(amino_acid(&c.to_string()).unwrap(), ("", c));
        }
    }

    #[test]
    fn test_lateral_chain() -> Result<(), Box<dyn Error>> {
        assert_eq!(lateral_chain("[A]")?, ("", vec![('A', None)]));
        assert_eq!(lateral_chain("[G]")?, ("", vec![('G', None)]));
        assert_eq!(lateral_chain("[AG]")?, ("", vec![('A', None), ('G', None)]));
        assert_eq!(
            lateral_chain("[AG(-Ac)]")?,
            (
                "",
                vec![
                    ('A', None),
                    ('G', Some(vec![Modification::Remove("Ac".to_string())]))
                ]
            )
        );
        assert_eq!(
            lateral_chain("[A(+Am)G(-Ac)]")?,
            (
                "",
                vec![
                    ('A', Some(vec![Modification::Add("Am".to_string())])),
                    ('G', Some(vec![Modification::Remove("Ac".to_string())]))
                ]
            )
        );
        assert_eq!(
            lateral_chain("[AG(-Ac,+Am)]")?,
            (
                "",
                vec![
                    ('A', None),
                    (
                        'G',
                        Some(vec![
                            Modification::Remove("Ac".to_string()),
                            Modification::Add("Am".to_string())
                        ])
                    )
                ]
            )
        );
        Ok(())
    }

    #[test]
    fn test_crosslink() -> Result<(), Box<dyn Error>> {
        assert_eq!(crosslink("~")?, ("", Crosslink::Glycosidic));
        assert_eq!(crosslink("=")?, ("", Crosslink::Ambiguous));
        assert_eq!(
            crosslink("=4>3=")?,
            (
                "",
                Crosslink::Explicit(vec![(4, BondDirection::DonatesTo, 3)])
            )
        );
        assert_eq!(
            crosslink("=3<4=")?,
            (
                "",
                Crosslink::Explicit(vec![(3, BondDirection::AcceptsFrom, 4)])
            )
        );
        assert_eq!(
            crosslink("=3<3:3>3=")?,
            (
                "",
                Crosslink::Explicit(vec![
                    (3, BondDirection::AcceptsFrom, 3),
                    (3, BondDirection::DonatesTo, 3)
                ])
            )
        );
        Ok(())
    }
}
