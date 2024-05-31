use miette::Diagnostic;
use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::{alpha1, alphanumeric1, char, satisfy, u32},
    combinator::{cut, map, not, opt, recognize},
    error::ErrorKind,
    multi::{many0, many1},
    sequence::{pair, preceded, separated_pair, terminated},
    IResult,
};
use nom_miette::{expect, wrap_err, LabeledErrorKind, LabeledParseError};
use polychem::{parsers::errors::PolychemErrorKind, Count, Polymer};
use thiserror::Error;

use crate::{
    AminoAcid, LateralChain, Monomer, Monosaccharide, Multimer, ResidueAbbr, UnbranchedAminoAcid,
};

// FIXME: This needs a test, and is currently just forwarding the monomer parser!!!
/// Multimer = Monomer , { Connection , Monomer } , [ Connection ] , [ { " " }- ,
///   ( Modifications , [ { " " }- , Crosslinks ]
///   | Crosslinks , [ { " " }- , Modifications ]
///   ) ] ;
pub fn multimer(i: &str) -> ParseResult<Multimer<ResidueAbbr>> {
    // FIXME: Don't just use the `monomer` parser! There are other things to parse here!
    let mut parser = map(monomer, |monomer| Multimer {
        monomers: vec![monomer],
        connections: Vec::new(),
        modifications: Vec::new(),
    });
    // FIXME: Wrap this error?
    parser(i)
}

// FIXME: Paste all of these EBNF comments into another file and make sure they are valid!
/// Monomer = Glycan , "-" , Peptide | Peptide | Glycan ;
fn monomer(i: &str) -> ParseResult<Monomer<ResidueAbbr>> {
    let glycan_and_peptide = map(
        separated_pair(glycan, char('-'), peptide),
        |(glycan, peptide)| Monomer { glycan, peptide },
    );

    let peptide_only = map(peptide, |peptide| Monomer {
        glycan: Vec::new(),
        peptide,
    });

    let glycan_only = map(glycan, |glycan| Monomer {
        glycan,
        peptide: Vec::new(),
    });
    let mut parser = alt((glycan_and_peptide, peptide_only, glycan_only));
    // FIXME: Wrap this error?
    parser(i)
}

// =

/// Glycan = { Monosaccharide }- ;
fn glycan(i: &str) -> ParseResult<Vec<Monosaccharide<ResidueAbbr>>> {
    let mut parser = many1(monosaccharide);
    // FIXME: Wrap this error?
    parser(i)
}

// FIXME: This is using the wrong amino acid parser — needs lateral chain support!
/// Peptide = { Amino Acid }- ;
fn peptide(i: &str) -> ParseResult<Vec<UnbranchedAminoAcid<ResidueAbbr>>> {
    // FIXME: Change to branched amino acid!
    let mut parser = many1(unbranched_amino_acid);
    // FIXME: Wrap this error?
    parser(i)
}

// =

// FIXME: Add modifications
/// Monosaccharide = lowercase , [ Modifications ] ;
fn monosaccharide(i: &str) -> ParseResult<Monosaccharide<ResidueAbbr>> {
    let mut parser = recognize(lowercase);
    // FIXME: Wrap this error?
    parser(i)
}

/// Amino Acid = Unbranched Amino Acid , [ Lateral Chain ] ;
fn amino_acid<'a, 'p, 's>(
    _polymer: &mut Polymer<'a, 'p>,
) -> impl FnMut(&'s str) -> ParseResult<AminoAcid<ResidueAbbr>> {
    |_| todo!()
}

// =

// /// Modifications = "(" , Any Modification ,
// ///   { { " " } , "," , { " " } , Any Modification } , ")" ;
// fn modifications<'a, 'p, 's>(
//     polymer: &mut Polymer<'a, 'p>,
// ) -> impl FnMut(&'s str) -> ParseResult<Vec<ModificationId>> {
//     let separator = delimited(space0, char(','), space0);
//     delimited(
//         char('('),
//         separated_list1(separator, any_modification(polymer, identifier)),
//         char(')'),
//     )
// }

// FIXME: Add modifications
/// Unbranched Amino Acid = [ lowercase ] , uppercase , [ Modifications ] ;
fn unbranched_amino_acid(i: &str) -> ParseResult<UnbranchedAminoAcid<ResidueAbbr>> {
    let mut parser = recognize(preceded(opt(lowercase), uppercase));
    // FIXME: Wrap this error?
    parser(i)
}

// NOTE: These are not meant to be links, it's just EBNF
#[allow(clippy::doc_link_with_quotes)]
/// Lateral Chain = "[" , [ "<" (* C-to-N *) | ">" (* N-to-C *) ] ,
///   { Unbranched Amino Acid }- , "]" ;
fn lateral_chain<'a, 'p, 's>(
    _polymer: &mut Polymer<'a, 'p>,
) -> impl FnMut(&'s str) -> ParseResult<LateralChain<ResidueAbbr>> {
    |_| todo!()
}

// =

/// Identifier = letter , { letter | digit | "_" } ;
fn identifier(i: &str) -> ParseResult<&str> {
    // PERF: Could maybe avoid allocations by using `many0_count` instead, but needs benchmarking
    let parser = recognize(pair(alpha1, many0(alt((alphanumeric1, tag("_"))))));
    wrap_err(parser, MuropeptideErrorKind::ExpectedIdentifier)(i)
}

// /// Any Modification = Named Modification | Offset Modification
// pub fn any_modification<'a, 'p, 's, K>(
//     polymer: &mut Polymer<'a, 'p>,
// ) -> impl FnMut(&'s str) -> ParseResult<ModificationId> {
//     alt((named_modification(polymer), offset_modification(polymer)))
// }

// // FIXME: I probably need to add a lot of `wrap_err`s around these parsers!
// /// Named Modification = [ Multiplier ] , Identifier
// pub fn named_modification<'a, 'p, 's, K>(
//     _polymer: &mut Polymer<'a, 'p>,
// ) -> impl FnMut(&'s str) -> ParseResult<ModificationId> {
//     |_| todo!()
// }

// /// Offset Modification = Offset Kind , [ Multiplier ] ,
// ///   Chemical Composition ;
// pub fn offset_modification<'a, 's, K>(
//     _polymer: &mut Polymer<'a, '_>,
// ) -> impl FnMut(&'s str) -> ParseResult<ModificationId> {
//     |_| todo!()
// }

/// Multiplier = Count , "x" ;
fn multiplier(i: &str) -> ParseResult<Count> {
    let mut parser = terminated(count, char('x'));
    // FIXME: Add error handling / reporting!
    parser(i)
}

/// uppercase
///   = "A" | "B" | "C" | "D" | "E" | "F" | "G"
///   | "H" | "I" | "J" | "K" | "L" | "M" | "N"
///   | "O" | "P" | "Q" | "R" | "S" | "T" | "U"
///   | "V" | "W" | "X" | "Y" | "Z"
///   ;
pub(crate) fn uppercase(i: &str) -> ParseResult<char> {
    let parser = satisfy(|c| c.is_ascii_uppercase());
    expect(parser, MuropeptideErrorKind::ExpectedUppercase)(i)
}

/// lowercase
///   = "a" | "b" | "c" | "d" | "e" | "f" | "g"
///   | "h" | "i" | "j" | "k" | "l" | "m" | "n"
///   | "o" | "p" | "q" | "r" | "s" | "t" | "u"
///   | "v" | "w" | "x" | "y" | "z"
///   ;
pub(crate) fn lowercase(i: &str) -> ParseResult<char> {
    let parser = satisfy(|c| c.is_ascii_lowercase());
    expect(parser, MuropeptideErrorKind::ExpectedLowercase)(i)
}

/// Count = digit - "0" , { digit } ;
pub(crate) fn count(i: &str) -> ParseResult<Count> {
    let not_zero = expect(
        cut(not(char('0'))),
        MuropeptideErrorKind::ExpectedNoLeadingZero,
    );
    let digits = expect(u32, MuropeptideErrorKind::ExpectedDigit);
    map(preceded(not_zero, digits), |c| Count::new(c).unwrap())(i)
}

type ParseResult<'a, O> = IResult<&'a str, O, LabeledParseError<'a, MuropeptideErrorKind>>;

#[derive(Clone, Eq, PartialEq, Debug, Diagnostic, Error)]
pub enum MuropeptideErrorKind {
    #[error("expected an ASCII letter, optionally followed by any number of ASCII letters, digits, and underscores")]
    ExpectedIdentifier,

    #[diagnostic(transparent)]
    #[error(transparent)]
    CompositionError(#[from] PolychemErrorKind),

    // FIXME: Update this help message to talk about PG structures, not chemical compositions!
    #[diagnostic(help(
        "a 0 value doesn't make sense here, if you've mistakenly included a leading zero, like \
        NH02, try just NH2 instead"
    ))]
    #[error("counts cannot start with 0")]
    ExpectedNoLeadingZero,

    #[error("expected an ASCII digit 1-9")]
    ExpectedDigit,

    #[error("expected an uppercase ASCII letter")]
    ExpectedUppercase,

    #[error("expected a lowercase ASCII letter")]
    ExpectedLowercase,

    #[diagnostic(help(
        "this is an internal error that you shouldn't ever see! If you have gotten this error, \
        then please report it as a bug!"
    ))]
    #[error("internal `nom` error: {0:?}")]
    NomError(ErrorKind),

    #[diagnostic(help(
        "check the unparsed region for errors, or remove it from the rest of the muropeptide"
    ))]
    #[error("could not interpret the full input as a valid muropeptide structure")]
    Incomplete,
}

impl LabeledErrorKind for MuropeptideErrorKind {
    fn label(&self) -> Option<&'static str> {
        Some(match self {
            // FIXME: Need to add branches for passing labels through the transparent errors?
            Self::Incomplete => "input was valid up until this point",
            Self::NomError(_) => "the region that triggered this bug!",
            _ => return None,
        })
    }
}

impl From<ErrorKind> for MuropeptideErrorKind {
    fn from(value: ErrorKind) -> Self {
        match value {
            ErrorKind::Eof => Self::Incomplete,
            kind => Self::NomError(kind),
        }
    }
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
        assert_eq!(count("1"), Ok(("", Count::new(1).unwrap())));
        assert_eq!(count("10"), Ok(("", Count::new(10).unwrap())));
        assert_eq!(count("422"), Ok(("", Count::new(422).unwrap())));
        assert_eq!(count("9999"), Ok(("", Count::new(9999).unwrap())));
        // Invalid Counts
        assert!(count("0").is_err());
        assert!(count("01").is_err());
        assert!(count("00145").is_err());
        assert!(count("H").is_err());
        assert!(count("p").is_err());
        assert!(count("+H").is_err());
        assert!(count("[H]").is_err());
        // Multiple Counts
        assert_eq!(count("1OH"), Ok(("OH", Count::new(1).unwrap())));
        assert_eq!(count("42HeH"), Ok(("HeH", Count::new(42).unwrap())));
    }

    #[test]
    fn test_identifier() {
        // Valid Identifiers
        assert_eq!(identifier("Ac"), Ok(("", "Ac")));
        assert_eq!(identifier("H2O"), Ok(("", "H2O")));
        assert_eq!(identifier("Anh"), Ok(("", "Anh")));
        assert_eq!(identifier("E2E"), Ok(("", "E2E")));
        assert_eq!(identifier("no_way"), Ok(("", "no_way")));
        assert_eq!(identifier("H"), Ok(("", "H")));
        assert_eq!(identifier("p"), Ok(("", "p")));
        // Invalid Identifiers
        assert!(identifier(" H2O").is_err());
        assert!(identifier("1").is_err());
        assert!(identifier("9999").is_err());
        assert!(identifier("0").is_err());
        assert!(identifier("00145").is_err());
        assert!(identifier("+H").is_err());
        assert!(identifier("[H]").is_err());
        assert!(identifier("Øof").is_err());
        assert!(identifier("2xAc").is_err());
        assert!(identifier("-Ac").is_err());
        assert!(identifier("_Ac").is_err());
        // Multiple Identifiers
        assert_eq!(identifier("OH-p"), Ok(("-p", "OH")));
        assert_eq!(identifier("HeH 2slow"), Ok((" 2slow", "HeH")));
        assert_eq!(identifier("Gefählt"), Ok(("ählt", "Gef")));
        // This is a weird unicode 6
        assert!('𝟨'.is_numeric());
        assert!(!'𝟨'.is_ascii_digit());
        assert_eq!(identifier("C2H𝟨O"), Ok(("𝟨O", "C2H")));
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

    // FIXME: Unfininshed! Needs modification support — same with unbranched_amino_acid!
    #[test]
    fn test_monosaccharide() {
        // Valid Monosaccharides
        assert_eq!(monosaccharide("g"), Ok(("", "g")));
        assert_eq!(monosaccharide("m"), Ok(("", "m")));
        // Invalid Monosaccharides
        assert!(monosaccharide("P").is_err());
        assert!(monosaccharide("EP").is_err());
        assert!(monosaccharide("1h").is_err());
        assert!(monosaccharide("+m").is_err());
        assert!(monosaccharide("-g").is_err());
        assert!(monosaccharide("[h]").is_err());
        // Multiple Monosaccharides
        assert_eq!(monosaccharide("gm"), Ok(("m", "g")));
        assert_eq!(monosaccharide("m-A"), Ok(("-A", "m")));
    }

    // FIXME: Unfininshed! Needs modification support — same with monosaccharide!
    #[test]
    fn test_unbranched_amino_acid() {
        // Valid Unbranched Amino Acids
        assert_eq!(unbranched_amino_acid("A"), Ok(("", "A")));
        assert_eq!(unbranched_amino_acid("E"), Ok(("", "E")));
        assert_eq!(unbranched_amino_acid("J"), Ok(("", "J")));
        assert_eq!(unbranched_amino_acid("yE"), Ok(("", "yE")));
        assert_eq!(unbranched_amino_acid("eK"), Ok(("", "eK")));
        // Invalid Unbranched Amino Acids
        assert!(unbranched_amino_acid("p").is_err());
        assert!(unbranched_amino_acid("eeP").is_err());
        assert!(unbranched_amino_acid("1H").is_err());
        assert!(unbranched_amino_acid("+M").is_err());
        assert!(unbranched_amino_acid("-G").is_err());
        assert!(unbranched_amino_acid("[H]").is_err());
        // Multiple Unbranched Amino Acids
        assert_eq!(unbranched_amino_acid("AEJA"), Ok(("EJA", "A")));
        assert_eq!(unbranched_amino_acid("EJA"), Ok(("JA", "E")));
        assert_eq!(unbranched_amino_acid("JA"), Ok(("A", "J")));
        assert_eq!(unbranched_amino_acid("yEJA"), Ok(("JA", "yE")));
        assert_eq!(unbranched_amino_acid("eK[GGGGG]"), Ok(("[GGGGG]", "eK")));
    }

    // FIXME: Add modification testing!
    // FIXME: Add lateral chain testing!
    #[test]
    fn test_peptide() {
        macro_rules! assert_peptide {
            ($input:literal, $output:literal, $residues:expr) => {
                let (rest, parsed_abbrs) = peptide($input).unwrap();
                assert_eq!(rest, $output);

                let residues = Vec::from($residues);
                assert_eq!(parsed_abbrs, residues);
            };
        }

        // Valid Peptides
        assert_peptide!("AEJA", "", ["A", "E", "J", "A"]);
        assert_peptide!("AyEJA", "", ["A", "yE", "J", "A"]);
        assert_peptide!("AE", "", ["A", "E"]);
        assert_peptide!("A", "", ["A"]);
        // Invalid Peptides
        assert!(peptide("y").is_err());
        assert!(peptide("yrE").is_err());
        assert!(peptide("-AEJA").is_err());
        assert!(peptide("[GGGGG]").is_err());
        assert!(peptide("gm-AEJA").is_err());
        assert!(peptide("(Am)").is_err());
        // Multiple Peptides
        assert_peptide!("AE=gm-AEJ", "=gm-AEJ", ["A", "E"]);
        assert_peptide!("AeeK", "eeK", ["A"]);
    }

    // FIXME: Add modification testing!
    #[test]
    fn test_glycan() {
        macro_rules! assert_glycan {
            ($input:literal, $output:literal, $residues:expr) => {
                let (rest, parsed_abbrs) = glycan($input).unwrap();
                assert_eq!(rest, $output);

                let residues = Vec::from($residues);
                assert_eq!(parsed_abbrs, residues);
            };
        }

        // Valid Glycans
        assert_glycan!("gmgm", "", ["g", "m", "g", "m"]);
        assert_glycan!("gm", "", ["g", "m"]);
        assert_glycan!("g", "", ["g"]);
        assert_glycan!("m", "", ["m"]);
        // Invalid Glycans
        assert!(glycan("Y").is_err());
        assert!(glycan("Ygm").is_err());
        assert!(glycan("-AEJA").is_err());
        assert!(glycan("[GGGGG]").is_err());
        assert!(glycan("EA=gm-AEJA").is_err());
        assert!(glycan("(Am)").is_err());
        // Multiple Glycans
        assert_glycan!("gm-AEJ", "-AEJ", ["g", "m"]);
        assert_glycan!("xAJgmK", "AJgmK", ["x"]);
    }

    #[test]
    fn test_monomer() {
        macro_rules! assert_monomer {
            ($input:literal, $output:literal, $glycan:expr, $peptide:expr) => {
                let (rest, parsed_monomer) = monomer($input).unwrap();
                assert_eq!(rest, $output);

                let glycan = Vec::<&str>::from($glycan);
                let peptide = Vec::<&str>::from($peptide);
                assert_eq!(parsed_monomer.glycan, glycan);
                assert_eq!(parsed_monomer.peptide, peptide);
            };
        }

        // Valid Monomers
        assert_monomer!("gmgm", "", ["g", "m", "g", "m"], []);
        assert_monomer!("AEJA", "", [], ["A", "E", "J", "A"]);
        assert_monomer!("gm-AE", "", ["g", "m"], ["A", "E"]);
        assert_monomer!("m-AyE", "", ["m"], ["A", "yE"]);
        assert_monomer!("m-A", "", ["m"], ["A"]);
        assert_monomer!("mA", "", [], ["mA"]);
        // Invalid Monomers
        assert!(monomer("+").is_err());
        assert!(monomer("~gm-").is_err());
        assert!(monomer("-AEJA").is_err());
        assert!(monomer("[GGGGG]").is_err());
        assert!(monomer("=gm-AEJA").is_err());
        assert!(monomer("(Am)").is_err());
        // Multiple Monomers
        assert_monomer!("gm-", "-", ["g", "m"], []);
        assert_monomer!("gm-AQK=gm-AQKAA", "=gm-AQKAA", ["g", "m"], ["A", "Q", "K"]);
        assert_monomer!("gm-AQK~gm-AQK", "~gm-AQK", ["g", "m"], ["A", "Q", "K"]);
        assert_monomer!("gm-AE (Am)", " (Am)", ["g", "m"], ["A", "E"]);
        assert_monomer!("xAJgmK", "gmK", [], ["xA", "J"]);
    }

    #[ignore]
    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_modifications() {
        // TODO: Restore from git!
        todo!();
    }

    // FIXME: Add a test that checks all of the errors using `assert_miette_snapshot`! Maybe make that a crate?
}
