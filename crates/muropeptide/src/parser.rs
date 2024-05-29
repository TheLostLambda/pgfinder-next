use std::cell::RefCell;

use miette::Diagnostic;
use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::{alpha1, alphanumeric1, char, space0},
    combinator::{opt, recognize},
    error::ErrorKind,
    multi::{many0, many1, separated_list1},
    sequence::{delimited, pair, preceded},
    IResult,
};
use nom_miette::{map_res, wrap_err, FromExternalError, LabeledErrorKind, LabeledParseError};
use polychem::{
    errors::PolychemError,
    parsers::{
        errors::PolychemErrorKind,
        modifications,
        primitives::{lowercase, uppercase},
    },
    ModificationId, Polymer,
};
use thiserror::Error;

use crate::{AminoAcid, LateralChain, Monomer, Monosaccharide, UnbranchedAminoAcid};

// FIXME: Need to think about if these should really live in another KDL config?
const PEPTIDE_BOND: &str = "Pep";
const GLYCOSIDIC_BOND: &str = "Gly";

// FIXME: A horrible hack that's needed to specify the lifetimes captured by `impl FnMut(...) -> ...` correctly. Once
// Rust 2024 is stabilized, however, this hack can be removed. Keep an eye on:
// https://github.com/rust-lang/rust/issues/117587
trait Captures<U> {}
impl<T: ?Sized, U> Captures<U> for T {}

// FIXME: Paste all of these EBNF comments into another file and make sure they are valid!
/// Monomer = Glycan , [ "-" , Peptide ] | Peptide ;
fn monomer<'a, 'p, 's>(
    _polymer: &mut Polymer<'a, 'p>,
) -> impl FnMut(&'s str) -> ParseResult<Monomer> {
    |_| todo!()
}

// =

/// Glycan = { Monosaccharide }- ;
fn glycan<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<Vec<Monosaccharide>> + Captures<(&'c (), &'a (), &'p ())> {
    let parser = many1(monosaccharide(polymer));
    map_res(parser, |residues| {
        let _ = polymer
            .borrow_mut()
            .bond_chain(GLYCOSIDIC_BOND, &residues)?;
        Ok(residues)
    })
}

// FIXME: This is using the wrong amino acid parser — needs lateral chain support!
/// Peptide = { Amino Acid }- ;
fn peptide<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<Vec<UnbranchedAminoAcid>> + Captures<(&'c (), &'a (), &'p ())>
{
    // FIXME: Change to branched amino acid!
    let parser = many1(unbranched_amino_acid(polymer));
    map_res(parser, |residues| {
        let _ = polymer.borrow_mut().bond_chain(PEPTIDE_BOND, &residues)?;
        Ok(residues)
    })
}

// =

// FIXME: Add modifications
/// Monosaccharide = lowercase , [ Modifications ] ;
fn monosaccharide<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<Monosaccharide> + Captures<(&'c (), &'a (), &'p ())> {
    let parser = recognize(lowercase);
    map_res(parser, |abbr| polymer.borrow_mut().new_residue(abbr))
}

/// Amino Acid = Unbranched Amino Acid , [ Lateral Chain ] ;
fn amino_acid<'a, 'p, 's>(
    _polymer: &mut Polymer<'a, 'p>,
) -> impl FnMut(&'s str) -> ParseResult<AminoAcid> {
    |_| todo!()
}

// =

/// Modifications = "(" , Any Modification ,
///   { { " " } , "," , { " " } , Any Modification } , ")" ;
fn modifications<'a, 'p, 's>(
    polymer: &mut Polymer<'a, 'p>,
) -> impl FnMut(&'s str) -> ParseResult<Vec<ModificationId>> {
    let separator = delimited(space0, char(','), space0);
    delimited(
        char('('),
        separated_list1(separator, modifications::any(polymer, identifier)),
        char(')'),
    )
}

// FIXME: Add modifications
/// Unbranched Amino Acid = [ lowercase ] , uppercase , [ Modifications ] ;
fn unbranched_amino_acid<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<UnbranchedAminoAcid> + Captures<(&'c (), &'a (), &'p ())> {
    let parser = recognize(preceded(opt(lowercase), uppercase));
    map_res(parser, |abbr| polymer.borrow_mut().new_residue(abbr))
}

// NOTE: These are not meant to be links, it's just EBNF
#[allow(clippy::doc_link_with_quotes)]
/// Lateral Chain = "[" , [ "<" (* C-to-N *) | ">" (* N-to-C *) ] ,
///   { Unbranched Amino Acid }- , "]" ;
fn lateral_chain<'a, 'p, 's>(
    _polymer: &mut Polymer<'a, 'p>,
) -> impl FnMut(&'s str) -> ParseResult<LateralChain> {
    |_| todo!()
}

// =

/// Identifier = letter , { letter | digit | "_" } ;
fn identifier(i: &str) -> ParseResult<&str> {
    // PERF: Could maybe avoid allocations by using `many0_count` instead, but needs benchmarking
    let parser = recognize(pair(alpha1, many0(alt((alphanumeric1, tag("_"))))));
    wrap_err(parser, MuropeptideErrorKind::ExpectedIdentifier)(i)
}

type ParseResult<'a, O> = IResult<&'a str, O, LabeledParseError<'a, MuropeptideErrorKind>>;

#[derive(Clone, Eq, PartialEq, Debug, Diagnostic, Error)]
pub enum MuropeptideErrorKind {
    #[error("expected an ASCII letter, optionally followed by any number of ASCII letters, digits, and underscores")]
    ExpectedIdentifier,

    // FIXME: Kill this and merge into the error below!
    #[diagnostic(transparent)]
    #[error(transparent)]
    PolychemError(Box<PolychemError>),

    #[diagnostic(transparent)]
    #[error(transparent)]
    CompositionError(#[from] PolychemErrorKind),

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

// FIXME: Can I get rid of this?
impl<'a> FromExternalError<'a, Box<PolychemError>> for MuropeptideErrorKind {
    const FATAL: bool = true;

    fn from_external_error(input: &'a str, e: Box<PolychemError>) -> LabeledParseError<'_, Self> {
        LabeledParseError::new(input, Self::PolychemError(e))
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
    use once_cell::sync::Lazy;
    use polychem::{AtomicDatabase, Massive, PolymerDatabase, Polymerizer};
    use rust_decimal::Decimal;
    use rust_decimal_macros::dec;

    use super::*;

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);
    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "polymer_database.kdl",
            include_str!("../tests/data/polymer_database.kdl"),
        )
        .unwrap()
    });

    static POLYMERIZER: Lazy<Polymerizer> = Lazy::new(|| Polymerizer::new(&ATOMIC_DB, &POLYMER_DB));

    #[ignore]
    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_modifications() {
        // TODO: Restore from git!
        todo!();
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

    // FIXME: Unfininshed! Needs modification support — same with unbranched_amino_acid!
    #[test]
    fn test_monosaccharide() {
        let polymer = RefCell::new(POLYMERIZER.new_polymer());

        let mut monosaccharide = monosaccharide(&polymer);
        macro_rules! assert_monosaccharide_name {
            ($input:literal, $output:literal, $name:literal) => {
                let (rest, id) = monosaccharide($input).unwrap();
                assert_eq!(
                    (rest, polymer.borrow().residue(id).unwrap().name()),
                    ($output, $name)
                );
            };
        }

        // Valid Monosaccharides
        assert_monosaccharide_name!("g", "", "N-Acetylglucosamine");
        assert_monosaccharide_name!("m", "", "N-Acetylmuramic Acid");
        // Invalid Monosaccharides
        assert!(monosaccharide("P").is_err());
        assert!(monosaccharide("EP").is_err());
        assert!(monosaccharide("1h").is_err());
        assert!(monosaccharide("+m").is_err());
        assert!(monosaccharide("-g").is_err());
        assert!(monosaccharide("[h]").is_err());
        // Non-Existent Monosaccharides
        assert!(monosaccharide("s").is_err());
        assert!(monosaccharide("f").is_err());
        // Multiple Monosaccharides
        assert_monosaccharide_name!("gm", "m", "N-Acetylglucosamine");
        assert_monosaccharide_name!("m-A", "-A", "N-Acetylmuramic Acid");
    }

    // FIXME: Unfininshed! Needs modification support — same with monosaccharide!
    #[test]
    fn test_unbranched_amino_acid() {
        let polymer = RefCell::new(POLYMERIZER.new_polymer());

        let mut unbranched_amino_acid = unbranched_amino_acid(&polymer);
        macro_rules! assert_unbranched_aa_name {
            ($input:literal, $output:literal, $name:literal) => {
                let (rest, id) = unbranched_amino_acid($input).unwrap();
                assert_eq!(
                    (rest, polymer.borrow().residue(id).unwrap().name()),
                    ($output, $name)
                );
            };
        }

        // Valid Unbranched Amino Acids
        assert_unbranched_aa_name!("A", "", "Alanine");
        assert_unbranched_aa_name!("E", "", "Glutamic Acid");
        assert_unbranched_aa_name!("J", "", "Diaminopimelic Acid");
        assert_unbranched_aa_name!("yE", "", "γ-Glutamate");
        assert_unbranched_aa_name!("eK", "", "ε-Lysine");
        // Invalid Unbranched Amino Acids
        assert!(unbranched_amino_acid("p").is_err());
        assert!(unbranched_amino_acid("eP").is_err());
        assert!(unbranched_amino_acid("1H").is_err());
        assert!(unbranched_amino_acid("+M").is_err());
        assert!(unbranched_amino_acid("-G").is_err());
        assert!(unbranched_amino_acid("[H]").is_err());
        // Non-Existent Unbranched Amino Acids
        assert!(unbranched_amino_acid("iA").is_err());
        assert!(unbranched_amino_acid("yK").is_err());
        // Multiple Unbranched Amino Acids
        assert_unbranched_aa_name!("AEJA", "EJA", "Alanine");
        assert_unbranched_aa_name!("EJA", "JA", "Glutamic Acid");
        assert_unbranched_aa_name!("JA", "A", "Diaminopimelic Acid");
        assert_unbranched_aa_name!("yEJA", "JA", "γ-Glutamate");
        assert_unbranched_aa_name!("eK[GGGGG]", "[GGGGG]", "ε-Lysine");
    }

    // FIXME: Add modification testing!
    // FIXME: Add lateral chain testing!
    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_peptide() {
        let polymer = RefCell::new(POLYMERIZER.new_polymer());

        let mut err_peptide = peptide(&polymer);
        macro_rules! assert_chain_residues_and_masses {
            ($input:literal, $output:literal, $residues:expr, $mono_mass:literal, $avg_mass:literal) => {
                let polymer = RefCell::new(POLYMERIZER.new_polymer());

                let (rest, parsed_ids) = peptide(&polymer)($input).unwrap();
                assert_eq!(rest, $output);

                let polymer = polymer.borrow();
                let parsed_ids: Vec<_> = parsed_ids
                    .into_iter()
                    .map(|id| polymer.residue(id).unwrap().name())
                    .collect();
                let residues = Vec::from($residues);
                assert_eq!(parsed_ids, residues);

                assert_eq!(Decimal::from(polymer.monoisotopic_mass()), dec!($mono_mass));
                assert_eq!(Decimal::from(polymer.average_mass()), dec!($avg_mass));
            };
        }

        // Valid Peptides
        assert_chain_residues_and_masses!(
            "AEJA",
            "",
            ["Alanine", "Glutamic Acid", "Diaminopimelic Acid", "Alanine"],
            461.21217759741,
            461.46756989305707095
        );
        assert_chain_residues_and_masses!(
            "AyEJA",
            "",
            ["Alanine", "γ-Glutamate", "Diaminopimelic Acid", "Alanine"],
            461.21217759741,
            461.46756989305707095
        );
        assert_chain_residues_and_masses!(
            "AE",
            "",
            ["Alanine", "Glutamic Acid"],
            218.09027155793,
            218.20748877514586040
        );
        assert_chain_residues_and_masses!(
            "A",
            "",
            ["Alanine"],
            89.04767846918,
            89.09330602867854225
        );
        // Invalid Peptides
        assert!(err_peptide("y").is_err());
        assert!(err_peptide("yrE").is_err());
        assert!(err_peptide("-AEJA").is_err());
        assert!(err_peptide("[GGGGG]").is_err());
        assert!(err_peptide("gm-AEJA").is_err());
        assert!(err_peptide("(Am)").is_err());
        // Non-Existent Peptide Residues
        assert!(err_peptide("AEJiA").is_err());
        assert!(err_peptide("AQyK").is_err());
        // Multiple Peptides
        assert_chain_residues_and_masses!(
            "AE=gm-AEJ",
            "=gm-AEJ",
            ["Alanine", "Glutamic Acid"],
            218.09027155793,
            218.20748877514586040
        );
        assert_chain_residues_and_masses!(
            "AeeK",
            "eeK",
            ["Alanine"],
            89.04767846918,
            89.09330602867854225
        );
    }

    // FIXME: Add modification testing!
    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_glycan() {
        let polymer = RefCell::new(POLYMERIZER.new_polymer());

        let mut err_glycan = glycan(&polymer);
        macro_rules! assert_chain_residues_and_masses {
            ($input:literal, $output:literal, $residues:expr, $mono_mass:literal, $avg_mass:literal) => {
                let polymer = RefCell::new(POLYMERIZER.new_polymer());

                let (rest, parsed_ids) = glycan(&polymer)($input).unwrap();
                assert_eq!(rest, $output);

                let polymer = polymer.borrow();
                let parsed_ids: Vec<_> = parsed_ids
                    .into_iter()
                    .map(|id| polymer.residue(id).unwrap().name())
                    .collect();
                let residues = Vec::from($residues);
                assert_eq!(parsed_ids, residues);

                assert_eq!(Decimal::from(polymer.monoisotopic_mass()), dec!($mono_mass));
                assert_eq!(Decimal::from(polymer.average_mass()), dec!($avg_mass));
            };
        }

        // Valid Glycans
        assert_chain_residues_and_masses!(
            "gmgm",
            "",
            [
                "N-Acetylglucosamine",
                "N-Acetylmuramic Acid",
                "N-Acetylglucosamine",
                "N-Acetylmuramic Acid"
            ],
            974.37031350523,
            974.91222678113779720
        );
        assert_chain_residues_and_masses!(
            "gm",
            "",
            ["N-Acetylglucosamine", "N-Acetylmuramic Acid"],
            496.19043909463,
            496.46375660678381490
        );
        assert_chain_residues_and_masses!(
            "g",
            "",
            ["N-Acetylglucosamine"],
            221.08993720530,
            221.20813124207411765
        );
        assert_chain_residues_and_masses!(
            "m",
            "",
            ["N-Acetylmuramic Acid"],
            293.11106657336,
            293.27091179713952985
        );
        // Invalid Glycans
        assert!(err_glycan("Y").is_err());
        assert!(err_glycan("Ygm").is_err());
        assert!(err_glycan("-AEJA").is_err());
        assert!(err_glycan("[GGGGG]").is_err());
        assert!(err_glycan("EA=gm-AEJA").is_err());
        assert!(err_glycan("(Am)").is_err());
        // Non-Existent Glycan Residues
        assert!(err_glycan("y").is_err());
        assert!(err_glycan("fp").is_err());
        // Multiple Glycans
        assert_chain_residues_and_masses!(
            "gm-AEJ",
            "-AEJ",
            ["N-Acetylglucosamine", "N-Acetylmuramic Acid"],
            496.19043909463,
            496.46375660678381490
        );
        assert_chain_residues_and_masses!("xAJgmK", "AJgmK", ["Unknown Monosaccharide"], 0.0, 0.0);
    }
    // FIXME: Add a test that checks all of the errors using `assert_miette_snapshot`! Maybe make that a crate?
}
