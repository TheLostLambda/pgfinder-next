use miette::Diagnostic;
use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::{alpha1, alphanumeric1, char, space0},
    combinator::{map, opt, recognize},
    error::ErrorKind,
    multi::{many0, separated_list1},
    sequence::{delimited, pair, terminated, tuple},
    IResult,
};
use nom_miette::{into, map_res, wrap_err, FromExternalError, LabeledErrorKind, LabeledParseError};
use polychem::{
    atoms::chemical_composition::{self, CompositionErrorKind},
    polymerizer::Polymerizer,
    AnyModification, ChemicalComposition, Count, Modification, NamedMod, OffsetKind, OffsetMod,
    PolychemError,
};
use thiserror::Error;

use crate::{AminoAcid, LateralChain, Monomer, Monosaccharide, UnbranchedAminoAcid};

// FIXME: Paste all of these EBNF comments into another file and make sure they are valid!
/// Monomer = Glycan , [ "-" , Peptide ] | Peptide ;
fn monomer<'a, 's>(
    _polymerizer: &mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<Monomer<'a>> {
    |_| todo!()
}

// =

/// Glycan = { Monosaccharide }- ;
fn glycan<'a, 's>(
    _polymerizer: &mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<Vec<Monosaccharide<'a>>> {
    |_| todo!()
}

/// Peptide = { Amino Acid }- ;
fn peptide<'a, 's>(
    _polymerizer: &mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<Vec<AminoAcid<'a>>> {
    |_| todo!()
}

// =

// FIXME: Don't know if it's a great idea to tie together the lifetimes of the chemical databases and the polymerizer
// instance here? Everything is using 'a...
// FIXME: Add modifications
/// Monosaccharide = lowercase , [ Modifications ] ;
fn monosaccharide<'a, 's>(
    polymerizer: &'a mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<Monosaccharide<'a>> {
    let parser = recognize(lowercase);
    map_res(parser, |abbr| polymerizer.residue(abbr))
}

/// Amino Acid = Unbranched Amino Acid , [ Lateral Chain ] ;
fn amino_acid<'a, 's>(
    _polymerizer: &mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<AminoAcid<'a>> {
    |_| todo!()
}

// =

/// Modifications = "(" ,
///   ( Predefined Modification
///   | Chemical Offset
///   ) , { { " " } , "," , { " " } ,
///   ( Predefined Modification
///   | Chemical Offset
///   ) } , ")" ;
fn modifications<'a, 's>(
    polymerizer: &Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<Vec<AnyModification<'a, 'a>>> {
    let separator = delimited(space0, char(','), space0);
    let any_mod = alt((
        map(predefined_modification(polymerizer), Modification::convert),
        map(chemical_offset(polymerizer), Modification::convert),
    ));
    delimited(char('('), separated_list1(separator, any_mod), char(')'))
}

// FIXME: Make private again!
// FIXME: Switch to a more efficient modification application API
/// Unbranched Amino Acid = uppercase , [ Modifications ] ;
pub fn unbranched_amino_acid<'a, 's>(
    polymerizer: &'a mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<UnbranchedAminoAcid<'a>> {
    let parser = pair(recognize(uppercase), opt(modifications(polymerizer)));
    map_res(parser, |(abbr, modifications)| {
        let mut amino_acid = polymerizer.residue(abbr)?;
        for modification in modifications.into_iter().flatten() {
            polymerizer.apply_modification(modification, &mut amino_acid)?;
        }
        Ok(amino_acid)
    })
}

/// Lateral Chain = "[" , [ "<" (* C-to-N *) | ">" (* N-to-C *) ] ,
///   { Unbranched Amino Acid }- , "]" ;
fn lateral_chain<'a, 's>(
    _polymerizer: &mut Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<LateralChain> {
    |_| todo!()
}

// =

// FIXME: I probably need to add a lot of `wrap_err`s around these parsers!
/// Predefined Modification = [ Multiplier ] , Identifier
fn predefined_modification<'a, 's>(
    polymerizer: &Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<Modification<NamedMod<'a, 'a>>> {
    let polymer_db = polymerizer.polymer_db();
    let named_mod = map_res(identifier, |abbr| NamedMod::new(polymer_db, abbr));
    let parser = pair(opt(multiplier), named_mod);
    map(parser, |(multiplier, named_mod)| {
        Modification::new(multiplier.unwrap_or(1), named_mod)
    })
}

/// Chemical Offset = Offset Kind , [ Multiplier ] ,
///   Chemical Composition ;
fn chemical_offset<'a, 's>(
    polymerizer: &Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<Modification<OffsetMod<'a>>> {
    let parser = tuple((
        offset_kind,
        opt(multiplier),
        chemical_composition(polymerizer),
    ));
    map(parser, |(kind, multiplier, composition)| {
        Modification::new(
            multiplier.unwrap_or(1),
            OffsetMod::new_with_composition(kind, composition),
        )
    })
}

// =

/// Multiplier = Count , "x" ;
fn multiplier(i: &str) -> ParseResult<Count> {
    let parser = terminated(count, char('x'));
    wrap_err(parser, MuropeptideErrorKind::ExpectedMultiplier)(i)
}

/// Identifier = letter , { letter | digit | "_" } ;
fn identifier(i: &str) -> ParseResult<&str> {
    // PERF: Could maybe avoid allocations by using `many0_count` instead, but needs benchmarking
    let parser = recognize(pair(alpha1, many0(alt((alphanumeric1, tag("_"))))));
    wrap_err(parser, MuropeptideErrorKind::ExpectedIdentifier)(i)
}

// Adapted parsers =

fn chemical_composition<'a, 's>(
    polymerizer: &Polymerizer<'a, 'a>,
) -> impl FnMut(&'s str) -> ParseResult<ChemicalComposition<'a>> {
    into(chemical_composition::chemical_composition(
        polymerizer.atomic_db(),
    ))
}

macro_rules! wrap_composition_parsers {
    ($($f:ident -> $t:ty),+ $(,)?) => {
        $(
            fn $f(i: &str) -> ParseResult<$t> {
                into(chemical_composition::$f)(i)
            }
        )+
    };
}

wrap_composition_parsers!(
    count -> Count,
    offset_kind -> OffsetKind,
    uppercase -> char,
    lowercase -> char,
);

type ParseResult<'a, O> = IResult<&'a str, O, LabeledParseError<'a, MuropeptideErrorKind>>;

#[derive(Clone, Eq, PartialEq, Debug, Diagnostic, Error)]
pub enum MuropeptideErrorKind {
    #[error("expected an ASCII letter, optionally followed by any number of ASCII letters, digits, and underscores")]
    ExpectedIdentifier,

    #[error("expected a count followed by 'x'")]
    ExpectedMultiplier,

    #[diagnostic(transparent)]
    #[error(transparent)]
    PolychemError(Box<PolychemError>),

    #[diagnostic(transparent)]
    #[error(transparent)]
    CompositionError(#[from] CompositionErrorKind),

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
    use polychem::{AtomicDatabase, Charge, Charged, Massive, PolymerDatabase};
    use rust_decimal::Decimal;
    use rust_decimal_macros::dec;

    use super::*;

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(|| {
        AtomicDatabase::from_kdl(
            "atomic_database.kdl",
            include_str!("../../polychem/data/atomic_database.kdl"),
        )
        .unwrap()
    });

    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::from_kdl(
            &ATOMIC_DB,
            "polymer_database.kdl",
            include_str!("../data/polymer_database.kdl"),
        )
        .unwrap()
    });

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
    // NOTE: Complexity comes from the macro confusing clippy — converting to a function, however, adds real complexity
    // since borrow-checker issues crop up when writing a closure that captures `chemical_offset`
    #[allow(clippy::cognitive_complexity)]
    fn test_chemical_offset() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut chemical_offset = chemical_offset(&mut polymerizer);
        macro_rules! assert_offset_mz {
            ($input:literal, $output:literal, $mass:expr, $charge:literal) => {
                let (rest, modification) = chemical_offset($input).unwrap();
                assert_eq!(rest, $output);
                assert_eq!(modification.monoisotopic_mass(), $mass);
                assert_eq!(modification.charge(), $charge);
            };
        }
        // Valid Chemical Offsets
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
        // Invalid Chemical Offsets
        assert!(chemical_offset(" ").is_err());
        assert!(chemical_offset("H2O").is_err());
        assert!(chemical_offset("(-H2O)").is_err());
        assert!(chemical_offset("+0xH2O").is_err());
        assert!(chemical_offset("2xH2O").is_err());
        assert!(chemical_offset("-2x3xH2O").is_err());
        assert!(chemical_offset("-2x+H2O").is_err());
        assert!(chemical_offset("+2[2H]").is_err());
        assert!(chemical_offset("-[H+p]O").is_err());
        assert!(chemical_offset("+NH2[100Tc]").is_err());
        // Multiple Chemical Offsets
        assert_offset_mz!("+[37Cl]5-2p10", "10", dec!(182.814960076758), -2);
        assert_offset_mz!("+[2H]2O*H2O", "*H2O", dec!(20.02311817581), 0);
        assert_offset_mz!("+NH2{100Tc", "{100Tc", dec!(16.01872406889), 0);
        assert_offset_mz!("+C11H12N2O2 H2O", " H2O", dec!(204.08987763476), 0);
    }

    #[test]
    fn test_predefined_modification() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut predefined_modification = predefined_modification(&mut polymerizer);
        macro_rules! assert_offset_mass {
            ($input:literal, $output:literal, $mass:expr) => {
                let (rest, modification) = predefined_modification($input).unwrap();
                assert_eq!(rest, $output);
                assert_eq!(modification.monoisotopic_mass(), $mass);
            };
        }
        // Valid Predefined Modifications
        assert_offset_mass!("Am", "", dec!(-0.98401558291));
        assert_offset_mass!("Ac", "", dec!(42.01056468403));
        assert_offset_mass!("Poly", "", dec!(77.95068082490));
        assert_offset_mass!("DeAc", "", dec!(-42.01056468403));
        assert_offset_mass!("Red", "", dec!(2.01565006446));
        assert_offset_mass!("Anh", "", dec!(-18.01056468403));
        assert_offset_mass!("1xAm", "", dec!(-0.98401558291));
        assert_offset_mass!("2xRed", "", dec!(4.03130012892));
        assert_offset_mass!("3xAnh", "", dec!(-54.03169405209));
        // Invalid Predefined Modifications
        assert!(predefined_modification(" H2O").is_err());
        assert!(predefined_modification("1").is_err());
        assert!(predefined_modification("9999").is_err());
        assert!(predefined_modification("0").is_err());
        assert!(predefined_modification("00145").is_err());
        assert!(predefined_modification("+H").is_err());
        assert!(predefined_modification("[H]").is_err());
        assert!(predefined_modification("Øof").is_err());
        assert!(predefined_modification("-Ac").is_err());
        assert!(predefined_modification("_Ac").is_err());
        assert!(predefined_modification("+Am").is_err());
        assert!(predefined_modification("-2xAm").is_err());
        assert!(predefined_modification("(Am)").is_err());
        assert!(predefined_modification("-4xH2O").is_err());
        assert!(predefined_modification("-2p").is_err());
        assert!(predefined_modification("+C2H2O-2e").is_err());
        assert!(predefined_modification("-3xC2H2O-2e").is_err());
        assert!(predefined_modification("+NH3+p").is_err());
        assert!(predefined_modification("+2xD2O").is_err());
        assert!(predefined_modification("-2x[2H]2O").is_err());
        // Non-Existent Predefined Modifications
        assert!(predefined_modification("Blue").is_err());
        assert!(predefined_modification("Hydro").is_err());
        assert!(predefined_modification("1xAm2").is_err());
        assert!(predefined_modification("2xR_ed").is_err());
        // Multiple Predefined Modifications
        assert_offset_mass!("Anh, Am", ", Am", dec!(-18.01056468403));
        assert_offset_mass!("1xAm)JAA", ")JAA", dec!(-0.98401558291));
    }

    #[test]
    fn test_modifications() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut modifications = modifications(&mut polymerizer);
        macro_rules! assert_offset_mz {
            ($input:literal, $output:literal, $mass:expr, $charge:literal) => {
                let (rest, modification) = modifications($input).unwrap();
                assert_eq!(rest, $output);
                let net_mass: Decimal = modification.iter().map(|m| m.monoisotopic_mass()).sum();
                assert_eq!(net_mass, $mass);
                let net_charge: Charge = modification.iter().map(|m| m.charge()).sum();
                assert_eq!(net_charge, $charge);
            };
        }
        // Valid Modifications
        assert_offset_mz!("(-H2O)", "", dec!(-18.01056468403), 0);
        assert_offset_mz!("(+2xH2O)", "", dec!(36.02112936806), 0);
        assert_offset_mz!("(-2p)", "", dec!(-2.014552933242), -2);
        assert_offset_mz!("(-3xC2H2O-2e)", "", dec!(-126.02840257263561), -6);
        assert_offset_mz!("(+[37Cl]5-2p)", "", dec!(182.814960076758), -2);
        assert_offset_mz!("(Red)", "", dec!(2.01565006446), 0);
        assert_offset_mz!("(Anh)", "", dec!(-18.01056468403), 0);
        assert_offset_mz!("(1xAm)", "", dec!(-0.98401558291), 0);
        assert_offset_mz!("(2xRed)", "", dec!(4.03130012892), 0);
        assert_offset_mz!("(-OH, +NH2)", "", dec!(-0.98401558291), 0);
        assert_offset_mz!("(-H2, +Ca)", "", dec!(37.94694079854), 0);
        assert_offset_mz!("(Anh, +H2O)", "", dec!(0), 0);
        assert_offset_mz!("(Anh,+H2O)", "", dec!(0), 0);
        assert_offset_mz!("(Anh   ,+H2O)", "", dec!(0), 0);
        assert_offset_mz!("(Anh  ,  +H2O)", "", dec!(0), 0);
        assert_offset_mz!("(2xAnh, +3xH2O)", "", dec!(18.01056468403), 0);
        assert_offset_mz!("(Anh, Anh, +3xH2O)", "", dec!(18.01056468403), 0);
        // There is a super small mass defect (13.6 eV, or ~1e-8 u) stored in the binding energy between a proton and
        // electon — that's why this result is slightly different from the one above!
        assert_offset_mz!("(-2p, +Ca-2e)", "", dec!(37.946940769939870), 0);
        assert_offset_mz!("(+2p, -2p, +Ca-2e)", "", dec!(39.961493703181870), 2);
        // Invalid Modifications
        assert!(modifications(" ").is_err());
        assert!(modifications("H2O").is_err());
        assert!(modifications("(-H2O").is_err());
        assert!(modifications("(+0xH2O)").is_err());
        assert!(modifications("(2xH2O)").is_err());
        assert!(modifications("(-2x3xH2O)").is_err());
        assert!(modifications("(-2x+H2O)").is_err());
        assert!(modifications("(+2[2H])").is_err());
        assert!(modifications("(-[H+p]O)").is_err());
        assert!(modifications("(+NH2[100Tc])").is_err());
        assert!(modifications("( H2O)").is_err());
        assert!(modifications("(1)").is_err());
        assert!(modifications("(9999)").is_err());
        assert!(modifications("(0)").is_err());
        assert!(modifications("(00145)").is_err());
        assert!(modifications("([H])").is_err());
        assert!(modifications("(Øof)").is_err());
        assert!(modifications("(-Ac)").is_err());
        assert!(modifications("(_Ac)").is_err());
        assert!(modifications("(+Am)").is_err());
        assert!(modifications("(-2xAm)").is_err());
        assert!(modifications("((Am))").is_err());
        assert!(modifications("(Anh +H2O)").is_err());
        assert!(modifications("(Anh; +H2O)").is_err());
        // Non-Existent Modifications
        assert!(modifications("(Blue)").is_err());
        assert!(modifications("(Hydro)").is_err());
        assert!(modifications("(1xAm2)").is_err());
        assert!(modifications("(2xR_ed)").is_err());
        // Multiple Modifications
        assert_offset_mz!("(+[37Cl]5-2p)10", "10", dec!(182.814960076758), -2);
        assert_offset_mz!("(+[2H]2O)*H2O", "*H2O", dec!(20.02311817581), 0);
        assert_offset_mz!("(+NH2){100Tc", "{100Tc", dec!(16.01872406889), 0);
        assert_offset_mz!("(+C11H12N2O2) H2O", " H2O", dec!(204.08987763476), 0);
        assert_offset_mz!(
            "(2xAnh, +3xH2O)AA=gm-AEJA",
            "AA=gm-AEJA",
            dec!(18.01056468403),
            0
        );
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
    fn test_monosaccharide() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut monosaccharide = monosaccharide(&mut polymerizer);
        macro_rules! assert_monosaccharide_name {
            ($input:literal, $output:literal, $name:literal) => {
                assert_eq!(
                    monosaccharide($input).map(|(r, e)| (r, e.name())),
                    Ok(($output, $name))
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
    #[ignore]
    #[test]
    fn test_unbranched_amino_acid() {
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut unbranched_amino_acid = unbranched_amino_acid(&mut polymerizer);
        macro_rules! assert_unbranched_aa_name {
            ($input:literal, $output:literal, $name:literal) => {
                assert_eq!(
                    unbranched_amino_acid($input).map(|(r, e)| (r, e.name())),
                    Ok(($output, $name))
                );
            };
        }
        // Valid Unbranched Amino Acids
        assert_unbranched_aa_name!("g", "", "N-Acetylglucosamine");
        assert_unbranched_aa_name!("m", "", "N-Acetylmuramic Acid");
        // Invalid Unbranched Amino Acids
        assert!(unbranched_amino_acid("P").is_err());
        assert!(unbranched_amino_acid("EP").is_err());
        assert!(unbranched_amino_acid("1h").is_err());
        assert!(unbranched_amino_acid("+m").is_err());
        assert!(unbranched_amino_acid("-g").is_err());
        assert!(unbranched_amino_acid("[h]").is_err());
        // Non-Existent Unbranched Amino Acids
        assert!(unbranched_amino_acid("s").is_err());
        assert!(unbranched_amino_acid("f").is_err());
        // Multiple Unbranched Amino Acids
        assert_unbranched_aa_name!("gm", "m", "N-Acetylglucosamine");
        assert_unbranched_aa_name!("m-A", "-A", "N-Acetylmuramic Acid");
    }

    // FIXME: Add a test that checks all of the errors using `assert_miette_snapshot`! Maybe make that a crate?
}
