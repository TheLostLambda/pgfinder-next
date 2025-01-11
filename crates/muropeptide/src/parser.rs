use std::{cell::RefCell, iter::zip};

use itertools::{EitherOrBoth, Itertools};
use miette::Diagnostic;
use nom::{
    IResult,
    branch::alt,
    bytes::complete::tag,
    character::complete::{alpha1, alphanumeric1, char, one_of, space0, space1},
    combinator::{cut, map, opt, recognize},
    error::ErrorKind,
    multi::{many0, many1, separated_list1},
    sequence::{delimited, pair, preceded, terminated, tuple},
};
use nom_miette::{FromExternalError, LabeledErrorKind, LabeledParseError, map_res, wrap_err};
use polychem::{
    Count, ModificationId, Polymer, Polymerizer,
    errors::PolychemError,
    parsers::{
        chemical_composition,
        errors::PolychemErrorKind,
        primitives::{count, lowercase, offset_kind, uppercase},
    },
};
use thiserror::Error;

use crate::{
    AminoAcid, Connection, CrosslinkDescriptor, CrosslinkDescriptors, LateralChain, Monomer,
    Monosaccharide, Muropeptide, PeptideDirection, Position, UnbranchedAminoAcid,
};

// FIXME: Need to think about if these should really live in another KDL config?
const PEPTIDE_BOND: &str = "Pep";
const GLYCOSIDIC_BOND: &str = "Gly";
const STEM_BOND: &str = "Stem";
const NTOC_BOND: &str = "NToC";
const CTON_BOND: &str = "CToN";
const CROSSLINK_BOND: &str = "Link";
const LAT_CROSSLINK_BOND: &str = "Lat-Link";

// FIXME: Paste all of these EBNF comments into another file and make sure they are valid!

/// Muropeptide = Monomer , { Connection , Monomer } , [ Connection ] , [ { " " }- ,
///   ( Modifications , [ { " " }- , Crosslinks ]
///   | Crosslinks , [ { " " }- , Modifications ]
///   ) ] ;
// FIXME: Very very incomplete!
// FIXME: Needs proper testing!
// FIXME: And disgustingly messy... Needs a major refactor... When you do, remove this allow!
#[allow(clippy::too_many_lines)]
pub fn muropeptide<'z, 'a, 'p, 's>(
    polymerizer: &'z Polymerizer<'a, 'p>,
) -> impl FnMut(&'s str) -> ParseResult<'s, Muropeptide<'a, 'p>> {
    move |i| {
        let polymer = RefCell::new(polymerizer.new_polymer());
        // FIXME: Perhaps there is a better way to shorten that `polymer` borrow...
        let (rest, (monomers, connections)) = {
            let multimer = map(
                tuple((
                    monomer(&polymer),
                    many0(pair(connection, monomer(&polymer))),
                )),
                |(monomer, multimers)| {
                    let mut monomers = Vec::with_capacity(1 + multimers.len());
                    // PERF: This often over-allocates by one element — is checking `circular_connection.is_some()`
                    // something that's worth the branch here? Need to benchmark that!
                    let mut connections = Vec::with_capacity(multimers.len() + 1);

                    monomers.push(monomer);
                    for (connection, monomer) in multimers {
                        connections.push(connection);
                        monomers.push(monomer);
                    }

                    (monomers, connections)
                },
            );
            let opt_crosslinks = opt(preceded(space1, crosslinks));

            let mut parser = map_res(
                pair(multimer, opt_crosslinks),
                |((monomers, mut connections), crosslinks)| {
                    let crosslink_connections =
                        connections.iter_mut().filter_map(|conn| match conn {
                            Connection::GlycosidicBond => None,
                            Connection::Crosslink(descriptors) | Connection::Both(descriptors) => {
                                Some(descriptors)
                            }
                        });
                    let crosslinks = crosslinks.into_iter().flatten();

                    for either_or_both in crosslink_connections.zip_longest(crosslinks) {
                        match either_or_both {
                            EitherOrBoth::Both(connection, crosslinks) => {
                                *connection = vec![crosslinks];
                            }
                            _ => return Err(ConstructionError::CrosslinkCountMismatch),
                        }
                    }

                    let mut last_glycan = None;
                    let mut last_peptide = None;
                    let monomer_pairs = monomers.iter().circular_tuple_windows();
                    let connection_refs = connections.iter();
                    for ((left, right), connection) in zip(monomer_pairs, connection_refs) {
                        if !left.glycan.is_empty() {
                            last_glycan = Some(&left.glycan);
                        }
                        if !left.peptide.is_empty() {
                            last_peptide = Some(&left.peptide);
                        }

                        let gly_bond = || -> Result<(), ConstructionError> {
                            let donor = *last_glycan
                                .and_then(|lg| lg.last())
                                .ok_or(ConstructionError::NoGlycan)?;
                            let acceptor =
                                *right.glycan.first().ok_or(ConstructionError::NoGlycan)?;
                            polymer
                                .borrow_mut()
                                .bond_residues(GLYCOSIDIC_BOND, donor, acceptor)?;
                            Ok(())
                        };
                        let crosslink = |descriptors| -> Result<(), ConstructionError> {
                            fn lookup_amino_acid(
                                peptide: &[AminoAcid],
                                idx: u8,
                            ) -> Result<&AminoAcid, ConstructionError> {
                                peptide
                                    .get(idx.saturating_sub(1) as usize)
                                    .ok_or(ConstructionError::NoStemResidue)
                            }
                            for &descriptor in descriptors {
                                // FIXME: A nasty, very lazy hack:
                                let mut lateral_acceptor = false;
                                // FIXME: Are lateral chains ever donors?
                                let lookup_donor = |peptide, idx| {
                                    lookup_amino_acid(peptide, idx).map(|aa| aa.residue)
                                };
                                // FIXME: Good lord... What a mess...
                                let mut lookup_acceptor =
                                    |peptide, idx| -> Result<_, ConstructionError> {
                                        let AminoAcid {
                                            residue,
                                            lateral_chain,
                                        } = lookup_amino_acid(peptide, idx)?;
                                        Ok(*if let Some(LateralChain { peptide, .. }) =
                                            lateral_chain
                                        {
                                            lateral_acceptor = true;
                                            // SAFETY: All `LateralChain`s carry non-empty vectors, so this call to
                                            // `.last()` should never fail!
                                            peptide.last().unwrap()
                                        } else {
                                            residue
                                        })
                                    };
                                // FIXME: NAMING!!!
                                let left_peptide =
                                    last_peptide.ok_or(ConstructionError::NoPeptide)?;
                                let right_peptide = &right.peptide;
                                // FIXME: Awful...
                                if right_peptide.is_empty() {
                                    return Err(ConstructionError::NoPeptide);
                                }
                                let (donor, acceptor) = match descriptor {
                                    CrosslinkDescriptor::DonorAcceptor(l_idx, r_idx) => (
                                        lookup_donor(left_peptide, l_idx)?,
                                        lookup_acceptor(right_peptide, r_idx)?,
                                    ),
                                    CrosslinkDescriptor::AcceptorDonor(l_idx, r_idx) => (
                                        lookup_donor(right_peptide, r_idx)?,
                                        lookup_acceptor(left_peptide, l_idx)?,
                                    ),
                                };

                                let bond = if lateral_acceptor {
                                    LAT_CROSSLINK_BOND
                                } else {
                                    CROSSLINK_BOND
                                };
                                polymer.borrow_mut().bond_residues(bond, donor, acceptor)?;
                            }
                            Ok(())
                        };
                        match connection {
                            Connection::GlycosidicBond => gly_bond()?,
                            Connection::Crosslink(descriptors) => crosslink(descriptors)?,
                            Connection::Both(descriptors) => {
                                gly_bond()?;
                                crosslink(descriptors)?;
                            }
                        }
                    }

                    Ok((monomers, connections))
                },
            );
            parser(i)?
        };

        let polymer = polymer.into_inner();
        Ok((rest, Muropeptide {
            polymer,
            monomers,
            connections,
        }))
    }
}

/// Monomer = Glycan , [ "-" , Peptide ] | Peptide ;
// FIXME: Make private again
pub fn monomer<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<'s, Monomer> {
    let optional_peptide = opt(preceded(char('-'), cut(peptide(polymer))));
    let glycan_and_peptide = map_res(
        pair(glycan(polymer), optional_peptide),
        |(glycan, peptide)| {
            if let Some(peptide) = peptide {
                // SAFETY: Both the `glycan` and `peptide` parsers ensure at least one residue is present, so `.last()` and
                // `.first()` will never return `None`!
                let donor = *glycan.last().unwrap();
                let acceptor = peptide.first().unwrap().residue;

                polymer
                    .borrow_mut()
                    .bond_residues(STEM_BOND, donor, acceptor)?;
                Ok(Monomer { glycan, peptide })
            } else {
                Ok(Monomer {
                    glycan,
                    peptide: Vec::new(),
                })
            }
        },
    );

    let just_peptide = map(peptide(polymer), |peptide| Monomer {
        glycan: Vec::new(),
        peptide,
    });

    // FIXME: Add a `map_res` wrapping this final parser
    alt((glycan_and_peptide, just_peptide))
}

// =

/// Glycan = { Monosaccharide }- ;
fn glycan<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<'s, Vec<Monosaccharide>> {
    let parser = many1(monosaccharide(polymer));
    map_res(parser, |residues| {
        polymer
            .borrow_mut()
            .bond_chain(GLYCOSIDIC_BOND, &residues)?;
        Ok(residues)
    })
}

/// Peptide = { Amino Acid }- ;
fn peptide<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<'s, Vec<AminoAcid>> {
    let parser = many1(amino_acid(polymer));
    map_res(parser, |residues| {
        let residue_ids = residues.iter().map(|aa| aa.residue);
        polymer.borrow_mut().bond_chain(PEPTIDE_BOND, residue_ids)?;
        Ok(residues)
    })
}

// =

/// Monosaccharide = lowercase , [ Modifications ] ;
fn monosaccharide<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<'s, Monosaccharide> {
    let parser = pair(recognize(lowercase), opt(modifications(polymer)));
    map_res(parser, |(abbr, modifications)| {
        let residue = polymer.borrow_mut().new_residue(abbr)?;
        for modification in modifications.into_iter().flatten() {
            polymer
                .borrow_mut()
                .localize_modification(modification, residue)?;
        }

        Ok(residue)
    })
}

// FIXME: Damn... This is messy... Need to sort that out!
/// Amino Acid = Unbranched Amino Acid , [ Lateral Chain ] ;
fn amino_acid<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<'s, AminoAcid> {
    let parser = pair(unbranched_amino_acid(polymer), opt(lateral_chain(polymer)));
    map_res(parser, |(residue, lateral_chain)| {
        if let Some(LateralChain { direction, peptide }) = &lateral_chain {
            let c_to_n = || -> polychem::Result<_> {
                polymer
                    .borrow_mut()
                    .bond_residues(CTON_BOND, peptide[0], residue)?;
                // FIXME: Replace with `.clone()` then `.reverse()`?
                Ok(peptide.iter().copied().rev().collect())
            };
            let n_to_c = || -> polychem::Result<_> {
                polymer
                    .borrow_mut()
                    .bond_residues(NTOC_BOND, residue, peptide[0])?;
                // FIXME: This feels silly... Maybe better when `bond_chain` takes an `impl IntoIterator`?
                Ok(peptide.clone())
            };

            let chain: Vec<_> = match direction {
                // FIXME: This default should be moved to a configuration file and not be hard-coded!
                // FIXME: Furthermore, this should really return several possible polymers instead of picking one.
                // That might end up being difficult to do inline with this parsing stuff...
                // FIXME: If the peptide direction is initially unspecified, then I need to make sure that I update it
                // after working out which of the bonding patterns work (for pretty-printing later)!
                PeptideDirection::Unspecified => c_to_n().or_else(|_| n_to_c())?,
                PeptideDirection::CToN => c_to_n()?,
                PeptideDirection::NToC => n_to_c()?,
            };
            polymer.borrow_mut().bond_chain(PEPTIDE_BOND, &chain)?;
        }
        Ok(AminoAcid {
            residue,
            lateral_chain,
        })
    })
}

// =

/// Modifications = "(" , Any Modification ,
///   { { " " } , "," , { " " } , Any Modification } , ")" ;
fn modifications<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<'s, Vec<ModificationId>> {
    let separator = delimited(space0, char(','), space0);
    delimited(
        char('('),
        separated_list1(separator, any_modification(polymer)),
        char(')'),
    )
}

/// Unbranched Amino Acid = [ lowercase ] , uppercase , [ Modifications ] ;
fn unbranched_amino_acid<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<'s, UnbranchedAminoAcid> {
    let abbr = recognize(preceded(opt(lowercase), uppercase));
    let parser = pair(abbr, opt(modifications(polymer)));
    map_res(parser, |(abbr, modifications)| {
        let residue = polymer.borrow_mut().new_residue(abbr)?;
        for modification in modifications.into_iter().flatten() {
            polymer
                .borrow_mut()
                .localize_modification(modification, residue)?;
        }

        Ok(residue)
    })
}

// NOTE: These are not meant to be links, it's just EBNF
#[allow(clippy::doc_link_with_quotes)]
/// Lateral Chain = "[" , Peptide Direction , { Unbranched Amino Acid }- , "]" ;
fn lateral_chain<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<'s, LateralChain> {
    let peptide = many1(unbranched_amino_acid(polymer));
    let parser = delimited(char('['), peptide, char(']'));
    map(parser, |peptide| LateralChain {
        direction: PeptideDirection::Unspecified,
        peptide,
    })
}

// NOTE: These are not meant to be links, it's just EBNF
#[allow(clippy::doc_link_with_quotes)]
/// Peptide Direction = [ "<" (* C-to-N *) | ">" (* N-to-C *) ] ;
fn peptide_direction(i: &str) -> ParseResult<PeptideDirection> {
    map(opt(one_of("<>")), |c| match c {
        Some('<') => PeptideDirection::CToN,
        Some('>') => PeptideDirection::NToC,
        None => PeptideDirection::Unspecified,
        _ => unreachable!(),
    })(i)
}
// =

/// Identifier = letter , { letter | digit | "_" } ;
fn identifier(i: &str) -> ParseResult<&str> {
    // PERF: Could maybe avoid allocations by using `many0_count` instead, but needs benchmarking
    let parser = recognize(pair(alpha1, many0(alt((alphanumeric1, tag("_"))))));
    wrap_err(parser, MuropeptideErrorKind::ExpectedIdentifier)(i)
}

/// Any Modification = Named Modification | Offset Modification
pub fn any_modification<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<'s, ModificationId> {
    alt((named_modification(polymer), offset_modification(polymer)))
}

// FIXME: I probably need to add a lot of `wrap_err`s around these parsers!
/// Named Modification = [ Multiplier ] , Identifier
pub fn named_modification<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<'s, ModificationId> {
    let parser = pair(opt(multiplier), identifier);
    map_res(parser, |(multiplier, named_mod)| {
        polymer
            .borrow_mut()
            .new_modification(multiplier.unwrap_or_default(), named_mod)
            .map_err(Into::into)
    })
}

/// Offset Modification = Offset Kind , [ Multiplier ] ,
///   Chemical Composition ;
pub fn offset_modification<'c, 'a, 'p, 's>(
    polymer: &'c RefCell<Polymer<'a, 'p>>,
) -> impl FnMut(&'s str) -> ParseResult<'s, ModificationId> {
    let chemical_composition = chemical_composition(polymer.borrow().atomic_db());
    let parser = tuple((offset_kind, opt(multiplier), chemical_composition));

    map_res(parser, |(kind, multiplier, composition)| {
        polymer
            .borrow_mut()
            .new_offset_with_composition(kind, multiplier.unwrap_or_default(), composition)
            .map_err(Into::into)
    })
}

/// Connection
///   = "=" (* Crosslink *)
///   | "~" (* Glycosidic Bond *)
///   | ( "~=" | "=~" ) (* Both *)
///   ;
// FIXME: Not in love with carrying these empty Vecs around just to be populated later... Perhaps this should be a
// different type?
fn connection(i: &str) -> ParseResult<Connection> {
    // NOTE: The order here must start with the tags for `Both` — backtracking and trying the single-char parsers if
    // the two-char ones don't apply
    // PERF: There is probably a way to remove this backtracking and do things in a single-pass, but I don't know if
    // that would have a noticeable performance impact...
    let connection = alt((tag("="), tag("~")));
    let mut parser = map(connection, |c| match c {
        "~" => Connection::GlycosidicBond,
        "=" => Connection::Crosslink(Vec::new()),
        _ => unreachable!(),
    });
    // FIXME: Add error handling / reporting!
    parser(i)
}

/// Multiplier = Count , "x" ;
fn multiplier(i: &str) -> ParseResult<Count> {
    let mut parser = terminated(count, char('x'));
    // FIXME: Add error handling / reporting!
    parser(i)
}

/// Crosslinks = "(" , Crosslink Descriptors ,
///   { { " " } , "," , { " " } , Crosslink Descriptors } , ")" ;
fn crosslinks(i: &str) -> ParseResult<CrosslinkDescriptors> {
    let separator = delimited(space0, char(','), space0);
    let mut parser = delimited(
        char('('),
        separated_list1(separator, crosslink_descriptor),
        char(')'),
    );
    // FIXME: Add error handling / reporting!
    parser(i)
}

/// Crosslink Descriptors = Crosslink Descriptor ,
///   { { " " } , "&" , { " " } , Crosslink Descriptor } ;
fn crosslink_descriptors(i: &str) -> ParseResult<CrosslinkDescriptors> {
    let separator = delimited(space0, char('&'), space0);
    let mut parser = separated_list1(separator, crosslink_descriptor);
    // FIXME: Add error handling / reporting!
    parser(i)
}

/// Crosslink Descriptor = position ,
///   ( "-" (* Donor-Acceptor *)
///   | "=" (* Acceptor=Donor *)
///   ) , position ;
fn crosslink_descriptor(i: &str) -> ParseResult<CrosslinkDescriptor> {
    let mut parser = map(tuple((position, one_of("-="), position)), |(left, kind, right)| match kind {
        '-' => CrosslinkDescriptor::DonorAcceptor,
        '=' => CrosslinkDescriptor::AcceptorDonor,
        _ => unreachable!()
    }(left, right));
    // FIXME: Add error handling / reporting!
    parser(i)
}

/// position = "1" | "2" | "3" | "4" | "5" ;
fn position(i: &str) -> ParseResult<Position> {
    // SAFETY: If a 1, 2, 3, 4, or 5 is parsed, then calling `.to_digit()` will always return `Some(...)`, and casting
    // via `as u8` is also guaranteed not to overflow or truncate, since `5` is the largest parsable digit
    #[allow(clippy::cast_possible_truncation)]
    let mut parser = map(one_of("12345"), |p| p.to_digit(10).unwrap() as u8);
    // FIXME: Add error handling / reporting!
    parser(i)
}

type ParseResult<'a, O> = IResult<&'a str, O, LabeledParseError<'a, MuropeptideErrorKind>>;

#[derive(Clone, Eq, PartialEq, Debug, Diagnostic, Error)]
pub enum ConstructionError {
    #[diagnostic(transparent)]
    #[error(transparent)]
    PolychemError(#[from] Box<PolychemError>),

    // FIXME: Eventually make this error either more informative (printing the number of each), or replace it!
    #[error(
        "the number of described crosslinks (e.g. 3-3) doesn't match the number of crosslink connectors (=)"
    )]
    CrosslinkCountMismatch,

    // FIXME: Eventually make this error either more informative, or replace it!
    #[error("failed to find two glycan chains to bond")]
    NoGlycan,

    // FIXME: Eventually make this error either more informative, or replace it!
    #[error("failed to find two peptide stems to bond")]
    NoPeptide,

    // FIXME: Eventually make this error either more informative, or replace it!
    #[error("failed to look up the stem residues as described by the crosslink descriptor")]
    NoStemResidue,
}

#[derive(Clone, Eq, PartialEq, Debug, Diagnostic, Error)]
pub enum MuropeptideErrorKind {
    #[error(
        "expected an ASCII letter, optionally followed by any number of ASCII letters, digits, and underscores"
    )]
    ExpectedIdentifier,

    // FIXME: Kill this and merge into the error below!
    #[diagnostic(transparent)]
    #[error(transparent)]
    ConstructionError(ConstructionError),

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
impl<'a> FromExternalError<'a, ConstructionError> for MuropeptideErrorKind {
    const FATAL: bool = true;

    fn from_external_error(input: &'a str, e: ConstructionError) -> LabeledParseError<'a, Self> {
        LabeledParseError::new(input, Self::ConstructionError(e))
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
    use polychem::{AtomicDatabase, Charged, Massive, PolymerDatabase, Polymerizer};
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
            ($input:literal, $output:literal, $count:literal) => {
                let (rest, count) = multiplier($input).unwrap();
                assert_eq!((rest, u32::from(count)), ($output, $count));
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

    #[test]
    fn test_position() {
        // Valid Positions
        assert_eq!(position("1"), Ok(("", 1)));
        assert_eq!(position("2"), Ok(("", 2)));
        assert_eq!(position("3"), Ok(("", 3)));
        assert_eq!(position("4"), Ok(("", 4)));
        assert_eq!(position("5"), Ok(("", 5)));
        // Invalid Positions
        assert!(position("0").is_err());
        assert!(position("6").is_err());
        assert!(position("60").is_err());
        assert!(position("8422").is_err());
        assert!(position("01").is_err());
        assert!(position("00145").is_err());
        assert!(position("H").is_err());
        assert!(position("p").is_err());
        assert!(position("+H").is_err());
        assert!(position("[H]").is_err());
        // Multiple Positions
        assert_eq!(position("10"), Ok(("0", 1)));
        assert_eq!(position("1OH"), Ok(("OH", 1)));
        assert_eq!(position("42HeH"), Ok(("2HeH", 4)));
    }

    #[test]
    fn test_crosslink_descriptor() {
        // Valid Crosslink Descriptors
        assert_eq!(
            crosslink_descriptor("4-3"),
            Ok(("", CrosslinkDescriptor::DonorAcceptor(4, 3)))
        );
        assert_eq!(
            crosslink_descriptor("3=4"),
            Ok(("", CrosslinkDescriptor::AcceptorDonor(3, 4)))
        );
        assert_eq!(
            crosslink_descriptor("4-2"),
            Ok(("", CrosslinkDescriptor::DonorAcceptor(4, 2)))
        );
        assert_eq!(
            crosslink_descriptor("2=4"),
            Ok(("", CrosslinkDescriptor::AcceptorDonor(2, 4)))
        );
        assert_eq!(
            crosslink_descriptor("1-3"),
            Ok(("", CrosslinkDescriptor::DonorAcceptor(1, 3)))
        );
        assert_eq!(
            crosslink_descriptor("3=1"),
            Ok(("", CrosslinkDescriptor::AcceptorDonor(3, 1)))
        );
        // Invalid Crosslink Descriptors
        assert!(crosslink_descriptor("0").is_err());
        assert!(crosslink_descriptor("4-").is_err());
        assert!(crosslink_descriptor("3=").is_err());
        assert!(crosslink_descriptor("4->3").is_err());
        assert!(crosslink_descriptor("6").is_err());
        assert!(crosslink_descriptor("60").is_err());
        assert!(crosslink_descriptor("8422").is_err());
        assert!(crosslink_descriptor("01").is_err());
        assert!(crosslink_descriptor("00145").is_err());
        assert!(crosslink_descriptor("H").is_err());
        assert!(crosslink_descriptor("p").is_err());
        assert!(crosslink_descriptor("+H").is_err());
        assert!(crosslink_descriptor("[H]").is_err());
        // Multiple Crosslink Descriptors
        assert_eq!(
            crosslink_descriptor("4-3 & 3=3"),
            Ok((" & 3=3", CrosslinkDescriptor::DonorAcceptor(4, 3)))
        );
        assert_eq!(
            crosslink_descriptor("3=3) (Am)"),
            Ok((") (Am)", CrosslinkDescriptor::AcceptorDonor(3, 3)))
        );
    }

    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_crosslink_descriptors() {
        let da43 = CrosslinkDescriptor::DonorAcceptor(4, 3);
        let ad34 = CrosslinkDescriptor::AcceptorDonor(3, 4);
        let da33 = CrosslinkDescriptor::DonorAcceptor(3, 3);
        let ad33 = CrosslinkDescriptor::AcceptorDonor(3, 3);
        let da42 = CrosslinkDescriptor::DonorAcceptor(4, 2);
        let ad24 = CrosslinkDescriptor::AcceptorDonor(2, 4);
        let da13 = CrosslinkDescriptor::DonorAcceptor(1, 3);
        let ad31 = CrosslinkDescriptor::AcceptorDonor(3, 1);

        // Valid Crosslink Descriptors
        assert_eq!(crosslink_descriptors("4-3"), Ok(("", vec![da43])));
        assert_eq!(crosslink_descriptors("3=4"), Ok(("", vec![ad34])));
        assert_eq!(crosslink_descriptors("4-2"), Ok(("", vec![da42])));
        assert_eq!(crosslink_descriptors("2=4"), Ok(("", vec![ad24])));
        assert_eq!(crosslink_descriptors("1-3"), Ok(("", vec![da13])));
        assert_eq!(crosslink_descriptors("3=1"), Ok(("", vec![ad31])));
        assert_eq!(crosslink_descriptors("4-3&3=3"), Ok(("", vec![da43, ad33])));
        assert_eq!(
            crosslink_descriptors("4-3    &3=3"),
            Ok(("", vec![da43, ad33]))
        );
        assert_eq!(
            crosslink_descriptors("4-3&     3=3"),
            Ok(("", vec![da43, ad33]))
        );
        assert_eq!(
            crosslink_descriptors("4-3 & 3=3"),
            Ok(("", vec![da43, ad33]))
        );
        assert_eq!(
            crosslink_descriptors("4-3 & 3=3 & 3-3 & 2=4"),
            Ok(("", vec![da43, ad33, da33, ad24]))
        );
        // Invalid Crosslink Descriptors
        assert!(crosslink_descriptors("").is_err());
        assert!(crosslink_descriptors("0").is_err());
        assert!(crosslink_descriptors("4-").is_err());
        assert!(crosslink_descriptors("3=").is_err());
        assert!(crosslink_descriptors("4->3").is_err());
        assert!(crosslink_descriptors("& 3=3").is_err());
        assert!(crosslink_descriptors("6").is_err());
        assert!(crosslink_descriptors("60").is_err());
        assert!(crosslink_descriptors("8422").is_err());
        assert!(crosslink_descriptors("01").is_err());
        assert!(crosslink_descriptors("00145").is_err());
        assert!(crosslink_descriptors("H").is_err());
        assert!(crosslink_descriptors("p").is_err());
        assert!(crosslink_descriptors("+H").is_err());
        assert!(crosslink_descriptors("[H]").is_err());
        // Multiple Crosslink Descriptors
        assert_eq!(
            crosslink_descriptors("4-3 & 3=3 & 7-8"),
            Ok((" & 7-8", vec![da43, ad33]))
        );
        assert_eq!(
            crosslink_descriptors("3=3) (Am)"),
            Ok((") (Am)", vec![ad33]))
        );
    }

    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_connection() {
        let gly = || Connection::GlycosidicBond;
        let link = || Connection::Crosslink(Vec::new());

        // Valid Connections
        assert_eq!(connection("="), Ok(("", link())));
        assert_eq!(connection("~"), Ok(("", gly())));
        // Invalid Connections
        assert!(connection("").is_err());
        assert!(connection("gm").is_err());
        assert!(connection("gm-AEJA").is_err());
        assert!(connection("-AEJA").is_err());
        assert!(connection("()").is_err());
        assert!(connection("(3-)").is_err());
        assert!(connection("(2-4").is_err());
        assert!(connection("2-4 & 3=3").is_err());
        assert!(connection("4-3 & 3=3 & 7-8").is_err());
        assert!(connection("0").is_err());
        assert!(connection("4-").is_err());
        assert!(connection("3=").is_err());
        assert!(connection("4->3").is_err());
        assert!(connection("& 3=3").is_err());
        assert!(connection("6").is_err());
        assert!(connection("60").is_err());
        assert!(connection("8422").is_err());
        assert!(connection("01").is_err());
        assert!(connection("00145").is_err());
        assert!(connection("H").is_err());
        assert!(connection("p").is_err());
        assert!(connection("+H").is_err());
        assert!(connection("[H]").is_err());
        // Multiple Connections
        assert_eq!(connection("=="), Ok(("=", link())));
        assert_eq!(connection("~~"), Ok(("~", gly())));
    }

    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_named_modification() {
        let polymer = RefCell::new(POLYMERIZER.new_polymer());

        let mut err_named_modification = named_modification(&polymer);
        macro_rules! assert_named_modification {
            ($input:literal, $output:literal, $multiplier:literal, $name:expr, $mass:literal) => {
                let polymer = RefCell::new(POLYMERIZER.new_polymer());

                let (rest, parsed_id) = named_modification(&polymer)($input).unwrap();
                assert_eq!(rest, $output);

                let polymer = polymer.borrow();
                let modification = polymer
                    .modification(parsed_id)
                    .unwrap()
                    .clone()
                    .unwrap_unlocalized();

                let multiplier = modification.multiplier();
                assert_eq!(u32::from(multiplier), $multiplier);

                let name = modification.kind().clone().unwrap_named().name();
                assert_eq!(name, $name);

                assert_eq!(Decimal::from(polymer.monoisotopic_mass()), dec!($mass));
            };
        }

        // Valid Named Modifications
        assert_named_modification!("Am", "", 1, "Amidation", -0.98401558291);
        assert_named_modification!("Ac", "", 1, "O-Acetylation", 42.01056468403);
        assert_named_modification!("Poly", "", 1, "Wall Polymer Linkage", 77.95068082490);
        assert_named_modification!("DeAc", "", 1, "De-N-Acetylation", -42.01056468403);
        assert_named_modification!("Red", "", 1, "Reduced", 2.01565006446);
        assert_named_modification!("Anh", "", 1, "1,6-Anhydro", -18.01056468403);
        assert_named_modification!("1xAm", "", 1, "Amidation", -0.98401558291);
        assert_named_modification!("2xRed", "", 2, "Reduced", 4.03130012892);
        assert_named_modification!("3xAnh", "", 3, "1,6-Anhydro", -54.03169405209);
        // Invalid Named Modifications
        assert!(err_named_modification(" H2O").is_err());
        assert!(err_named_modification("1").is_err());
        assert!(err_named_modification("9999").is_err());
        assert!(err_named_modification("0").is_err());
        assert!(err_named_modification("00145").is_err());
        assert!(err_named_modification("+H").is_err());
        assert!(err_named_modification("[H]").is_err());
        assert!(err_named_modification("Øof").is_err());
        assert!(err_named_modification("-Ac").is_err());
        assert!(err_named_modification("_Ac").is_err());
        assert!(err_named_modification("+Am").is_err());
        assert!(err_named_modification("-2xAm").is_err());
        assert!(err_named_modification("(Am)").is_err());
        assert!(err_named_modification("-4xH2O").is_err());
        assert!(err_named_modification("-2p").is_err());
        assert!(err_named_modification("+C2H2O-2e").is_err());
        assert!(err_named_modification("-3xC2H2O-2e").is_err());
        assert!(err_named_modification("+NH3+p").is_err());
        assert!(err_named_modification("+2xD2O").is_err());
        assert!(err_named_modification("-2x[2H]2O").is_err());
        // Non-Existent Named Modifications
        assert!(err_named_modification("Blue").is_err());
        assert!(err_named_modification("Hydro").is_err());
        assert!(err_named_modification("1xAm2").is_err());
        assert!(err_named_modification("2xR_ed").is_err());
        // Multiple Named Modifications
        assert_named_modification!("Anh, Am", ", Am", 1, "1,6-Anhydro", -18.01056468403);
        assert_named_modification!("1xAm)JAA", ")JAA", 1, "Amidation", -0.98401558291);
    }

    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_offset_modification() {
        use polychem::OffsetKind::{Add, Remove};
        let polymer = RefCell::new(POLYMERIZER.new_polymer());

        let mut err_offset_modification = offset_modification(&polymer);
        macro_rules! assert_offset_modification {
            ($input:literal, $output:literal, $kind:ident, $multiplier:literal, $mass:literal, $charge:literal) => {
                let polymer = RefCell::new(POLYMERIZER.new_polymer());

                let (rest, parsed_id) = offset_modification(&polymer)($input).unwrap();
                assert_eq!(rest, $output);

                let polymer = polymer.borrow();
                let modification = polymer
                    .modification(parsed_id)
                    .unwrap()
                    .clone()
                    .unwrap_unlocalized();

                let kind = modification.kind().clone().unwrap_offset().kind();
                assert_eq!(kind, $kind);

                let multiplier = modification.multiplier();
                assert_eq!(u32::from(multiplier), $multiplier);

                assert_eq!(Decimal::from(polymer.monoisotopic_mass()), dec!($mass));
                assert_eq!(i64::from(polymer.charge()), $charge);
            };
        }

        // Valid Offset Modifications
        assert_offset_modification!("+H2O", "", Add, 1, 18.01056468403, 0);
        assert_offset_modification!("-H2O", "", Remove, 1, -18.01056468403, 0);
        assert_offset_modification!("+2xH2O", "", Add, 2, 36.02112936806, 0);
        assert_offset_modification!("-4xH2O", "", Remove, 4, -72.04225873612, 0);
        assert_offset_modification!("-2p", "", Remove, 1, -2.014552933242, -2);
        assert_offset_modification!("+H", "", Add, 1, 1.00782503223, 0);
        assert_offset_modification!("+C2H2O-2e", "", Add, 1, 42.009467524211870, 2);
        assert_offset_modification!("-3xC2H2O-2e", "", Remove, 3, -126.02840257263561, -6);
        assert_offset_modification!("+NH3+p", "", Add, 1, 18.033825567741, 1);
        assert_offset_modification!("+2xD2O", "", Add, 2, 40.04623635162, 0);
        assert_offset_modification!("-2x[2H]2O", "", Remove, 2, -40.04623635162, 0);
        assert_offset_modification!("+[37Cl]5-2p", "", Add, 1, 182.814960076758, -2);
        assert_offset_modification!("-NH2[99Tc]", "", Remove, 1, -114.92497486889, 0);
        // Invalid Offset Modifications
        assert!(err_offset_modification(" ").is_err());
        assert!(err_offset_modification("H2O").is_err());
        assert!(err_offset_modification("(-H2O)").is_err());
        assert!(err_offset_modification("+0xH2O").is_err());
        assert!(err_offset_modification("2xH2O").is_err());
        assert!(err_offset_modification("-2x3xH2O").is_err());
        assert!(err_offset_modification("-2x+H2O").is_err());
        assert!(err_offset_modification("+2[2H]").is_err());
        assert!(err_offset_modification("-[H+p]O").is_err());
        assert!(err_offset_modification("+NH2[100Tc]").is_err());
        // Multiple Offset Modifications
        assert_offset_modification!("+[37Cl]5-2p10", "10", Add, 1, 182.814960076758, -2);
        assert_offset_modification!("+[2H]2O*H2O", "*H2O", Add, 1, 20.02311817581, 0);
        assert_offset_modification!("+NH2{100Tc", "{100Tc", Add, 1, 16.01872406889, 0);
        assert_offset_modification!("+C11H12N2O2 H2O", " H2O", Add, 1, 204.08987763476, 0);
    }

    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_modifications() {
        let polymer = RefCell::new(POLYMERIZER.new_polymer());

        let mut err_modifications = modifications(&polymer);
        macro_rules! assert_modifications {
            ($input:literal, $output:literal, $count:literal, $mass:literal, $charge:literal) => {
                let polymer = RefCell::new(POLYMERIZER.new_polymer());

                let (rest, parsed_ids) = modifications(&polymer)($input).unwrap();
                assert_eq!(rest, $output);
                assert_eq!(parsed_ids.len(), $count);

                let polymer = polymer.borrow();
                assert_eq!(Decimal::from(polymer.monoisotopic_mass()), dec!($mass));
                assert_eq!(i64::from(polymer.charge()), $charge);
            };
        }
        // Valid Modifications
        assert_modifications!("(-H2O)", "", 1, -18.01056468403, 0);
        assert_modifications!("(+2xH2O)", "", 1, 36.02112936806, 0);
        assert_modifications!("(-2p)", "", 1, -2.014552933242, -2);
        assert_modifications!("(-3xC2H2O-2e)", "", 1, -126.02840257263561, -6);
        assert_modifications!("(+[37Cl]5-2p)", "", 1, 182.814960076758, -2);
        assert_modifications!("(Red)", "", 1, 2.01565006446, 0);
        assert_modifications!("(Anh)", "", 1, -18.01056468403, 0);
        assert_modifications!("(1xAm)", "", 1, -0.98401558291, 0);
        assert_modifications!("(2xRed)", "", 1, 4.03130012892, 0);
        assert_modifications!("(-OH, +NH2)", "", 2, -0.98401558291, 0);
        assert_modifications!("(Anh, +H2O)", "", 2, 0, 0);
        assert_modifications!("(Anh,+H2O)", "", 2, 0, 0);
        assert_modifications!("(Anh   ,+H2O)", "", 2, 0, 0);
        assert_modifications!("(Anh  ,  +H2O)", "", 2, 0, 0);
        assert_modifications!("(2xAnh, +3xH2O)", "", 2, 18.01056468403, 0);
        assert_modifications!("(Anh, Anh, +3xH2O)", "", 3, 18.01056468403, 0);
        assert_modifications!("(-H2, +Ca)", "", 2, 37.94694079854, 0);
        // NOTE: There is a super small mass defect (13.6 eV, or ~1e-8 u) stored in the binding energy between a proton
        // and electon — that's why this result is slightly different from the one above!
        assert_modifications!("(-2p, +Ca-2e)", "", 2, 37.946940769939870, 0);
        assert_modifications!("(+2p, -2p, +Ca-2e)", "", 3, 39.961493703181870, 2);
        // Invalid Modifications
        assert!(err_modifications(" ").is_err());
        assert!(err_modifications("H2O").is_err());
        assert!(err_modifications("(-H2O").is_err());
        assert!(err_modifications("(+0xH2O)").is_err());
        assert!(err_modifications("(2xH2O)").is_err());
        assert!(err_modifications("(-2x3xH2O)").is_err());
        assert!(err_modifications("(-2x+H2O)").is_err());
        assert!(err_modifications("(+2[2H])").is_err());
        assert!(err_modifications("(-[H+p]O)").is_err());
        assert!(err_modifications("(+NH2[100Tc])").is_err());
        assert!(err_modifications("( H2O)").is_err());
        assert!(err_modifications("(1)").is_err());
        assert!(err_modifications("(9999)").is_err());
        assert!(err_modifications("(0)").is_err());
        assert!(err_modifications("(00145)").is_err());
        assert!(err_modifications("([H])").is_err());
        assert!(err_modifications("(Øof)").is_err());
        assert!(err_modifications("(-Ac)").is_err());
        assert!(err_modifications("(_Ac)").is_err());
        assert!(err_modifications("(+Am)").is_err());
        assert!(err_modifications("(-2xAm)").is_err());
        assert!(err_modifications("((Am))").is_err());
        assert!(err_modifications("(Anh +H2O)").is_err());
        assert!(err_modifications("(Anh; +H2O)").is_err());
        // Non-Existent Modifications
        assert!(err_modifications("(Blue)").is_err());
        assert!(err_modifications("(Hydro)").is_err());
        assert!(err_modifications("(1xAm2)").is_err());
        assert!(err_modifications("(2xR_ed)").is_err());
        // Multiple Modifications
        assert_modifications!("(+[37Cl]5-2p)10", "10", 1, 182.814960076758, -2);
        assert_modifications!("(+[2H]2O)*H2O", "*H2O", 1, 20.02311817581, 0);
        assert_modifications!("(+NH2){100Tc", "{100Tc", 1, 16.01872406889, 0);
        assert_modifications!("(+C11H12N2O2) H2O", " H2O", 1, 204.08987763476, 0);
        assert_modifications!(
            "(2xAnh, +3xH2O)AA=gm-AEJA",
            "AA=gm-AEJA",
            2,
            18.01056468403,
            0
        );
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
                    .map(|id| polymer.residue(id.residue).unwrap().name())
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

    // FIXME: Add modification testing!
    #[test]
    #[allow(clippy::cognitive_complexity)]
    #[allow(clippy::too_many_lines)]
    fn test_monomer() {
        let polymer = RefCell::new(POLYMERIZER.new_polymer());

        let mut err_monomer = monomer(&polymer);
        macro_rules! assert_monomer_residues_and_masses {
            ($input:literal, $output:literal, $glycan:expr, $peptide:expr, $mono_mass:literal, $avg_mass:literal) => {
                let polymer = RefCell::new(POLYMERIZER.new_polymer());

                let (rest, Monomer { glycan, peptide }) = monomer(&polymer)($input).unwrap();
                assert_eq!(rest, $output);

                let polymer = polymer.borrow();
                let glycan: Vec<_> = glycan
                    .into_iter()
                    .map(|id| polymer.residue(id).unwrap().name())
                    .collect();
                let peptide: Vec<_> = peptide
                    .into_iter()
                    .map(|id| polymer.residue(id.residue).unwrap().name())
                    .collect();
                assert_eq!(glycan, Vec::<&str>::from($glycan));
                assert_eq!(peptide, Vec::<&str>::from($peptide));

                assert_eq!(Decimal::from(polymer.monoisotopic_mass()), dec!($mono_mass));
                assert_eq!(Decimal::from(polymer.average_mass()), dec!($avg_mass));
            };
        }

        // Valid Monomers
        assert_monomer_residues_and_masses!(
            "gmgm",
            "",
            [
                "N-Acetylglucosamine",
                "N-Acetylmuramic Acid",
                "N-Acetylglucosamine",
                "N-Acetylmuramic Acid"
            ],
            [],
            974.37031350523,
            974.91222678113779720
        );
        assert_monomer_residues_and_masses!(
            "gm",
            "",
            ["N-Acetylglucosamine", "N-Acetylmuramic Acid"],
            [],
            496.19043909463,
            496.46375660678381490
        );
        assert_monomer_residues_and_masses!(
            "g",
            "",
            ["N-Acetylglucosamine"],
            [],
            221.08993720530,
            221.20813124207411765
        );
        assert_monomer_residues_and_masses!(
            "m",
            "",
            ["N-Acetylmuramic Acid"],
            [],
            293.11106657336,
            293.27091179713952985
        );
        assert_monomer_residues_and_masses!(
            "AEJA",
            "",
            [],
            ["Alanine", "Glutamic Acid", "Diaminopimelic Acid", "Alanine"],
            461.21217759741,
            461.46756989305707095
        );
        assert_monomer_residues_and_masses!(
            "AyEJA",
            "",
            [],
            ["Alanine", "γ-Glutamate", "Diaminopimelic Acid", "Alanine"],
            461.21217759741,
            461.46756989305707095
        );
        assert_monomer_residues_and_masses!(
            "AE",
            "",
            [],
            ["Alanine", "Glutamic Acid"],
            218.09027155793,
            218.20748877514586040
        );
        assert_monomer_residues_and_masses!(
            "A",
            "",
            [],
            ["Alanine"],
            89.04767846918,
            89.09330602867854225
        );
        assert_monomer_residues_and_masses!(
            "gm-AEJA",
            "",
            ["N-Acetylglucosamine", "N-Acetylmuramic Acid"],
            ["Alanine", "Glutamic Acid", "Diaminopimelic Acid", "Alanine"],
            939.39205200801,
            939.91604006741105325
        );
        assert_monomer_residues_and_masses!(
            "gm-AyEJA",
            "",
            ["N-Acetylglucosamine", "N-Acetylmuramic Acid"],
            ["Alanine", "γ-Glutamate", "Diaminopimelic Acid", "Alanine"],
            939.39205200801,
            939.91604006741105325
        );
        // Invalid Monomers
        assert!(err_monomer("-AEJA").is_err());
        assert!(err_monomer("[GGGGG]").is_err());
        assert!(err_monomer("(Am)").is_err());
        // Non-Existent Monomer Residues & Bonds
        assert!(err_monomer("y").is_err());
        assert!(err_monomer("fp").is_err());
        assert!(err_monomer("AEJiA").is_err());
        assert!(err_monomer("AQyK").is_err());
        assert!(err_monomer("g-A").is_err());
        // Multiple Monomers
        assert_monomer_residues_and_masses!(
            "gm,AEJ",
            ",AEJ",
            ["N-Acetylglucosamine", "N-Acetylmuramic Acid"],
            [],
            496.19043909463,
            496.46375660678381490
        );
        assert_monomer_residues_and_masses!(
            "xAJgmK",
            "AJgmK",
            ["Unknown Monosaccharide"],
            [],
            0.0,
            0.0
        );
        assert_monomer_residues_and_masses!(
            "AE=gm-AEJ",
            "=gm-AEJ",
            [],
            ["Alanine", "Glutamic Acid"],
            218.09027155793,
            218.20748877514586040
        );
        assert_monomer_residues_and_masses!(
            "AeeK",
            "eeK",
            [],
            ["Alanine"],
            89.04767846918,
            89.09330602867854225
        );
        assert_monomer_residues_and_masses!(
            "gm-AE=gm-AEJA",
            "=gm-AEJA",
            ["N-Acetylglucosamine", "N-Acetylmuramic Acid"],
            ["Alanine", "Glutamic Acid"],
            696.27014596853,
            696.65595894949984270
        );
    }

    // FIXME: Add a test that checks all of the errors using `assert_miette_snapshot`! Maybe make that a crate?
}
