//! Responsible for parsing strings into meaningful `Muropeptide` structures
mod parser;
mod smiles_database;

use std::fmt::{self, Display, Formatter};

use itertools::Itertools;
use miette::Diagnostic;
use nom_miette::{LabeledError, final_parser};
use parser::{MuropeptideErrorKind, muropeptide};
use std::sync::LazyLock;
// FIXME: Blocks need separating and reordering!
use polychem::{
    AtomicDatabase, AverageMass, BondId, Charged, GroupState, Massive, ModificationInfo,
    MonoisotopicMass, Polymer, PolymerDatabase, Polymerizer, ResidueId, errors::PolychemError,
};
use smiles_database::SmilesDatabase;
use smithereens::Dissociable;
use thiserror::Error;

// FIXME: Need to think about if these should really live in another KDL config?
const PEPTIDE_BOND: &str = "Pep";
const GLYCOSIDIC_BOND: &str = "Gly";
const STEM_BOND: &str = "Stem";
const NTOC_BOND: &str = "NToC";
const CTON_BOND: &str = "CToN";
const CROSSLINK_BOND: &str = "Link";
const LAT_CROSSLINK_BOND: &str = "Lat-Link";

// FIXME: These need more thought / are a temporary hack!
static ATOMIC_DB: LazyLock<AtomicDatabase> = LazyLock::new(AtomicDatabase::default);
pub static POLYMER_DB: LazyLock<PolymerDatabase> = LazyLock::new(|| {
    PolymerDatabase::new(
        &ATOMIC_DB,
        "polymer_database.kdl",
        include_str!("../data/polymer_database.kdl"),
    )
    .unwrap()
});

// FIXME: This maybe shouldn't be here long term? Needs some thought...
pub static POLYMERIZER: LazyLock<Polymerizer> =
    LazyLock::new(|| Polymerizer::new(&ATOMIC_DB, &POLYMER_DB));

// FIXME: This maybe shouldn't be here long term? Needs some thought...
pub static SMILES_DB: LazyLock<SmilesDatabase> = LazyLock::new(|| {
    SmilesDatabase::new(
        "smiles_database.kdl",
        include_str!("../data/smiles_database.kdl"),
    )
    .unwrap()
});

const AUTO_MODS: [&str; 1] = ["Red"];

#[derive(Debug)]
pub struct Muropeptide<'a, 'p> {
    polymer: Polymer<'a, 'p>,
    monomers: Vec<Monomer>,
    connections: Vec<Connection>,
}

#[derive(Debug, Clone)]
struct Monomer {
    glycan: Vec<Monosaccharide>,
    peptide: Vec<AminoAcid>,
}

type Monosaccharide = ResidueId;

#[derive(Debug, Clone)]
struct AminoAcid {
    residue: UnbranchedAminoAcid,
    lateral_chain: Option<LateralChain>,
}

// FIXME: This should store the `BondId` of each connection, so we know when it's removed!
#[derive(Clone, Debug, Eq, PartialEq)]
enum Connection {
    GlycosidicBond,
    Crosslink(CrosslinkDescriptors),
    Both(CrosslinkDescriptors),
}

#[derive(Debug, Clone)]
struct LateralChain {
    direction: PeptideDirection,
    peptide: Vec<UnbranchedAminoAcid>,
}

type CrosslinkDescriptors = Vec<CrosslinkDescriptor>;

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum CrosslinkDescriptor {
    DonorAcceptor(Position, Position),
    AcceptorDonor(Position, Position),
}

#[derive(Copy, Clone, Debug)]
enum PeptideDirection {
    Unspecified,
    CToN,
    NToC,
}

type UnbranchedAminoAcid = ResidueId;

type Position = u8;

// ===

#[derive(Clone, Debug, Error, Diagnostic)]
pub enum Error {
    #[diagnostic(transparent)]
    #[error(transparent)]
    ParseError(#[from] LabeledError<MuropeptideErrorKind>),

    #[diagnostic(transparent)]
    #[error(transparent)]
    PolychemError(#[from] Box<PolychemError>),
}
// FIXME: Maybe `Box` this error?
pub type Result<T> = std::result::Result<T, Error>;

impl<'a, 'p> Muropeptide<'a, 'p> {
    // FIXME: Messy stuff here! Needs a look! Especially the error type!
    pub fn new(polymerizer: &Polymerizer<'a, 'p>, structure: impl AsRef<str>) -> Result<Self> {
        let mut muropeptide = final_parser(muropeptide(polymerizer))(structure.as_ref())?;

        // FIXME: This is a temporary, hard-coded version of this feature. It should eventually be moved to a
        // configuration file where the list of automatic modifications can be set by the user!
        for abbr in AUTO_MODS {
            muropeptide.polymer.modify_polymer(abbr)?;
        }

        Ok(muropeptide)
    }

    // FIXME: I hate this... Get rid of it.
    #[must_use]
    pub fn oligomerization_state(&self) -> usize {
        self.monomers.len()
    }
}

impl Massive for Muropeptide<'_, '_> {
    fn monoisotopic_mass(&self) -> MonoisotopicMass {
        self.polymer.monoisotopic_mass()
    }

    fn average_mass(&self) -> AverageMass {
        self.polymer.average_mass()
    }
}

impl Charged for Muropeptide<'_, '_> {
    fn charge(&self) -> polychem::Charge {
        self.polymer.charge()
    }
}

impl<'s, 'a: 's, 'p: 's> Dissociable<'s, 'a, 'p> for Muropeptide<'a, 'p> {
    fn polymer(&self) -> &Polymer<'a, 'p> {
        &self.polymer
    }

    // FIXME: Fucking hideous.
    fn new_fragment(
        &self,
        fragmented_polymer: Polymer<'a, 'p>,
        lost_residues: Vec<ResidueId>,
        _broken_bonds: Vec<BondId>,
    ) -> Self {
        // FIXME: Obviously incomplete!
        let monomers = self
            .monomers
            .iter()
            .map(|Monomer { glycan, peptide }| {
                // FIXME: Likely inefficient, with linear search and not HashSets? Also this should be a function...
                let glycan: Vec<_> = glycan
                    .iter()
                    .copied()
                    .filter(|id| !lost_residues.contains(id))
                    .collect();
                let mut lateral_peptides = Vec::new();
                let peptide: Vec<_> = peptide
                    .iter()
                    .filter_map(|aa| {
                        let residue = aa.residue;
                        if !lost_residues.contains(&residue) {
                            let lateral_chain = aa.lateral_chain.as_ref().and_then(|chain| {
                                let peptide: Vec<_> = chain
                                    .peptide
                                    .iter()
                                    .copied()
                                    .filter(|id| !lost_residues.contains(id))
                                    .collect();
                                if peptide.is_empty() {
                                    None
                                } else {
                                    Some(LateralChain {
                                        direction: chain.direction,
                                        peptide,
                                    })
                                }
                            });
                            Some(AminoAcid {
                                residue,
                                lateral_chain,
                            })
                        } else if let Some(LateralChain { peptide, .. }) = &aa.lateral_chain {
                            let peptide: Vec<_> = peptide
                                .iter()
                                .copied()
                                .filter(|id| !lost_residues.contains(id))
                                .collect();
                            if !peptide.is_empty() {
                                lateral_peptides.push(peptide);
                            }
                            None
                        } else {
                            None
                        }
                    })
                    .collect();
                if peptide.is_empty() && glycan.is_empty() && !lateral_peptides.is_empty() {
                    // FIXME: Probably not true eventually...
                    assert_eq!(lateral_peptides.len(), 1);
                    let peptide = lateral_peptides[0]
                        .iter()
                        .map(|&residue| AminoAcid {
                            residue,
                            lateral_chain: None,
                        })
                        .collect();
                    Monomer { glycan, peptide }
                } else {
                    Monomer { glycan, peptide }
                }
            })
            .collect();
        let connections = self.connections.clone();
        Self {
            polymer: fragmented_polymer,
            monomers,
            connections,
        }
    }
}

// FIXME: Oh god.
impl Display for Muropeptide<'_, '_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let connections = self.connections.clone();
        let mut cross_links = Vec::new();

        for (i, left_monomer) in self.monomers.iter().enumerate() {
            if left_monomer.glycan.is_empty() && left_monomer.peptide.is_empty() {
                continue;
            }
            display_monomer(f, &self.polymer, left_monomer)?;
            if let Some(right_monomer) = self.monomers.get(i + 1) {
                if right_monomer.glycan.is_empty() && right_monomer.peptide.is_empty() {
                    continue;
                }
                if let Some(connection) = connections.get(i) {
                    match connection {
                        Connection::GlycosidicBond => write!(f, "~")?,
                        Connection::Crosslink(descriptors) => {
                            write!(f, "=")?;
                            assert_eq!(descriptors.len(), 1);
                            cross_links.push(descriptors[0].to_string());
                        }
                        Connection::Both(descriptors) => {
                            write!(f, "~=")?;
                            assert_eq!(descriptors.len(), 1);
                            cross_links.push(descriptors[0].to_string());
                        }
                    }
                }
            }
        }

        let cross_links = cross_links.join(", ");
        if !cross_links.is_empty() {
            write!(f, " ({cross_links})")?;
        }

        Ok(())
    }
}

impl Display for CrosslinkDescriptor {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::DonorAcceptor(l, r) => write!(f, "{l}-{r}"),
            Self::AcceptorDonor(l, r) => write!(f, "{l}={r}"),
        }
    }
}

// FIXME: Change all of these individual functions into a trait, then implement it for all of the sub-components. The
// trait could be called something like `DisplayMoiety` and could be a bit like the `ValidateInto` trait?
fn display_monomer(f: &mut Formatter, polymer: &Polymer, monomer: &Monomer) -> fmt::Result {
    let Monomer { glycan, peptide } = monomer;
    for &monosaccharide in glycan {
        display_residue(f, polymer, monosaccharide)?;
    }
    if !glycan.is_empty() && !peptide.is_empty() {
        write!(f, "-")?;
    }
    for amino_acid in peptide {
        display_residue(f, polymer, amino_acid.residue)?;
        if let Some(chain) = &amino_acid.lateral_chain {
            write!(f, "[")?;
            for &residue in &chain.peptide {
                display_residue(f, polymer, residue)?;
            }
            write!(f, "]")?;
        }
    }
    // TODO: Bring this back once I have more ion types, but this just noise for now (when it's all protons)
    // let unlocalised_modifications = polymer
    //     .modification_refs()
    //     .filter_map(|info| {
    //         if let ModificationInfo::Unlocalized(modification) = info {
    //             Some(modification.to_string())
    //         } else {
    //             None
    //         }
    //     })
    //     .join(", ");
    // if !unlocalised_modifications.is_empty() {
    //     write!(f, " ({unlocalised_modifications})")?;
    // }
    Ok(())
}

// FIXME: Sickening.
fn display_residue(f: &mut Formatter, polymer: &Polymer, residue: ResidueId) -> fmt::Result {
    // FIXME: What to do about the unwrap()? Is that fine here?
    let residue = polymer.residue(residue).unwrap();
    let abbr = residue.abbr();
    let named_mods = residue.functional_groups().filter_map(|(_, gs)| {
        if let &GroupState::Modified(id) = gs {
            let ModificationInfo::Named(named_mod, _) = polymer.modification(id).unwrap() else {
                unreachable!();
            };
            Some(named_mod.abbr().to_owned())
        } else {
            None
        }
    });
    let offset_mods = residue.offset_modifications().map(|id| {
        let ModificationInfo::Offset(modification, _) = polymer.modification(id).unwrap() else {
            unreachable!();
        };
        modification.to_string()
    });
    let modifications = named_mods
        .chain(offset_mods)
        // FIXME: Awful hacks to format things in a way Steph likes...
        .map(|m| {
            m.replace("Red", "r")
                .replace("-H2O", "")
                .replace("-2xH2O", "")
        })
        .filter(|m| !m.is_empty())
        .join(", ");

    if modifications.is_empty() {
        write!(f, "{abbr}")
    } else {
        write!(f, "{abbr}({modifications})")
    }
}

// OPEN QUESTIONS =============================================================
// 1) Which direction do lateral chains run off from mDAP? (from the amine!)
// 2) What should I do when I have several B-ion (O+) termini?

// Findings ===================================================================
// 1) B-ions were off by an electron mass — they are slighly lighter now!

// Things To Do ===============================================================
// 1) Add support for ambiguous modifications: `gm-AEJA (-Ac)` — maybe (-2Ac)
// 2) Add a converter for generating Windows-legal filenames!
// 3) Add in rules for checking legal crosslinks and modifications!
// 4) Allow lateral chains to be nested (for M. luteus, ffs...)
// 5) Add lactyl ions for residues attached to the MurNAc
// 6) Check for immonium ions!
// 7) Generate ions with charges greater than 1 (up to 3)! Output m/z values!

// Other Thoughts / Notes =====================================================
// I'm going to need to either revisit the use of `~` for glycosidic bonds, or
// (and this is probably preferable) add a converter to output escaped `~`
// characters that are compatable with Excel!
//
// I should add the lactyl-group to my list of sugar moieties as `l`. This will
// let me build graphs for and search for lactyl peptides. I can probably
// insert a lactyl-group linker automatically whenever I'm attaching a glycan
// to a peptide. That way the user doesn't need to include an implied `l` every
// time, but the graph can still be fragmented properly and the outputted
// structures can explicitly include the lactyl group if the rest of the glycan
// is lost from the peptide.
//
// Better yet, I could have the `-` represent the lactyl? So things like `gm-`
// or `-AEJA` could indicate where the lactyl ends up?

// Some Rules? ================================================================
// Could fill this tool with rules like that to apply, saying how each
// returned structure was validated (MS2 or a biochemical rule) so that you get
// a result of confirmed structures and a set of explicitly listed assumptions
// for each search!
//
// 1) Structures like `gm-AEJH-gm-AEJ` must be 3-3 crosslinks because the `H`
// in position four has been exchanged (need the concept of a parent stem) and
// that means there was no `gm-AEJHX` parent to form a 4-3 crosslink with! This
// could also rule out structures like `gm-AEJDP`!

// FIXME: These masses need to be checked against the mass_calc databases Steph has vetted!
// FIXME: These moieties need to store a list of side-chain functional groups (for checking
// which crosslinks are legal and which modifications are allowed where!)
// FIXME: Replace these hard-coded maps with CSV(?) files!
// Masses computed using https://mstools.epfl.ch/info/
// Checked against http://www.matrixscience.com/help/aa_help.html
// Used to update https://github.com/Mesnage-Org/rhizobium-pg-pipeline/blob/7f3a322624c027f5c42b796c6a1c0a1d7d81dbb0/Data/Constants/masses_table.csv
// pub static RESIDUES: phf::Map<char, Moiety> = phf_map! {
//     'A' => Moiety::new("A", "Alanine", dec!(71.037114)),
//     'C' => Moiety::new("C", "Cysteine", dec!(103.009185)),
//     'D' => Moiety::new("D", "Aspartic Acid", dec!(115.026943)),
//     'E' => Moiety::new("E", "Glutamic Acid", dec!(129.042593)),
//     'F' => Moiety::new("F", "Phenylalanine", dec!(147.068414)),
//     'G' => Moiety::new("G", "Glycine", dec!(57.021464)),
//     'H' => Moiety::new("H", "Histidine", dec!(137.058912)),
//     'I' => Moiety::new("I", "Isoleucine", dec!(113.084064)),
//     'J' => Moiety::new("J", "Diaminopimelic Acid", dec!(172.084792)),
//     'K' => Moiety::new("K", "Lysine", dec!(128.094963)),
//     'L' => Moiety::new("L", "Leucine", dec!(113.084064)),
//     'M' => Moiety::new("M", "Methionine", dec!(131.040485)),
//     'N' => Moiety::new("N", "Asparagine", dec!(114.042927)),
//     'O' => Moiety::new("O", "Ornithine", dec!(114.079313)),
//     'P' => Moiety::new("P", "Proline", dec!(97.052764)),
//     'Q' => Moiety::new("Q", "Glutamine", dec!(128.058578)),
//     'R' => Moiety::new("R", "Arginine", dec!(156.101111)),
//     'S' => Moiety::new("S", "Serine", dec!(87.032028)),
//     'T' => Moiety::new("T", "Threonine", dec!(101.047678)),
//     'U' => Moiety::new("U", "Selenocysteine", dec!(150.953636)),
//     'V' => Moiety::new("V", "Valine", dec!(99.068414)),
//     'W' => Moiety::new("W", "Tryptophan", dec!(186.079313)),
//     'Y' => Moiety::new("Y", "Tyrosine", dec!(163.063329)),
//     'g' => Moiety::new("g", "GlcNAc", dec!(203.079373)),
//     'm' => Moiety::new("m", "MurNAc", dec!(275.100502)),
//     // These are placeholder residues, allowing modifications to entirely determine the residue mass
//     'X' => Moiety::new("X", "Unknown Amino Acid", dec!(0.000000)),
//     'x' => Moiety::new("x", "Unknown Monosaccharide", dec!(0.000000)),
// };
// // Used to update https://github.com/Mesnage-Org/rhizobium-pg-pipeline/blob/7f3a322624c027f5c42b796c6a1c0a1d7d81dbb0/Data/Constants/mods_table.csv
// pub static MODIFICATIONS: phf::Map<&str, Moiety> = phf_map! {
//     "+" => Moiety::new("+", "Proton", dec!(1.007276)),
//     "-" => Moiety::new("-", "Electron", dec!(0.000549)),
//     "H" => Moiety::new("H", "Hydrogen", dec!(1.007825)),
//     "OH" => Moiety::new("OH", "Hydroxy", dec!(17.002740)),
//     "Ac" => Moiety::new("Ac", "Acetyl", dec!(42.010565)),
//     // "Anh" => Moiety::new("Anh", "Anhydro", dec!()),
//     // FIXME: Should this be -H2O, -Anh? +Anh? Mutally exclusive with the +2H
//     // from MurNAc reduction! -20.026215. I should extend the language to not
//     // require a + or - before every modification!
//     // FIXME: Need to deal with modifications that can only apply to certain
//     // residues! Most of them involved adding *and* removing groups! -0.984016
//     // Lac: OCH3COCH -> 72.021129 (Needed to find A3)
// };
// // FIXME: Consider merging this with the rest of the RESIDUES table!
// pub static BRANCH_RESIDUES: phf::Map<&str, SidechainGroup> = phf_map! {
//     "E" => SidechainGroup::Carboxyl,
//     "J" => SidechainGroup::Amine, // FIXME: Is actually `Both`!
//     "K" => SidechainGroup::Amine,
// };
