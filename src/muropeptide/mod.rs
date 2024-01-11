//! Responsible for parsing strings into meaningful `Muropeptide` structures

pub mod parser;

use crate::polymers::{Modification, Residue};

pub struct Muropeptide {
    monomers: Vec<Monomer>,
    connections: Vec<Connection>,
    modifications: Vec<Modification>,
}

struct Monomer {
    glycan: Vec<Monosaccharide>,
    peptide: Vec<AminoAcid>,
}

type Connection = Vec<ConnectionKind>;

type Monosaccharide = Residue;

struct AminoAcid {
    residue: Residue,
    lateral_chain: Option<LateralChain>,
}

enum ConnectionKind {
    GlycosidicBond,
    Crosslink(Vec<CrosslinkDescriptor>),
}

struct LateralChain {
    direction: PeptideDirection,
    peptide: Vec<UnbranchedAminoAcid>,
}

enum CrosslinkDescriptor {
    DonorAcceptor(Position, Position),
    AcceptorDonor(Position, Position),
}

enum PeptideDirection {
    Unspecified,
    CToN,
    NToC,
}

type UnbranchedAminoAcid = Residue;

type Position = u8;
