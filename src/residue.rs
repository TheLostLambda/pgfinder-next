/// Responsible for parsing strings into meaningful `Peptidoglycan` structures
mod peptidoglycan {
    use super::chemistry::{Modification, Residue};

    pub struct Muropeptide {
        monomers: Vec<Monomer>,
        connections: Vec<Connection>,
        modifications: Vec<Modification>,
    }

    // NOTE: Use `Option` fields when you expect them to usually be filled
    struct Monomer {
        glycan: Option<Vec<Monosaccharide>>,
        peptide: Option<Vec<AminoAcid>>,
    }

    type Monosaccharide = Residue;

    // NOTE: Convert to `enum` variants when `Option` fields would be regularly
    // unfilled, but keep as `struct` if the `enum` variants share data that's
    // useful regardless of the variant
    struct AminoAcid {
        residue: Residue,
        lateral_chain: Option<LateralChain>,
    }

    struct LateralChain {
        direction: PeptideDirection,
        peptide: Vec<UnbranchedAminoAcid>,
    }

    enum PeptideDirection {
        Unspecified,
        CToN,
        NToC,
    }

    type UnbranchedAminoAcid = Residue;

    type Connection = Vec<ConnectionKind>;

    enum ConnectionKind {
        GlycosidicBond,
        Crosslink(Vec<CrosslinkDescriptor>),
    }

    enum CrosslinkDescriptor {
        DonorAcceptor(u8, u8),
        AcceptorDonor(u8, u8),
    }
}

/// The building blocks of all molecular graphs
mod chemistry {
    use rust_decimal::Decimal;

    pub struct Residue {
        id: usize,
        moiety: Moiety,
        modifications: Vec<Modification>,
    }

    struct Moiety {
        abbr: String,
        name: String,
        mass: Decimal,
    }

    // TODO: This needs an `impl` for getting its mass
    pub struct Modification {
        multiplier: u32,
        kind: ModificationKind,
    }

    enum ModificationKind {
        Predefined(Moiety),
        ChemicalOffset(ChemicalFormula),
    }

    type ChemicalFormula = Vec<(Element, u32)>;

    // QUESTION: Is it ever important to store average mass? Or just monoisotopic?
    type Element = Moiety;
}
