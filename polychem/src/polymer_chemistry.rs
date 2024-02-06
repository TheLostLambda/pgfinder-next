use crate::ChemicalComposition;

// FIXME: Check that field names here line up with those in `lib.rs`!
struct PolymerChemistryKdl {
    residues: ResiduesKdl,
}

struct ResiduesKdl {
    types: Vec<ResidueTypeKdl>,
    residues: Vec<ResidueKdl>,
}

struct ResidueTypeKdl {
    name: String,
    functional_groups: Vec<FunctionalGroupKdl>,
}

struct ResidueKdl {
    residue_type: String,
    abbr: String,
    name: String,
    composition: ChemicalCompositionKdl,
    functional_groups: Vec<FunctionalGroupKdl>,
}

struct FunctionalGroupKdl {
    name: String,
    location: String,
}

#[derive(Debug, Default)]
struct ChemicalCompositionKdl(ChemicalComposition);
