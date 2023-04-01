// FIXME: This should probably be made private at some point!
pub mod parsers;

// FIXME: Need a new module for types that are universal! Like `Modification`!
use parsers::peptidoglycan::{Modification, Modifications};
use petgraph::Graph;
use phf::phf_map;
use rust_decimal::Decimal;
use rust_decimal_macros::dec;

// FIXME: These masses need to be checked against the mass_calc databases Steph has vetted!
// Masses computed using https://mstools.epfl.ch/info/
// Checked against http://www.matrixscience.com/help/aa_help.html
// Used to update https://github.com/Mesnage-Org/rhizobium-pg-pipeline/blob/7f3a322624c027f5c42b796c6a1c0a1d7d81dbb0/Data/Constants/masses_table.csv
static RESIDUES: phf::Map<char, Moiety> = phf_map! {
    'A' => Moiety::new("A", "Alanine", dec!(71.037114)),
    'C' => Moiety::new("C", "Cysteine", dec!(103.009185)),
    'D' => Moiety::new("D", "Aspartic Acid", dec!(115.026943)),
    'E' => Moiety::new("E", "Glutamic Acid", dec!(129.042593)),
    'F' => Moiety::new("F", "Phenylalanine", dec!(147.068414)),
    'G' => Moiety::new("G", "Glycine", dec!(57.021464)),
    'H' => Moiety::new("H", "Histidine", dec!(137.058912)),
    'I' => Moiety::new("I", "Isoleucine", dec!(113.084064)),
    'J' => Moiety::new("J", "Diaminopimelic Acid", dec!(172.084792)),
    'K' => Moiety::new("K", "Lysine", dec!(128.094963)),
    'L' => Moiety::new("L", "Leucine", dec!(113.084064)),
    'M' => Moiety::new("M", "Methionine", dec!(131.040485)),
    'N' => Moiety::new("N", "Asparagine", dec!(114.042927)),
    'P' => Moiety::new("P", "Proline", dec!(97.052764)),
    'Q' => Moiety::new("Q", "Glutamine", dec!(128.058578)),
    'R' => Moiety::new("R", "Arginine", dec!(156.101111)),
    'S' => Moiety::new("S", "Serine", dec!(87.032028)),
    'T' => Moiety::new("T", "Threonine", dec!(101.047678)),
    'U' => Moiety::new("U", "Selenocysteine", dec!(150.953636)),
    'V' => Moiety::new("V", "Valine", dec!(99.068414)),
    'W' => Moiety::new("W", "Tryptophan", dec!(186.079313)),
    'Y' => Moiety::new("Y", "Tyrosine", dec!(163.063329)),
    'g' => Moiety::new("g", "GlcNAc", dec!(203.079373)),
    'm' => Moiety::new("m", "MurNAc Alditol", dec!(277.116152)),
    // These are placeholder residues, allowing modifications to entirely determine the residue mass
    'X' => Moiety::new("X", "Unknown Amino Acid", dec!(0.0)),
    'x' => Moiety::new("x", "Unknown Monosaccharide", dec!(0.0)),
};
// Used to update https://github.com/Mesnage-Org/rhizobium-pg-pipeline/blob/7f3a322624c027f5c42b796c6a1c0a1d7d81dbb0/Data/Constants/mods_table.csv
static MODIFICATIONS: phf::Map<&str, Moiety> = phf_map! {
    "+" => Moiety::new("+", "Proton", dec!(1.007276)),
    "H" => Moiety::new("H", "Hydrogen", dec!(1.007825)),
    "OH" => Moiety::new("OH", "Hydroxy", dec!(17.002740)),
};

// pub struct Molecule {
//     name: String,
//     graph: Graph,
// }

pub struct Residue {
    id: usize,
    moiety: Moiety,
    mods: Modifications,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Moiety {
    abbr: &'static str,
    name: &'static str,
    mass: Decimal,
}

impl Moiety {
    // FIXME: Replace this with a real implementation
    pub const fn new(abbr: &'static str, name: &'static str, mass: Decimal) -> Self {
        Self { abbr, name, mass }
    }
}

#[cfg(test)]
mod tests {
    use std::error::Error;

    use super::*;

    #[test]
    fn it_works() -> Result<(), Box<dyn Error>> {
        Ok(())
    }
}
