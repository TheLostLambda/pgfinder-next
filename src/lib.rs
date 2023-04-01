// FIXME: This should probably be made private at some point!
pub mod parser;

use std::iter::repeat;

use petgraph::Graph;
use phf::phf_map;
use rust_decimal::Decimal;
use rust_decimal_macros::dec;
use serde::{Deserialize, Serialize};

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
    'm' => Moiety::new("m", "MurNAc", dec!(275.100502)),
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

// FIXME: Cache or memoise calculations of things like mass
#[derive(Clone, Debug)]
pub struct Peptidoglycan {
    name: String,
    graph: Graph<Residue, ()>,
}

impl Peptidoglycan {
    // FIXME: Replace `String` with a proper error type!
    // FIXME: Sopping code... Needs some DRYing!
    pub fn new(structure: &str) -> Result<Self, String> {
        // FIXME: Handle this error properly!
        let (_, (monomers, crosslinks)) = parser::multimer(structure).unwrap();
        // FIXME: Get rid of this explicit type
        let mut graph: Graph<Residue, ()> = Graph::new();
        for (glycan, peptide) in monomers {
            let glycan: Vec<_> = glycan
                .into_iter()
                .map(|monosaccharide| graph.add_node(monosaccharide.into()))
                .collect();
            // FIXME: This will crash with single nodes
            for w in glycan.windows(2) {
                graph.add_edge(w[0], w[1], ());
            }
            // Add H to the "N-terminal"
            graph[glycan[0]]
                .modifications
                .push(Modification::Add(MODIFICATIONS["H"]));

            // Add H2 to the reduced end of the glycan chain
            // FIXME: Can GlcNAc be reduced too? Or only MurNAc?
            graph[glycan[glycan.len() - 1]]
                .modifications
                .extend(repeat(Modification::Add(MODIFICATIONS["H"])).take(2));

            if let Some(peptide) = peptide {
                let peptide: Vec<_> = peptide
                    .into_iter()
                    // FIXME: Actually do something with the lateral chain!
                    .map(|(abbr, modifcations, _lateral_chain)| {
                        graph.add_node((abbr, modifcations).into())
                    })
                    .collect();

                // FIXME: Needs to validate that this is a MurNAc
                graph.add_edge(glycan[glycan.len()-1], peptide[0], ());
                // FIXME: This will crash with single nodes
                for w in peptide.windows(2) {
                    graph.add_edge(w[0], w[1], ());
                }
                // Add OH to the "C-terminal"
                graph[peptide[peptide.len() - 1]]
                    .modifications
                    .push(Modification::Add(MODIFICATIONS["OH"]));
            } else {
                // Add OH to the "C-terminal"
                graph[glycan[glycan.len() - 1]]
                    .modifications
                    .push(Modification::Add(MODIFICATIONS["OH"]));
            }
        }
        Ok(Self {
            name: structure.to_string(),
            graph,
        })
    }

    pub fn monoisotopic_mass(&self) -> Decimal {
        self.graph
            .raw_nodes()
            .iter()
            .map(|residue| residue.weight.monoisotopic_mass())
            .sum()
    }
}

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub enum Modification {
    Add(Moiety),
    Remove(Moiety),
}

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Residue {
    moiety: Moiety,
    modifications: Vec<Modification>,
}

impl Residue {
    pub fn monoisotopic_mass(&self) -> Decimal {
        self.modifications
            .iter()
            .fold(self.moiety.mass, |mass, modification| match modification {
                Modification::Add(m) => mass + m.mass,
                Modification::Remove(m) => mass - m.mass,
            })
    }
}

impl From<(char, Option<parser::Modifications>)> for Residue {
    fn from((abbr, modifications): (char, Option<parser::Modifications>)) -> Self {
        let moiety = RESIDUES[&abbr];
        let modifications = modifications
            .iter()
            .flatten()
            // FIXME: This Modification conversion could be improved using `strum`?
            .map(|modification| match modification {
                parser::Modification::Add(abbr) => Modification::Add(MODIFICATIONS[abbr]),
                parser::Modification::Remove(abbr) => Modification::Remove(MODIFICATIONS[abbr]),
            })
            .collect();
        Self {
            moiety,
            modifications,
        }
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Moiety {
    abbr: &'static str,
    name: &'static str,
    mass: Decimal,
}

impl Moiety {
    pub const fn new(abbr: &'static str, name: &'static str, mass: Decimal) -> Self {
        Self { abbr, name, mass }
    }
}

#[cfg(test)]
mod tests {
    use std::error::Error;

    use petgraph::dot::{Config, Dot};

    use super::*;

    #[test]
    fn it_works() -> Result<(), Box<dyn Error>> {
        let pg = dbg!(Peptidoglycan::new("gm-AEJYW")?);
        println!("{:?}", Dot::with_config(&pg.graph, &[Config::EdgeNoLabel]));
        dbg!(dbg!(pg.monoisotopic_mass()).round_dp(4));
        panic!();
        Ok(())
    }
}
