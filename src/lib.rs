// FIXME: This should probably be made private at some point!
pub mod parser;

use std::{collections::HashSet, hash::Hash, iter::repeat};

use memoize::memoize;
use parser::{LateralChain, Modifications};
use petgraph::{stable_graph::NodeIndex, Graph};
use phf::phf_map;
use polars::prelude::*;
use rust_decimal::Decimal;
use rust_decimal_macros::dec;
use serde::{Deserialize, Serialize};

// OPEN QUESTIONS =============================================================
// 1) Which direction do lateral chains run off from mDAP? (from the amine!)
// 2) What should I do when I have several B-ion (O+) termini?

// Findings ===================================================================
// 1) B-ions were off by an electron mass — they are slighly lighter now!

// Things To Do ===============================================================
// 1) Add support for ambiguous modifications: `gm-AEJA (-Ac)`
// 2) Add a converter for generating Windows-legal filenames!
// 3) Add in rules for checking legal crosslinks and modifications!
// 4) Allow lateral chains to be nested (for M. luteus, ffs...)
// 5) Add lactyl ions for residues attached to the MurNAc
// 6) Check for immonium ions!

// FIXME: These masses need to be checked against the mass_calc databases Steph has vetted!
// FIXME: These moieties need to store a list of side-chain functional groups (for checking
// which crosslinks are legal and which modifications are allowed where!)
// FIXME: Replace these hard-coded maps with CSV(?) files!
// Masses computed using https://mstools.epfl.ch/info/
// Checked against http://www.matrixscience.com/help/aa_help.html
// Used to update https://github.com/Mesnage-Org/rhizobium-pg-pipeline/blob/7f3a322624c027f5c42b796c6a1c0a1d7d81dbb0/Data/Constants/masses_table.csv
pub static RESIDUES: phf::Map<char, Moiety> = phf_map! {
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
    'O' => Moiety::new("O", "Ornithine", dec!(114.079313)),
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
    'X' => Moiety::new("X", "Unknown Amino Acid", dec!(0.000000)),
    'x' => Moiety::new("x", "Unknown Monosaccharide", dec!(0.000000)),
};
// Used to update https://github.com/Mesnage-Org/rhizobium-pg-pipeline/blob/7f3a322624c027f5c42b796c6a1c0a1d7d81dbb0/Data/Constants/mods_table.csv
pub static MODIFICATIONS: phf::Map<&str, Moiety> = phf_map! {
    "+" => Moiety::new("+", "Proton", dec!(1.007276)),
    "-" => Moiety::new("-", "Electron", dec!(0.000549)),
    "H" => Moiety::new("H", "Hydrogen", dec!(1.007825)),
    "OH" => Moiety::new("OH", "Hydroxy", dec!(17.002740)),
    "Ac" => Moiety::new("Ac", "Acetyl", dec!(42.010565)),
    // "Anh" => Moiety::new("Anh", "Anhydro", dec!()),
    // FIXME: Should this be -H2O, -Anh? +Anh? Mutally exclusive with the +2H
    // from MurNAc reduction! -20.026215
    // FIXME: Need to deal with modifications that can only apply to certain
    // residues! Most of them involved adding *and* removing groups! -0.984016
    // Lac: OCH3COCH -> 72.021129 (Needed to find A3)
};
// FIXME: Consider merging this with the rest of the RESIDUES table!
pub static BRANCH_RESIDUES: phf::Map<&str, SidechainGroup> = phf_map! {
    "E" => SidechainGroup::Carboxyl,
    "J" => SidechainGroup::Amine, // FIXME: Is actually `Both`!
    "K" => SidechainGroup::Amine,
};

// FIXME: Add a parsed data validation step! Use Crepe to enforce some rules!
// FIXME: Cache or memoise calculations of things like mass
#[derive(Clone, Debug)]
pub struct Peptidoglycan {
    // FIXME: Replace the `pub`s with accessor methods?
    pub name: String,
    pub graph: Graph<Residue, ()>,
}

impl Peptidoglycan {
    // FIXME: Replace `String` with a proper error type!
    // FIXME: Sopping code... Needs some DRYing!
    pub fn new(structure: &str) -> Result<Self, String> {
        // FIXME: Handle this error properly!
        let (_, (monomers, crosslinks)) = parser::multimer(structure).unwrap();
        let mut residue_id = 0;
        let mut graph = Graph::new();
        for (glycan, peptide) in monomers {
            // Build the glycan chain
            let mut last_monosaccharide = graph.add_node(Residue::new(&mut residue_id, &glycan[0]));

            // Add H to the "N-terminal"
            graph[last_monosaccharide]
                .modifications
                .push(Modification::Add(MODIFICATIONS["H"]));

            // Add edges
            for monosaccharide in &glycan[1..] {
                let monosaccharide = graph.add_node(Residue::new(&mut residue_id, monosaccharide));
                graph.add_edge(last_monosaccharide, monosaccharide, ());
                last_monosaccharide = monosaccharide;
            }

            // Add H2 to the reduced end of the glycan chain
            // FIXME: Can GlcNAc be reduced too? Or only MurNAc?
            // FIXME: This doesn't work for Anhydro!
            graph[last_monosaccharide]
                .modifications
                .extend(repeat(Modification::Add(MODIFICATIONS["H"])).take(2));

            if let Some(peptide) = peptide {
                let (abbr, modifications, lateral_chain) = &peptide[0];
                let mut last_amino_acid =
                    graph.add_node(Residue::new(&mut residue_id, (abbr, modifications)));

                build_lateral_chain(&mut graph, &mut residue_id, last_amino_acid, lateral_chain);

                // Join the glycan chain and peptide stem
                // FIXME: Needs to validate that this is a MurNAc
                graph.add_edge(last_monosaccharide, last_amino_acid, ());

                for (abbr, modifications, lateral_chain) in &peptide[1..] {
                    let amino_acid =
                        graph.add_node(Residue::new(&mut residue_id, (abbr, modifications)));
                    graph.add_edge(last_amino_acid, amino_acid, ());
                    last_amino_acid = amino_acid;
                    build_lateral_chain(
                        &mut graph,
                        &mut residue_id,
                        last_amino_acid,
                        lateral_chain,
                    );
                }

                // Add OH to the "C-terminal"
                graph[last_amino_acid]
                    .modifications
                    .push(Modification::Add(MODIFICATIONS["OH"]));
            } else {
                // Add OH to the "C-terminal"
                graph[last_monosaccharide]
                    .modifications
                    .push(Modification::Add(MODIFICATIONS["OH"]));
            }
        }

        // FIXME: Still very wet, needs breaking into more subfunctions!
        fn build_lateral_chain(
            graph: &mut Graph<Residue, ()>,
            residue_id: &mut usize,
            last_amino_acid: NodeIndex,
            lateral_chain: &Option<LateralChain>,
        ) {
            let sidechain_group = BRANCH_RESIDUES.get(graph[last_amino_acid].moiety.abbr);
            if let (Some(sidechain_group), Some(lateral_chain)) = (sidechain_group, lateral_chain) {
                let mut last_lateral_amino_acid =
                    graph.add_node(Residue::new(residue_id, &lateral_chain[0]));

                // FIXME: Dripping wet...
                match sidechain_group {
                    SidechainGroup::Amine => {
                        graph.add_edge(last_lateral_amino_acid, last_amino_acid, ());
                        // Remove H from the branch-point residue
                        graph[last_amino_acid]
                            .modifications
                            .push(Modification::Remove(MODIFICATIONS["H"]));

                        for lateral_amino_acid in &lateral_chain[1..] {
                            let lateral_amino_acid =
                                graph.add_node(Residue::new(residue_id, lateral_amino_acid));
                            graph.add_edge(lateral_amino_acid, last_lateral_amino_acid, ());
                            last_lateral_amino_acid = lateral_amino_acid;
                        }

                        // Add H to the "N-terminal"
                        graph[last_lateral_amino_acid]
                            .modifications
                            .push(Modification::Add(MODIFICATIONS["H"]));
                    }
                    SidechainGroup::Carboxyl => {
                        graph.add_edge(last_amino_acid, last_lateral_amino_acid, ());
                        // Remove OH from the branch-point residue
                        graph[last_amino_acid]
                            .modifications
                            .push(Modification::Remove(MODIFICATIONS["OH"]));

                        for lateral_amino_acid in &lateral_chain[1..] {
                            let lateral_amino_acid =
                                graph.add_node(Residue::new(residue_id, lateral_amino_acid));
                            graph.add_edge(last_lateral_amino_acid, lateral_amino_acid, ());
                            last_lateral_amino_acid = lateral_amino_acid;
                        }

                        // Add OH to the "C-terminal"
                        graph[last_lateral_amino_acid]
                            .modifications
                            .push(Modification::Add(MODIFICATIONS["OH"]));
                    }
                }
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

// FIXME: Make the cache shared between threads (and speed up fragment hash functions!)
#[memoize]
pub fn fragment(f: Fragment) -> HashSet<Fragment> {
    let Fragment { structure, termini } = f;
    if structure.raw_nodes().is_empty() {
        return HashSet::new();
    }
    let mut fragments = HashSet::new();
    for edge in structure.raw_edges() {
        let from = edge.source();
        let to = edge.target();
        let (b, y) = expand_from(&structure, from, to);
        let mut b_graph = structure.clone();
        b_graph.retain_nodes(|_, n| b.contains(&n));
        let mut b_termini = termini.clone();
        b_termini.retain(|(k, _)| {
            b.iter()
                .map(|&i| structure[i].id)
                .collect::<Vec<_>>()
                .contains(k)
        });
        b_termini.insert((structure[from].id, TerminalIon::B));
        fragments.insert(Fragment {
            structure: b_graph,
            termini: b_termini,
        });
        let mut y_graph = structure.clone();
        y_graph.retain_nodes(|_, n| y.contains(&n));
        let mut termini = termini.clone();
        termini.retain(|(k, _)| {
            y.iter()
                .map(|&i| structure[i].id)
                .collect::<Vec<_>>()
                .contains(k)
        });
        termini.insert((structure[to].id, TerminalIon::Y));
        fragments.insert(Fragment {
            structure: y_graph,
            termini,
        });
    }
    let internal: Vec<_> = fragments.iter().flat_map(|f| fragment(f.clone())).collect();
    fragments.extend(internal);
    fragments
}

fn expand_from(
    graph: &Graph<Residue, ()>,
    from: NodeIndex,
    to: NodeIndex,
) -> (Vec<NodeIndex>, Vec<NodeIndex>) {
    let mut visited = vec![to];
    let mut unexplored = vec![from];
    while let Some(node) = unexplored.pop() {
        visited.push(node);
        unexplored.extend(
            graph
                .neighbors_undirected(node)
                .filter(|n| !visited.contains(n)),
        );
    }
    visited.remove(0);
    let left = visited;
    // Just swap from and to?
    let mut visited = vec![from];
    let mut unexplored = vec![to];
    while let Some(node) = unexplored.pop() {
        visited.push(node);
        unexplored.extend(
            graph
                .neighbors_undirected(node)
                .filter(|n| !visited.contains(n)),
        );
    }
    visited.remove(0);
    (left, visited)
}

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub enum Modification {
    Add(Moiety),
    Remove(Moiety),
}

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub enum SidechainGroup {
    Amine,
    Carboxyl,
    // Both, // FIXME: Out of action until I know how mDAP lateral chains work!
}

// TODO: Add support for more ion types!
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub enum TerminalIon {
    B,
    Y,
}

#[derive(Clone, Debug)]
pub struct Fragment {
    structure: Graph<Residue, ()>,
    termini: HashSet<(usize, TerminalIon)>,
}

impl Hash for Fragment {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let mut nodes = self.structure.raw_nodes().to_vec();
        nodes.sort_by_key(|n| n.weight.id);
        for node in nodes {
            node.weight.hash(state);
        }
    }
}

impl PartialEq for Fragment {
    fn eq(&self, other: &Self) -> bool {
        let mut self_nodes = self
            .structure
            .raw_nodes()
            .iter()
            .map(|n| n.weight.clone())
            .collect::<Vec<_>>();
        let mut other_nodes = other
            .structure
            .raw_nodes()
            .iter()
            .map(|n| n.weight.clone())
            .collect::<Vec<_>>();
        self_nodes.sort_by_key(|n| n.id);
        other_nodes.sort_by_key(|n| n.id);
        self_nodes == other_nodes
    }
}

impl Eq for Fragment {}

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Ion {
    structure: String,
    mass: Decimal,
    charge: usize, // FIXME: Is this the right type?
}

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Residue {
    id: usize,
    moiety: Moiety,
    modifications: Vec<Modification>,
}

impl Residue {
    pub fn new(id: &mut usize, residue: impl Into<Residue>) -> Self {
        *id += 1;
        Self {
            id: *id,
            ..residue.into()
        }
    }
    pub fn monoisotopic_mass(&self) -> Decimal {
        self.modifications
            .iter()
            .fold(self.moiety.mass, |mass, modification| match modification {
                Modification::Add(m) => mass + m.mass,
                Modification::Remove(m) => mass - m.mass,
            })
    }
}

impl From<&(char, Option<Modifications>)> for Residue {
    fn from((abbr, modifications): &(char, Option<Modifications>)) -> Self {
        Residue::from((abbr, modifications))
    }
}

impl From<(&char, &Option<Modifications>)> for Residue {
    fn from((abbr, modifications): (&char, &Option<Modifications>)) -> Self {
        let moiety = RESIDUES[abbr];
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
            id: 0,
            moiety,
            modifications,
        }
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Moiety {
    pub abbr: &'static str,
    pub name: &'static str,
    pub mass: Decimal,
}

impl Moiety {
    pub const fn new(abbr: &'static str, name: &'static str, mass: Decimal) -> Self {
        Self { abbr, name, mass }
    }
}

// FIXME: Should this be a `From` impl?
pub fn fragments_to_df(fragments: &HashSet<Fragment>) -> DataFrame {
    let mut ion_types = Vec::new();
    let mut structures = Vec::new();
    let mut ion_masses = Vec::new();
    for Fragment { structure, termini } in fragments {
        let mut residues: Vec<_> = structure
            .raw_nodes()
            .iter()
            .cloned()
            .map(|n| n.weight)
            .collect();
        residues.sort_by_key(|n| n.id);
        let mut name = String::new();
        let mut mass = Decimal::default();
        for residue in residues {
            // Build up structure name
            if !name.is_empty() {
                name.push('-');
            }
            name.push_str(residue.moiety.abbr);
            name.push_str(&residue.id.to_string());
            // Build up fragment mass
            mass += residue.monoisotopic_mass();
        }
        let mut termini: Vec<_> = termini.iter().map(|(_, t)| t).collect();
        // FIXME: Would all of the sorting be faster with sort_unstable?
        termini.sort();
        if termini
            .iter()
            .all(|&t| [TerminalIon::B, TerminalIon::Y].contains(t))
        {
            // FIXME: Should termini be HashMap<TerminalIon, Vec<usize>>? So the length of the
            // vector is the number of termini of that type? And the val is a list of ids?
            let b = termini.iter().filter(|&&&t| t == TerminalIon::B).count();
            let y = termini.iter().filter(|&&&t| t == TerminalIon::Y).count();
            let ion_type = match b {
                0 => {
                    mass += MODIFICATIONS["+"].mass;
                    "Y".repeat(y)
                }
                1 => {
                    mass -= MODIFICATIONS["-"].mass;
                    ["B", &"Y".repeat(y)].concat()
                }
                n => {
                    // "???".to_string()
                    panic!(
                        "Sorry mate, no clue how to handle multple ({n}) B termini... {termini:?}"
                    )
                }
            };
            mass += MODIFICATIONS["H"].mass * Decimal::from(y);
            ion_types.push(ion_type);
            let ion_mass = mass.round_dp(4).to_string();
            ion_masses.push(ion_mass);
            structures.push(name);
        } else {
            panic!("Wat? {termini:?}");
        }
    }
    let terminal_count: Vec<_> = ion_types.iter().map(|s| s.len() as i32).collect();
    let float_masses: Vec<f64> = ion_masses.iter().map(|s| s.parse().unwrap()).collect();
    // FIXME: Should actually think about handling errors here...
    let mut df =
        df!("Termini" => terminal_count, "Float Mass" => float_masses, "Ion Type" => ion_types, "Ion Mass" => ion_masses, "Structure" => structures).unwrap();
    // FIXME: Is there any performance gained by doing this in-place?
    df.sort_in_place(["Termini", "Ion Type", "Float Mass", "Structure"], false);
    df.drop_in_place("Termini");
    df.drop_in_place("Float Mass");
    df
}

impl From<Peptidoglycan> for Fragment {
    fn from(value: Peptidoglycan) -> Self {
        Self {
            structure: value.graph,
            termini: Default::default(),
        }
    }
}

// FIXME: Keep writing messy code, then write tests, then refactor
#[cfg(test)]
mod tests {
    use insta::{assert_debug_snapshot, assert_snapshot};
    use std::error::Error;

    use petgraph::dot::Dot;

    use super::*;

    #[test]
    fn linear_monomer_graph() -> Result<(), Box<dyn Error>> {
        let pg = Peptidoglycan::new("gm-AEJA")?;
        assert_eq!(pg.monoisotopic_mass().round_dp(4), dec!(941.4077));
        assert_debug_snapshot!(Dot::new(&pg.graph));
        Ok(())
    }

    #[test]
    fn linear_monomer_fragments() -> Result<(), Box<dyn Error>> {
        let fragments = fragment(Peptidoglycan::new("gm-AEJA")?.into());
        assert_eq!(fragments.len(), 20);

        let mut csv = Vec::new();
        CsvWriter::new(&mut csv).finish(&mut fragments_to_df(&fragments))?;
        assert_snapshot!(String::from_utf8(csv)?);
        Ok(())
    }
}
