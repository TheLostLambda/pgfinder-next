// FIXME: This should probably be made private at some point!
mod muropeptide;
mod polychem;

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
