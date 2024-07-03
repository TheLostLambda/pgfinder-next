use std::fmt::{self, Display, Formatter};

use itertools::Itertools;
use polychem::{BondInfo, Polymer, ResidueGroup, ResidueId};

// FIXME: Consider using newtype? Especially if this is made public!
type BondAbbr<'p> = &'p str;

#[derive(Copy, Clone, Debug)]
enum Terminal<'p> {
    Donor(BondAbbr<'p>),
    Acceptor(BondAbbr<'p>),
}

impl Display for Terminal<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Terminal::Donor(abbr) => write!(f, "{abbr}-D"),
            Terminal::Acceptor(abbr) => write!(f, "{abbr}-A"),
        }
    }
}

#[derive(Clone, Debug)]
struct Residue<'p> {
    id: ResidueId,
    terminals: Vec<Terminal<'p>>,
}

impl<'p> Residue<'p> {
    const fn new(id: ResidueId) -> Self {
        let terminals = Vec::new();
        Self { id, terminals }
    }
}

#[derive(Copy, Clone, Debug)]
struct Bond<'p> {
    donor: ResidueId,
    abbr: BondAbbr<'p>,
    acceptor: ResidueId,
}

impl<'p> Bond<'p> {
    // FIXME: Should this be a From impl instead / as well?
    const fn new(bond_info: &BondInfo<'_, 'p>) -> Self {
        let &BondInfo(ResidueGroup(donor, _), ref bond, ResidueGroup(acceptor, _)) = bond_info;
        let abbr = bond.abbr();
        Self {
            donor,
            abbr,
            acceptor,
        }
    }
}

#[derive(Clone, Debug)]
struct Fragment<'p> {
    depth: u32,
    residues: Vec<Residue<'p>>,
    bonds: Vec<Bond<'p>>,
}

// FIXME: This was written way too quickly... Take a look back at this...
impl Display for Fragment<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        writeln!(f, "digraph {{")?;
        for Residue { id, terminals } in &self.residues {
            let terminals = terminals.iter().map(ToString::to_string).join(", ");
            writeln!(f, r#"  {id} [label="{id}\n[{terminals}]"]"#)?;
        }
        for Bond {
            donor,
            abbr,
            acceptor,
        } in &self.bonds
        {
            writeln!(f, r#"  {donor} -> {acceptor} [label="{abbr}"]"#)?;
        }
        writeln!(f, "}}")?;
        Ok(())
    }
}

// DESIGN: This `Sized` bound could be moved to the method-level with `where Self: Sized`, which would allow this trait
// to be implemented for unsized types, just without those methods, but I don't think this trait makes much sense
// without a `.fragment()` method (which currently requires a `Self: Sized`), so I've elected to make the whole trait
// `Sized` for now
pub trait Dissociable: Sized {
    #[must_use]
    fn polymer(&self) -> &Polymer;
    #[must_use]
    fn new_fragment(&self, fragmented_polymer: Polymer) -> Self;
    #[must_use]
    fn fragment(&self) -> Vec<Self> {
        todo!()
    }
    // FIXME: Only for testing! Remove!
    fn dbg_fragment(&self) {
        eprintln!("{}", Fragment::new(self.polymer()));
    }
}

impl<'p> Fragment<'p> {
    fn new<'r>(polymer: &'r Polymer<'_, 'p>) -> Self {
        let depth = 0;
        let residues: Vec<Residue<'p>> = polymer.residue_ids().map(Residue::new).collect();
        let bonds: Vec<Bond<'p>> = polymer.bond_refs().map(Bond::new).collect();
        Self {
            depth,
            residues,
            bonds,
        }
    }
}
// SEE NOTES FROM APRIL 8TH!
// use scc::HashCache, or HashIndex, or TreeIndex
// or use DashMap?

// FIXME: Tests need adding!
