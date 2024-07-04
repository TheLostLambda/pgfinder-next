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

#[derive(Clone, Debug, Default)]
struct Residue<'p> {
    terminals: Vec<Terminal<'p>>,
    bonds: Vec<Bond<'p>>,
}

#[derive(Copy, Clone, Debug)]
struct Bond<'p> {
    // FIXME: The naming `end` doesn't fit well with the `Terminal` here...
    end: Terminal<'p>,
    target: usize,
}

impl<'p> Bond<'p> {
    const fn donating_to(abbr: BondAbbr<'p>, acceptor: usize) -> Self {
        let end = Terminal::Donor(abbr);
        Self {
            end,
            target: acceptor,
        }
    }

    const fn accepting_from(abbr: BondAbbr<'p>, donor: usize) -> Self {
        let end = Terminal::Acceptor(abbr);
        Self { end, target: donor }
    }
}

// PERF: Not sold on this "tombstone" approach with `Option` — might be better to re-index the sub-graphs so their
// vectors can be shrunk?
#[derive(Clone, Debug)]
struct Fragment<'p>(Vec<Option<Residue<'p>>>);

// FIXME: This was written way too quickly... Take a look back at this...
// FIXME: Also maybe should die, as it's only useful for debugging?
impl Display for Fragment<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        writeln!(f, "digraph {{")?;
        for (id, residue) in self.0.iter().enumerate() {
            let Some(Residue { terminals, bonds }) = residue else {
                continue;
            };

            let terminals = terminals.iter().map(ToString::to_string).join(", ");
            writeln!(f, r#"  {id} [label="{id}\n[{terminals}]"]"#)?;

            // NOTE: The `.filter()` here prevents the double-printing of edges
            for Bond { end, target } in bonds.iter().filter(|b| id < b.target) {
                match end {
                    Terminal::Donor(abbr) => {
                        writeln!(f, r#"  {id} -> {target} [label="{abbr}"]"#)?;
                    }
                    Terminal::Acceptor(abbr) => {
                        writeln!(f, r#"  {target} -> {id} [label="{abbr}"]"#)?;
                    }
                }
            }
        }
        writeln!(f, "}}")?;
        Ok(())
    }
}

#[derive(Clone, Debug)]
struct NodeMapping(Vec<ResidueId>);

impl NodeMapping {
    fn new(polymer: &Polymer) -> Self {
        // FIXME: Needs a second look?
        // PERF: Is there some way to collect and sort at the same time? Would that even be faster?
        let mut node_mapping: Vec<_> = polymer.residue_ids().collect();
        node_mapping.sort_unstable();
        Self(node_mapping)
    }

    // NOTE: Keeping `id` as a ref, since it needs to be one for `binary_search()` and because it comes from
    // `bond_refs()` to begin with!
    #[allow(clippy::trivially_copy_pass_by_ref)]
    fn index(&self, id: &ResidueId) -> usize {
        // SAFETY: Panics if the `id` isn't found
        self.0.binary_search(id).unwrap()
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
        let polymer = self.polymer();
        let node_mapping = NodeMapping::new(polymer);
        let fragment = Fragment::new(&node_mapping, polymer);
        eprintln!("{fragment}");
    }
}

impl<'p> Fragment<'p> {
    // FIXME: This currently assumes that the `polymer` provided consists of a single connected component. If this is
    // not the case, I should report an error to the user!
    fn new(node_mapping: &NodeMapping, polymer: &Polymer<'_, 'p>) -> Self {
        let mut residues = vec![Some(Residue::default()); node_mapping.0.len()];

        for BondInfo(ResidueGroup(donor, _), bond, ResidueGroup(acceptor, _)) in polymer.bond_refs()
        {
            let abbr = bond.abbr();
            let donor = node_mapping.index(donor);
            let acceptor = node_mapping.index(acceptor);

            let mut push_bond =
                // FIXME: Shocking failure of type inference here with that `usize`...
                |residue: usize, bond| residues[residue].as_mut().unwrap().bonds.push(bond);

            push_bond(donor, Bond::donating_to(abbr, acceptor));
            push_bond(acceptor, Bond::accepting_from(abbr, donor));
        }

        Self(residues)
    }
}

// impl<'p> Index<usize> for Fragment<'p> {
//     type Output = Residue<'p>;

//     fn index(&self, index: usize) -> &Self::Output {
//         self.0[index].as_ref().unwrap()
//     }
// }

// SEE NOTES FROM APRIL 8TH!
// use DashMap or quick-cache for a global fragment cache

// FIXME: Tests need adding!
