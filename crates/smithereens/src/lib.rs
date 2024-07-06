use std::{
    fmt::{self, Display, Formatter},
    ops::{Index, IndexMut},
};

use ahash::{HashSet, HashSetExt};
use itertools::Itertools;
use polychem::{BondInfo, Polymer, ResidueGroup, ResidueId};

// FIXME: Consider using newtype? Especially if this is made public!
type BondAbbr<'p> = &'p str;

#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
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

#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
struct Residue<'p> {
    terminals: Vec<Terminal<'p>>,
    bonds: Vec<Bond<'p>>,
}

#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
struct Bond<'p> {
    // FIXME: The naming `end` doesn't fit well with the `Terminal` here...
    end: Terminal<'p>,
    target: NodeId,
}

impl<'p> Bond<'p> {
    const fn donating_to(abbr: BondAbbr<'p>, acceptor: NodeId) -> Self {
        let end = Terminal::Donor(abbr);
        Self {
            end,
            target: acceptor,
        }
    }

    const fn accepting_from(abbr: BondAbbr<'p>, donor: NodeId) -> Self {
        let end = Terminal::Acceptor(abbr);
        Self { end, target: donor }
    }
}

// PERF: Not sold on this "tombstone" approach with `Option` — might be better to re-index the sub-graphs so their
// vectors can be shrunk?
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
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

type NodeId = usize;

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
    fn index(&self, id: &ResidueId) -> NodeId {
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
        for piece in fragment.fragment(Some(1)) {
            eprintln!("{piece}");
        }
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
                // FIXME: Shocking failure of type inference here with that `NodeId`...
                |node: NodeId, bond| residues[node].as_mut().unwrap().bonds.push(bond);

            push_bond(donor, Bond::donating_to(abbr, acceptor));
            push_bond(acceptor, Bond::accepting_from(abbr, donor));
        }

        Self(residues)
    }

    fn fragment(&self, max_depth: Option<usize>) -> HashSet<Self> {
        let mut processing_queue = vec![self.clone()];
        // PERF: Any clever `with_capacity` pre-allocation I could do?
        let mut fragments = HashSet::new();

        let mut depth = 0;
        while let Some(next) = processing_queue.pop() {
            // FIXME: Can I avoid that clone?
            if fragments.insert(next.clone()) {
                continue;
            }

            if depth < max_depth.unwrap_or(usize::MAX) {
                depth += 1;
                fragments.extend(next.cut_each_bond().map(|(_, piece)| piece));
            }
        }

        fragments
    }

    // FIXME: Clarify type with some aliases?
    // FIXME: Remove the `+ '_` once Rust 2024 is released!
    fn cut_each_bond(&self) -> impl Iterator<Item = (NodeId, Self)> + '_ {
        self.0
            .iter()
            .enumerate()
            .filter_map(|(node, opt_residue)| opt_residue.as_ref().map(|residue| (node, residue)))
            .flat_map(|(node, residue)| residue.bonds.iter().map(move |bond| (node, bond.target)))
            // NOTE: Ensures that the same bonds aren't cut twice
            .filter(|(a, b)| a < b)
            .map(|(a, b)| {
                let mut fragment = self.clone();
                let mut remove_edge = |from, to| {
                    swap_remove_first(&mut fragment[from].bonds, |bond| bond.target == to).unwrap()
                };

                let terminal_a = remove_edge(a, b).end;
                let terminal_b = remove_edge(b, a).end;
                fragment[a].terminals.push(terminal_a);
                fragment[b].terminals.push(terminal_b);

                (a, fragment)
            })
    }
}

impl<'p> Index<NodeId> for Fragment<'p> {
    type Output = Residue<'p>;

    fn index(&self, index: NodeId) -> &Self::Output {
        self.0[index].as_ref().unwrap()
    }
}

impl<'p> IndexMut<NodeId> for Fragment<'p> {
    fn index_mut(&mut self, index: NodeId) -> &mut Self::Output {
        self.0[index].as_mut().unwrap()
    }
}

// FIXME: Should this just unwrap things here? Assuming there will always be one match?
fn swap_remove_first<T>(vec: &mut Vec<T>, f: impl FnMut(&T) -> bool) -> Option<T> {
    vec.iter().position(f).map(|i| vec.swap_remove(i))
}

// SEE NOTES FROM APRIL 8TH!
// use DashMap or quick-cache for a global fragment cache

// FIXME: Tests need adding!
