#![type_length_limit = "18554191"]

use std::{
    cmp::{self, Ordering},
    convert::identity,
    fmt::{self, Display, Formatter},
    hash::{Hash, Hasher},
    mem,
    ops::{Index, IndexMut},
};

use ahash::{HashSet, HashSetExt};
use derive_more::IsVariant;
use itertools::Itertools;
use polychem::{BondId, BondInfo, OffsetKind, Polymer, ResidueGroup, ResidueId};

// FIXME: Consider using newtype? Especially if this is made public!
type BondAbbr<'p> = &'p str;

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, IsVariant, Debug)]
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
    // FIXME: Since these terminal lists are typically super small, I'm using a `Vec` of sorted tuples instead of a
    // proper map type, but this needs real benchmarking to see if that makes sense!
    terminals: Vec<(Terminal<'p>, u32)>,
    bonds: Vec<Bond<'p>>,
}

#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
struct Bond<'p> {
    // FIXME: Not in love with how `BondId` tracking is done — does this perform well?
    id: BondId,
    // FIXME: The naming `end` doesn't fit well with the `Terminal` here...
    end: Terminal<'p>,
    target: NodeId,
}

impl<'p> Bond<'p> {
    const fn donating_to(id: BondId, abbr: BondAbbr<'p>, acceptor: NodeId) -> Self {
        let end = Terminal::Donor(abbr);
        Self {
            id,
            end,
            target: acceptor,
        }
    }

    const fn accepting_from(id: BondId, abbr: BondAbbr<'p>, donor: NodeId) -> Self {
        let end = Terminal::Acceptor(abbr);
        Self {
            id,
            end,
            target: donor,
        }
    }
}

// PERF: Not sold on this "tombstone" approach with `Option` — might be better to re-index the sub-graphs so their
// vectors can be shrunk?
#[derive(Clone, Debug)]
struct Fragment<'p> {
    broken_bonds: Vec<BondId>,
    residues: Vec<Option<Residue<'p>>>,
}

// NOTE: Whilst we do care about tracking which bonds are broken to generate a fragment, two fragments should be
// considered the same if their residues (and bonds) are the same, regardless of either fragment's history. These `Eq`
// and `Hash` implementations ignore the `broken_bonds` field when checking equality!
impl Eq for Fragment<'_> {}
impl PartialEq for Fragment<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.residues == other.residues
    }
}

impl Hash for Fragment<'_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.residues.hash(state);
    }
}

// FIXME: This was written way too quickly... Take a look back at this...
// FIXME: Also maybe should die, as it's only useful for debugging?
impl Display for Fragment<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        writeln!(f, "digraph {{")?;
        for (id, residue) in self.residues.iter().enumerate() {
            let Some(Residue { terminals, bonds }) = residue else {
                continue;
            };

            let terminals = terminals
                .iter()
                .map(|(terminal, count)| format!("{count}x{terminal}"))
                .join(", ");
            writeln!(f, r#"  {id} [label="{id}\n[{terminals}]"]"#)?;

            // NOTE: The `.filter()` here prevents the double-printing of edges
            for Bond { end, target, .. } in bonds.iter().filter(|b| id < b.target) {
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
        // PERF: Linear search might be faster here for small lists!
        // SAFETY: Panics if the `id` isn't found
        self.0.binary_search(id).unwrap()
    }
}

// DESIGN: This `Sized` bound could be moved to the method-level with `where Self: Sized`, which would allow this trait
// to be implemented for unsized types, just without those methods, but I don't think this trait makes much sense
// without a `.fragment()` method (which currently requires a `Self: Sized`), so I've elected to make the whole trait
// `Sized` for now
pub trait Dissociable: Sized {
    // TODO: Add a method for checking if Self's polymer is currently in one piece or if it's already disconnected
    #[must_use]
    fn polymer(&self) -> &Polymer;
    // FIXME: Naming?
    #[must_use]
    fn new_fragment(
        &self,
        fragmented_polymer: Polymer,
        lost_residues: Vec<ResidueId>,
        broken_bonds: Vec<BondId>,
    ) -> Self;
    #[must_use]
    fn fragment(&self) -> Vec<Self> {
        todo!()
    }
    // FIXME: Only for testing! Remove!
    fn dbg_fragment(&self) {
        let polymer = self.polymer();
        let node_mapping = NodeMapping::new(polymer);
        let fragment = Fragment::new(&node_mapping, polymer);
        eprintln!("\nPieces: {}", fragment.fragment(None).len());
    }
}

impl<'p> Fragment<'p> {
    // FIXME: This currently assumes that the `polymer` provided consists of a single connected component. If this is
    // not the case, I should report an error to the user!
    fn new(node_mapping: &NodeMapping, polymer: &Polymer<'_, 'p>) -> Self {
        let mut residues = vec![Some(Residue::default()); node_mapping.0.len()];

        for (id, BondInfo(ResidueGroup(donor, _), bond, ResidueGroup(acceptor, _))) in
            polymer.bonds()
        {
            let abbr = bond.abbr();
            let donor = node_mapping.index(donor);
            let acceptor = node_mapping.index(acceptor);

            let mut push_bond =
                // FIXME: Shocking failure of type inference here with that `NodeId`...
                |node: NodeId, bond| residues[node].as_mut().unwrap().bonds.push(bond);

            push_bond(donor, Bond::donating_to(id, abbr, acceptor));
            push_bond(acceptor, Bond::accepting_from(id, abbr, donor));
        }

        Self {
            broken_bonds: Vec::new(),
            residues,
        }
    }

    // FIXME: Naming?
    fn build_fragment_ion<'a>(
        self,
        node_mapping: &NodeMapping,
        polymer: &Polymer<'a, 'p>,
    ) -> Polymer<'a, 'p> {
        let mut fragmented_polymer = polymer.clone();

        for (id, opt_residue) in self.residues.into_iter().enumerate() {
            let residue_id = node_mapping.0[id];
            if let Some(residue) = opt_residue {
                // FIXME: This is hacky, hard-coded logic for PG — this will need to be generalized for other molecules!
                let donor_count: u32 = residue
                    .terminals
                    .into_iter()
                    .filter_map(|(t, c)| t.is_donor().then_some(c))
                    .sum();
                if donor_count != 0 {
                    // SAFETY: This should only panic if `residue_id` doesn't exist (which it should), or if `H2O` isn't a
                    // valid chemical formula (but it is)
                    fragmented_polymer
                        .offset_residue(OffsetKind::Remove, donor_count, "H2O", residue_id)
                        .unwrap();
                }
            } else {
                // FIXME: Need to build up a list of the lost residues and broken bonds to return!
                fragmented_polymer.remove_residue(residue_id);
            }
        }

        for id in self.broken_bonds {
            // FIXME: Need to build up a list of the lost residues and broken bonds to return!
            fragmented_polymer.remove_bond(id);
        }

        // FIXME: Hard-coded to generate the 1+ ions!
        // SAFETY: `p` is a valid formula, so this shouldn't panic
        fragmented_polymer
            .new_offset(OffsetKind::Add, 1, "p")
            .unwrap();
        fragmented_polymer
    }

    // FIXME: I'm pretty sure this assumption holds, but does a fragmentation depth equalling the degree / valency of
    // the graph mean 100% of fragments will be generated at that depth? If so, this saves a lot of time wasted on
    // "deeper" fragmentations that have no chance of generating a unique fragment!
    fn degree(&self) -> usize {
        self.residues
            .iter()
            .flatten()
            .map(|residue| residue.bonds.len())
            .max()
            .unwrap_or_default()
    }

    // PERF: This could be memoized, but I don't know how much that would gain me if I'm not adding sub-fragmentations
    // to the cache... It would only help if the user has duplicate structures they are asking to fragment!
    fn fragment(&self, max_depth: Option<usize>) -> HashSet<Self> {
        let mut processing = vec![self.clone()];
        let mut processing_queue = Vec::new();
        // PERF: Any clever `with_capacity` pre-allocation I could do?
        let mut fragments = HashSet::new();

        let max_depth = max_depth.unwrap_or(usize::MAX);
        let useful_depth = cmp::min(self.degree(), max_depth);

        for depth in 0..=useful_depth {
            while let Some(next) = processing.pop() {
                // FIXME: Can I avoid this clone?
                if fragments.insert(next.clone()) && depth < useful_depth {
                    processing_queue.extend(next.cut_each_bond().flat_map(divide_fragment));
                }
            }
            if processing_queue.is_empty() {
                break;
            }
            // FIXME: These names could be improved a lot — also, is this a double buffering?
            mem::swap(&mut processing, &mut processing_queue);
        }

        fragments
    }

    // FIXME: Clarify type with some aliases?
    // FIXME: Remove the `+ '_` once Rust 2024 is released!
    // FIXME: Should this really be a method?
    fn cut_each_bond(&self) -> impl Iterator<Item = (NodeId, Self)> + '_ {
        self.residues
            .iter()
            .enumerate()
            .filter_map(|(node, opt_residue)| opt_residue.as_ref().map(|residue| (node, residue)))
            .flat_map(|(node, residue)| residue.bonds.iter().map(move |bond| (node, bond.target)))
            // FIXME: Perhaps here is where I can eventually filter out any unfragmentable bond types!
            // NOTE: Ensures that the same bonds aren't cut twice
            .filter(|(a, b)| a < b)
            .map(|(a, b)| {
                let mut fragment = self.clone();
                let mut remove_edge = |from, to| {
                    swap_remove_first(&mut fragment[from].bonds, |bond| bond.target == to).unwrap()
                };

                let Bond {
                    id: bond_id,
                    end: terminal_a,
                    ..
                } = remove_edge(a, b);
                let terminal_b = remove_edge(b, a).end;

                // FIXME: This has some pretty bad complexity because the vector may need to be shifted, but I need to
                // make sure that, regardless of the order bonds are broken, two terminals lists are considered the
                // same if they contain the same items!
                sorted_insert(&mut fragment.broken_bonds, bond_id);
                // FIXME: These have the same insertion-shifting issue!
                insert_terminal(&mut fragment[a].terminals, terminal_a);
                insert_terminal(&mut fragment[b].terminals, terminal_b);

                (a, fragment)
            })
    }
}

// FIXME: This is performing a linear search (instead of a binary one) since these vectors should be super small most
// of the time. That said, I've *not* benchmarked things properly!!!
fn insert_terminal<'p>(terminals: &mut Vec<(Terminal<'p>, u32)>, terminal: Terminal<'p>) {
    // FIXME: Probably a more clever way to do this...
    let len = terminals.len();
    for i in (0..len).rev() {
        match terminals[i].0.cmp(&terminal) {
            Ordering::Equal => terminals[i].1 += 1,
            Ordering::Less => terminals.insert(i + 1, (terminal, 1)),
            Ordering::Greater => continue,
        }
        return;
    }
    terminals.insert(0, (terminal, 1));
}

// FIXME: This is performing a linear search (instead of a binary one) since these vectors should be super small most
// of the time. That said, I've *not* benchmarked things properly!!! Try `partition_point()` when benchmarking the
// binary search approach!
fn sorted_insert<T: PartialOrd>(vec: &mut Vec<T>, element: T) {
    if let Some(index) = vec.iter().rposition(|x| x <= &element) {
        vec.insert(index + 1, element);
    } else {
        vec.insert(0, element);
    }
}

// FIXME: Should this really be a bare function?
// FIXME: Do I need the lifetimes in `Fragment`?
fn divide_fragment((cut_node, mut fragment): (NodeId, Fragment)) -> Vec<Fragment> {
    // FIXME: Is this at all a good idea? It's more space for faster lookups?
    let mut visited = vec![false; fragment.residues.len()];
    let mut stack = vec![cut_node];

    while let Some(next) = stack.pop() {
        visited[next] = true;
        stack.extend(
            fragment[next]
                .bonds
                .iter()
                .map(|bond| bond.target)
                .filter(|&node| !visited[node]),
        );
    }

    // FIXME: Does this special case actually save me any time or effort?
    if visited.iter().copied().all(identity) {
        vec![fragment]
    } else {
        let mut rest = fragment.clone();
        for (node, visited) in visited.into_iter().enumerate() {
            if visited {
                rest.residues[node] = None;
            } else {
                fragment.residues[node] = None;
            }
        }
        vec![fragment, rest]
    }
}

// FIXME: Honestly, these might be making the code more confusing... Consider removing or replacing with a method?
impl<'p> Index<NodeId> for Fragment<'p> {
    type Output = Residue<'p>;

    fn index(&self, index: NodeId) -> &Self::Output {
        self.residues[index].as_ref().unwrap()
    }
}

// FIXME: Honestly, these might be making the code more confusing... Consider removing or replacing with a method?
impl<'p> IndexMut<NodeId> for Fragment<'p> {
    fn index_mut(&mut self, index: NodeId) -> &mut Self::Output {
        self.residues[index].as_mut().unwrap()
    }
}

// FIXME: Should this just unwrap things here? Assuming there will always be one match?
fn swap_remove_first<T>(vec: &mut Vec<T>, f: impl FnMut(&T) -> bool) -> Option<T> {
    vec.iter().position(f).map(|i| vec.swap_remove(i))
}

// SEE NOTES FROM APRIL 8TH!
// use DashMap or quick-cache for a global fragment cache

// FIXME: Tests need adding!
