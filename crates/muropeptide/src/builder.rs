use polychem::{Polymer, ResidueId};

use crate::{Monomer, Multimer, ResidueAbbr};

// FIXME: Need to think about if these should really live in another KDL config?
const PEPTIDE_BOND: &str = "Pep";
const GLYCOSIDIC_BOND: &str = "Gly";
const STEM_BOND: &str = "Stem";

#[derive(Debug)]
pub enum BuilderError {}
type BuilderResult<T> = Result<T, BuilderError>;

pub trait Build<'m, 'a, 'p> {
    type Into;

    fn build(self, polymer: &'m mut Polymer<'a, 'p>) -> BuilderResult<Self::Into>;
}

impl<'m, 'a, 'p> Build<'m, 'a, 'p> for Multimer<ResidueAbbr<'_>> {
    type Into = Multimer<ResidueId>;

    fn build(self, polymer: &'m mut Polymer<'a, 'p>) -> BuilderResult<Self::Into> {
        Ok(Multimer {
            monomers: self.monomers.build(polymer)?,
            // FIXME: Replace these with real implementations!
            connections: Vec::new(),
            modifications: Vec::new(),
        })
    }
}

// FIXME: Maybe get rid of this and move the mapping into `Multimer`? Or `Monomer`?
impl<'m, 'a, 'p> Build<'m, 'a, 'p> for Vec<Monomer<ResidueAbbr<'_>>> {
    type Into = Vec<Monomer<ResidueId>>;

    fn build(self, polymer: &'m mut Polymer<'a, 'p>) -> BuilderResult<Self::Into> {
        self.into_iter().map(|m| m.build(polymer)).collect()
    }
}

impl<'m, 'a, 'p> Build<'m, 'a, 'p> for Monomer<ResidueAbbr<'_>> {
    type Into = Monomer<ResidueId>;

    fn build(self, polymer: &'m mut Polymer<'a, 'p>) -> BuilderResult<Self::Into> {
        // FIXME: Move this to it's own build block once `peptide` includes lateral chains and has its own type...
        // FIXME: Replace this `.unwrap()` with `?`
        let glycan: Vec<_> = polymer
            .new_chain(GLYCOSIDIC_BOND, &self.glycan)
            .unwrap()
            .0
            .collect();

        // FIXME: Move this to it's own build block once `peptide` includes lateral chains and has its own type...
        // FIXME: Replace this `.unwrap()` with `?`
        let peptide: Vec<_> = polymer
            .new_chain(PEPTIDE_BOND, &self.peptide)
            .unwrap()
            .0
            .collect();

        if let (Some(&donor), Some(&acceptor)) = (glycan.last(), peptide.first()) {
            // FIXME: Properly handle this error!
            polymer.bond_residues(STEM_BOND, donor, acceptor).unwrap();
        }

        Ok(Monomer { glycan, peptide })
    }
}
