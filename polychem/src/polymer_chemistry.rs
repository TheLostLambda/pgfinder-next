use std::collections::HashMap;

use knuffel::{
    span::{Span, Spanned},
    Decode,
};
use miette::{Diagnostic, NamedSource, Result, SourceSpan};
use thiserror::Error;

use crate::{atomic_database::AtomicDatabase, ChemicalComposition, FunctionalGroup};

// FIXME: Might want to change how this is structured down the line...
// FIXME: Accessors, not public fields!
#[derive(Clone, PartialEq, Eq, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
struct PolymerChemistry {
    bonds: HashMap<String, BondDescription>,
    modifications: HashMap<String, ModificationDescription>,
    residues: HashMap<String, ResidueDescription>,
}

#[derive(Clone, PartialEq, Eq, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
struct BondDescription {
    from: Target,
    to: Target,
    lost: ChemicalComposition,
}

#[derive(Clone, PartialEq, Eq, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
struct ModificationDescription {
    name: String,
    lost: ChemicalComposition,
    gained: ChemicalComposition,
    targets: Vec<Target>,
}

#[derive(Clone, PartialEq, Eq, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
struct ResidueDescription {
    name: String,
    composition: ChemicalComposition,
    functional_groups: Vec<FunctionalGroup>,
}

#[derive(Clone, PartialEq, Eq, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
struct Target {
    functional_group: String,
    location: Option<String>,
    residue: Option<String>,
}

struct ProtoChemistryError(SourceSpan, ChemistryErrorKind);

impl ProtoChemistryError {
    fn finalize(self, file_name: impl AsRef<str>, kdl: impl AsRef<str>) -> ChemistryError {
        let Self(span, kind) = self;
        let kdl = NamedSource::new(file_name, kdl.as_ref().to_string());
        ChemistryError { kdl, span, kind }.into()
    }
}

impl PolymerChemistry {
    pub fn from_kdl(
        db: &AtomicDatabase,
        file_name: impl AsRef<str>,
        kdl: impl AsRef<str>,
    ) -> Result<Self> {
        let parsed_chemistry: PolymerChemistryKdl =
            knuffel::parse(file_name.as_ref(), kdl.as_ref())?;
        parsed_chemistry
            .validate(db)
            .map_err(|e| e.finalize(file_name, kdl).into())
    }
}
type ChemResult<T> = Result<T, ProtoChemistryError>;

impl PolymerChemistryKdl {
    fn validate(self, db: &AtomicDatabase) -> ChemResult<PolymerChemistry> {
        Ok(PolymerChemistry {
            bonds: self.bonds.validate(db)?,
            modifications: self.modifications.validate(db)?,
            residues: self.residues.validate(db)?,
        })
    }
}

impl ResiduesKdl {
    fn validate(self, db: &AtomicDatabase) -> ChemResult<HashMap<String, ResidueDescription>> {
        let types: HashMap<_, _> = self.types.into_iter().map(ResidueTypeEntry::from).collect();
        self.residues
            .into_iter()
            // FIXME: Big hmm for the ? operator here...
            .map(|r| Ok((r.abbr.clone(), r.validate(db, &types)?)))
            .collect()
    }
}

impl ResidueKdl {
    fn validate(
        self,
        db: &AtomicDatabase,
        types: &HashMap<String, Vec<FunctionalGroup>>,
    ) -> ChemResult<ResidueDescription> {
        let mut functional_groups: Vec<_> = self
            .functional_groups
            .into_iter()
            .map(FunctionalGroup::from)
            .collect();
        functional_groups.extend(types[&self.residue_type].clone());
        Ok(ResidueDescription {
            name: self.name,
            composition: ChemicalCompositionKdl::validate_optional(self.composition, db)?,
            functional_groups,
        })
    }
}

impl ChemicalCompositionKdl {
    fn validate(self, db: &AtomicDatabase) -> ChemResult<ChemicalComposition> {
        ChemicalComposition::new(db, &*self.0)
            .map_err(|e| ProtoChemistryError(self.0.span().clone().into(), e.into()))
    }

    fn validate_optional(
        value: Option<Self>,
        db: &AtomicDatabase,
    ) -> ChemResult<ChemicalComposition> {
        // FIXME: Good lord, maybe just if let?
        Ok(value
            .map(|c| c.validate(db))
            .transpose()?
            .unwrap_or_default())
    }
}

type ResidueTypeEntry = (String, Vec<FunctionalGroup>);

impl From<ResidueTypeKdl> for ResidueTypeEntry {
    fn from(value: ResidueTypeKdl) -> Self {
        let functional_groups = value
            .functional_groups
            .into_iter()
            .map(FunctionalGroup::from)
            .collect();
        (value.name, functional_groups)
    }
}

impl From<FunctionalGroupKdl> for FunctionalGroup {
    fn from(value: FunctionalGroupKdl) -> Self {
        Self {
            name: value.name,
            location: value.location,
        }
    }
}

impl ModificationsKdl {
    fn validate(self, db: &AtomicDatabase) -> ChemResult<HashMap<String, ModificationDescription>> {
        // FIXME: DRY this out...
        self.modifications
            .into_iter()
            .map(|m| Ok((m.abbr.clone(), m.validate(db)?)))
            .collect()
    }
}

impl ModificationKdl {
    fn validate(self, db: &AtomicDatabase) -> ChemResult<ModificationDescription> {
        Ok(ModificationDescription {
            name: self.name,
            lost: ChemicalCompositionKdl::validate_optional(self.lost, db)?,
            gained: ChemicalCompositionKdl::validate_optional(self.gained, db)?,
            targets: self.targets.into_iter().map(Target::from).collect(),
        })
    }
}
impl BondsKdl {
    fn validate(self, db: &AtomicDatabase) -> ChemResult<HashMap<String, BondDescription>> {
        // FIXME: This is a bit repetitive with the above...
        self.bonds
            .into_iter()
            .map(|b| Ok((b.name.clone(), b.validate(db)?)))
            .collect()
    }
}

impl BondKdl {
    fn validate(self, db: &AtomicDatabase) -> ChemResult<BondDescription> {
        Ok(BondDescription {
            from: self.from.into(),
            to: self.to.into(),
            lost: self.lost.validate(db)?,
        })
    }
}

impl From<TargetKdl> for Target {
    fn from(value: TargetKdl) -> Self {
        Self {
            functional_group: value.functional_group,
            location: value.location,
            residue: value.residue,
        }
    }
}

// FIXME: Check that field names here line up with those in `lib.rs`!
#[derive(Decode, Debug)]
#[knuffel(span_type=Span)]
struct PolymerChemistryKdl {
    #[knuffel(child)]
    bonds: BondsKdl,
    #[knuffel(child)]
    modifications: ModificationsKdl,
    #[knuffel(child)]
    residues: ResiduesKdl,
}

#[derive(Decode, Debug)]
#[knuffel(span_type=Span)]
struct BondsKdl {
    #[knuffel(children)]
    bonds: Vec<BondKdl>,
}

#[derive(Decode, Debug)]
#[knuffel(span_type=Span)]
struct ChemicalCompositionKdl(#[knuffel(argument)] Spanned<String, Span>);

#[derive(Decode, Debug)]
#[knuffel(span_type=Span)]
struct BondKdl {
    #[knuffel(node_name)]
    name: String,
    #[knuffel(child)]
    from: TargetKdl,
    #[knuffel(child)]
    to: TargetKdl,
    #[knuffel(child)]
    lost: ChemicalCompositionKdl,
}

#[derive(Decode, Debug)]
struct TargetKdl {
    #[knuffel(argument)]
    functional_group: String,
    #[knuffel(property(name = "at"))]
    location: Option<String>,
    #[knuffel(property(name = "of"))]
    residue: Option<String>,
}

#[derive(Decode, Debug)]
#[knuffel(span_type=Span)]
struct ModificationsKdl {
    #[knuffel(children)]
    modifications: Vec<ModificationKdl>,
}

#[derive(Decode, Debug)]
#[knuffel(span_type=Span)]
struct ModificationKdl {
    #[knuffel(node_name)]
    abbr: String,
    #[knuffel(argument)]
    name: String,
    #[knuffel(child)]
    lost: Option<ChemicalCompositionKdl>,
    #[knuffel(child)]
    gained: Option<ChemicalCompositionKdl>,
    #[knuffel(children(name = "targeting"))]
    targets: Vec<TargetKdl>,
}

#[derive(Decode, Debug)]
#[knuffel(span_type=Span)]
struct ResiduesKdl {
    #[knuffel(child, unwrap(children))]
    types: Vec<ResidueTypeKdl>,
    #[knuffel(children)]
    residues: Vec<ResidueKdl>,
}

#[derive(Decode, Debug)]
struct ResidueTypeKdl {
    #[knuffel(node_name)]
    name: String,
    #[knuffel(children)]
    functional_groups: Vec<FunctionalGroupKdl>,
}

// FIXME: Be sure to write tests to check this errors for missing composition keys but works for nulled ones — also
// check that the `lost` and `gained` keys are the opposite — null isn't allowed, but they can be left out entirely
// NOTE: This forces composition to have a `null` value instead of being a totally optional key
// #[derive(Decode, Debug)]
// #[knuffel(span_type=Span)]
// struct NullOr<T: Decode<Span>>(#[knuffel(argument)] Option<T>);

#[derive(Decode, Debug)]
#[knuffel(span_type=Span)]
struct ResidueKdl {
    #[knuffel(node_name)]
    residue_type: String,
    #[knuffel(argument)]
    abbr: String,
    #[knuffel(argument)]
    name: String,
    #[knuffel(child)]
    composition: Option<ChemicalCompositionKdl>,
    #[knuffel(children(name = "functional-group"))]
    functional_groups: Vec<FunctionalGroupKdl>,
}

#[derive(Decode, Debug)]
struct FunctionalGroupKdl {
    #[knuffel(argument)]
    name: String,
    #[knuffel(property(name = "at"))]
    location: String,
}

// FIXME: Might also want to pull out this error-handling into it's own crate...
#[derive(Debug, Error, Diagnostic)]
#[error("failed to validate polymer chemistry file")]
struct ChemistryError {
    #[source_code]
    kdl: NamedSource<String>,
    #[label("{}", kind.label())]
    span: SourceSpan,
    #[source]
    #[diagnostic_source]
    kind: ChemistryErrorKind,
}

// FIXME: Need to implement something like that labeled error trait...
#[derive(Debug, Clone, Error, Diagnostic)]
enum ChemistryErrorKind {
    #[error(transparent)]
    #[diagnostic(transparent)]
    Composition(#[from] crate::Error),
}

impl ChemistryErrorKind {
    const fn label(&self) -> &'static str {
        match self {
            Self::Composition(_) => "invalid chemical composition",
        }
    }
}

#[cfg(test)]
mod tests {
    use indoc::indoc;
    use insta::{assert_debug_snapshot, assert_ron_snapshot};
    use miette::Result;
    use once_cell::sync::Lazy;

    use crate::{
        atomic_database::AtomicDatabase, polymer_chemistry::ResidueDescription,
        testing_tools::assert_miette_snapshot,
    };

    use super::{PolymerChemistry, PolymerChemistryKdl, ResiduesKdl};

    static DB: Lazy<AtomicDatabase> = Lazy::new(|| {
        AtomicDatabase::from_kdl(
            "atomic_database.kdl",
            include_str!("../atomic_database.kdl"),
        )
        .unwrap()
    });

    const KDL: &str = include_str!("../muropeptide_chemistry.kdl");

    // FIXME: Be sure to test that miette errors from the composition parser get displayed! `[100Tc]` etc...
    #[test]
    fn parse_muropeptide_chemistry() {
        let chemistry: PolymerChemistryKdl =
            knuffel::parse("muropeptide_chemistry.rs", KDL).unwrap();
        assert_debug_snapshot!(chemistry);
    }

    #[test]
    fn build_muropeptide_chemistry() {
        let chemistry = PolymerChemistry::from_kdl(&DB, "muropeptide_chemistry.rs", KDL).unwrap();
        assert_ron_snapshot!(chemistry, {
            ".bonds, .modifications, .residues" => insta::sorted_redaction(),
            ".**.isotopes" => insta::sorted_redaction()
        });
    }

    #[test]
    fn parse_residues_types_without_groups() -> Result<()> {
        let kdl = indoc! {r#"
            types {
                AminoAcid
            }
            AminoAcid "G" "Glycine" {
                composition "C2H5NO2"
            }
        "#};
        let residues: Result<ResiduesKdl, _> = knuffel::parse("test", kdl);
        assert!(residues.is_ok());
        assert_debug_snapshot!(residues.unwrap());
        Ok(())
    }

    #[test]
    fn parse_residues_without_types() -> Result<()> {
        let kdl = indoc! {r#"
            types {
            //     AminoAcid
            }
            AminoAcid "G" "Glycine" {
                composition "C2H5NO2"
            }
        "#};
        let residues: ResiduesKdl = knuffel::parse("test", kdl).unwrap();
        // let residues: HashMap<String, ResidueDescription> = residues.ctx_try_into(&*DB)?;
        // assert!(residues.is_err());
        // assert_miette_snapshot!(residues);
        Ok(())
    }
}
