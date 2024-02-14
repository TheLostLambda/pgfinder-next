// FIXME: Create a sub-module with this parser / validator + the Polymerizer + the TargetIndex modules
use std::collections::HashMap;

use knuffel::{
    span::{Span, Spanned},
    Decode,
};
use miette::{Diagnostic, LabeledSpan, NamedSource, Result};
use thiserror::Error;

use crate::{
    atomic_database::AtomicDatabase, chemical_targets::Target, ChemicalComposition, FunctionalGroup,
};

// FIXME: Accessors, not public fields!
type Bonds = HashMap<String, BondDescription>;
type Modifications = HashMap<String, ModificationDescription>;
type Residues = HashMap<String, ResidueDescription>;

#[derive(Clone, PartialEq, Eq, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
struct PolymerChemistry {
    bonds: Bonds,
    modifications: Modifications,
    residues: Residues,
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

trait ValidateInto<'a, T> {
    type Context: 'a;

    fn validate(self, ctx: Self::Context) -> ChemResult<T>;
}

impl<'a> ValidateInto<'a, ChemicalComposition> for Option<ChemicalCompositionKdl> {
    type Context = &'a AtomicDatabase;

    fn validate(self, ctx: Self::Context) -> ChemResult<ChemicalComposition> {
        Ok(if let Some(c) = self {
            c.validate(ctx)?
        } else {
            ChemicalComposition::default()
        })
    }
}
impl PolymerChemistry {
    pub fn from_kdl(
        db: &AtomicDatabase,
        file_name: impl AsRef<str>,
        text: impl AsRef<str>,
    ) -> Result<Self> {
        let parsed_chemistry: PolymerChemistryKdl =
            knuffel::parse(file_name.as_ref(), text.as_ref())?;
        parsed_chemistry
            .validate(db)
            .map_err(|e| e.finalize(file_name, text).into())
    }
}

type ChemResult<T> = Result<T, ChemistryErrorKind>;

impl<'a> ValidateInto<'a, PolymerChemistry> for PolymerChemistryKdl {
    type Context = &'a AtomicDatabase;

    fn validate(self, ctx: Self::Context) -> ChemResult<PolymerChemistry> {
        let residues = self.residues.validate(ctx)?;
        Ok(PolymerChemistry {
            bonds: self.bonds.validate(ctx)?,
            modifications: self.modifications.validate((ctx, &residues))?,
            residues,
        })
    }
}

impl<'a> ValidateInto<'a, Residues> for ResiduesKdl {
    type Context = &'a AtomicDatabase;

    fn validate(self, ctx: Self::Context) -> ChemResult<Residues> {
        let types = self.types.validate(())?;
        self.residues
            .into_iter()
            // FIXME: Big hmm for the ? operator here...
            .map(|r| Ok((r.abbr.clone(), r.validate((ctx, &types))?)))
            .collect()
    }
}

impl<'a> ValidateInto<'a, HashMap<String, Vec<FunctionalGroupKdl>>> for Vec<ResidueTypeKdl> {
    type Context = ();

    fn validate(self, _ctx: Self::Context) -> ChemResult<HashMap<String, Vec<FunctionalGroupKdl>>> {
        // FIXME: Is there a more functional way to do this?
        let mut types = HashMap::with_capacity(self.len());
        for residue_type in self {
            if let Some((first_span, _)) = types.insert(
                residue_type.name.clone(),
                (residue_type.span, residue_type.functional_groups),
            ) {
                return Err(ChemistryErrorKind::DuplicateResidueType(
                    first_span,
                    residue_type.span,
                    residue_type.name,
                ));
            }
        }
        Ok(types.into_iter().map(|(k, (_, v))| (k, v)).collect())
    }
}

impl<'a> ValidateInto<'a, ResidueDescription> for ResidueKdl {
    type Context = (
        &'a AtomicDatabase,
        &'a HashMap<String, Vec<FunctionalGroupKdl>>,
    );

    fn validate(self, ctx: Self::Context) -> ChemResult<ResidueDescription> {
        let type_groups = ctx
            .1
            .get(&self.residue_type)
            .ok_or_else(|| ChemistryErrorKind::UndefinedResidueType(self.span, self.residue_type))?
            .clone();
        // FIXME: Abstract this logic out into a function — it's used in a couple of places, the dedeuplication
        let mut functional_group_spans = HashMap::new();
        for functional_group_kdl in type_groups.into_iter().chain(self.functional_groups) {
            let span = functional_group_kdl.span;
            let functional_group = FunctionalGroup::from(functional_group_kdl);
            if let Some(first_defined_span) =
                functional_group_spans.insert(functional_group.clone(), span)
            {
                return Err(ChemistryErrorKind::DuplicateFunctionalGroup(
                    first_defined_span,
                    span,
                    functional_group.name,
                    functional_group.location,
                ));
            }
        }
        let functional_groups = functional_group_spans
            .into_keys()
            .map(FunctionalGroup::from)
            .collect();
        Ok(ResidueDescription {
            name: self.name,
            composition: self.composition.validate(ctx.0)?,
            functional_groups,
        })
    }
}

// FIXME: Might also want to pull out this error-handling into it's own crate...
#[derive(Debug, Error)]
#[error("failed to validate polymer chemistry file")]
struct ChemistryError {
    kdl: NamedSource<String>,
    #[source]
    kind: ChemistryErrorKind,
}

impl Diagnostic for ChemistryError {
    fn source_code(&self) -> Option<&dyn miette::SourceCode> {
        Some(&self.kdl)
    }

    fn labels(&self) -> Option<Box<dyn Iterator<Item = miette::LabeledSpan> + '_>> {
        Some(Box::new(self.kind.labels().into_iter().map(|(s, l)| {
            LabeledSpan::new_with_span(Some(l.to_string()), s)
        })))
    }

    fn diagnostic_source(&self) -> Option<&dyn Diagnostic> {
        Some(&self.kind)
    }
}

#[derive(Debug, Clone, Error, Diagnostic)]
enum ChemistryErrorKind {
    #[error("polymer chemistry file contained an invalid chemical composition")]
    Composition(
        Span,
        // FIXME: I've checked with cargo expand — I need both source and diagnostic source. Make I've got this
        // #[source] added globally — wherever `diagnostic_source` is defined!
        #[source]
        #[diagnostic_source]
        crate::Error,
    ),
    #[error("the residue type {1:?} is undefined")]
    #[diagnostic(help("double-check for typos, or add {1:?} to the types section"))]
    UndefinedResidueType(Span, String),
    #[error("the residue type {2:?} has already been defined")]
    #[diagnostic(help("consider consolidating duplicate types or picking a new type name"))]
    DuplicateResidueType(Span, Span, String),
    #[error("the functional group {2:?} has already been defined at {3:?}")]
    #[diagnostic(help("double-check for typos, or remove the duplicate functional group"))]
    DuplicateFunctionalGroup(Span, Span, String, String),
}

impl ChemistryErrorKind {
    fn labels(&self) -> Vec<(Span, &'static str)> {
        match self {
            Self::Composition(s, _) => vec![(*s, "invalid chemical composition")],
            Self::UndefinedResidueType(s, _) => vec![(*s, "undefined residue type")],
            Self::DuplicateResidueType(s1, s2, _)
            | Self::DuplicateFunctionalGroup(s1, s2, _, _) => {
                vec![(*s1, "first defined here"), (*s2, "then again here")]
            }
        }
    }

    fn finalize(self, file_name: impl AsRef<str>, kdl: impl AsRef<str>) -> ChemistryError {
        let kdl = NamedSource::new(file_name, kdl.as_ref().to_string());
        ChemistryError { kdl, kind: self }
    }
}

impl<'a> ValidateInto<'a, ChemicalComposition> for ChemicalCompositionKdl {
    type Context = &'a AtomicDatabase;

    fn validate(self, ctx: Self::Context) -> ChemResult<ChemicalComposition> {
        ChemicalComposition::new(ctx, &self)
            .map_err(|e| ChemistryErrorKind::Composition(*self.span(), e))
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

impl<'a> ValidateInto<'a, Modifications> for ModificationsKdl {
    type Context = (&'a AtomicDatabase, &'a Residues);

    fn validate(self, ctx: Self::Context) -> ChemResult<Modifications> {
        // FIXME: This is a bit repetitive with the above...
        // FIXME: This is where I should build an index of functional groups, then pass that down to the validation
        // function for ModificationKdl
        self.modifications
            .into_iter()
            .map(|m| Ok((m.abbr.clone(), m.validate(ctx.0)?)))
            .collect()
    }
}

impl<'a> ValidateInto<'a, ModificationDescription> for ModificationKdl {
    type Context = &'a AtomicDatabase;

    fn validate(self, ctx: Self::Context) -> ChemResult<ModificationDescription> {
        Ok(ModificationDescription {
            name: self.name,
            lost: self.lost.validate(ctx)?,
            gained: self.gained.validate(ctx)?,
            targets: self.targets.into_iter().map(Target::from).collect(),
        })
    }
}

impl<'a> ValidateInto<'a, Bonds> for BondsKdl {
    type Context = &'a AtomicDatabase;

    fn validate(self, ctx: Self::Context) -> ChemResult<Bonds> {
        // FIXME: This is a bit repetitive with the above...
        self.bonds
            .into_iter()
            .map(|b| Ok((b.name.clone(), b.validate(ctx)?)))
            .collect()
    }
}

impl<'a> ValidateInto<'a, BondDescription> for BondKdl {
    type Context = &'a AtomicDatabase;

    fn validate(self, ctx: Self::Context) -> ChemResult<BondDescription> {
        Ok(BondDescription {
            from: self.from.into(),
            to: self.to.into(),
            lost: self.lost.validate(ctx)?,
        })
    }
}

impl From<TargetKdl> for Target {
    fn from(value: TargetKdl) -> Self {
        Self::new(value.functional_group, value.location, value.residue)
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

type ChemicalCompositionKdl = Spanned<String, Span>;

#[derive(Decode, Debug)]
#[knuffel(span_type=Span)]
struct BondKdl {
    #[knuffel(node_name)]
    name: String,
    #[knuffel(child)]
    from: TargetKdl,
    #[knuffel(child)]
    to: TargetKdl,
    #[knuffel(child, unwrap(argument))]
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
    #[knuffel(child, unwrap(argument))]
    lost: Option<ChemicalCompositionKdl>,
    #[knuffel(child, unwrap(argument))]
    gained: Option<ChemicalCompositionKdl>,
    #[knuffel(children(name = "targeting", non_empty))]
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
#[knuffel(span_type=Span)]
struct ResidueTypeKdl {
    #[knuffel(span)]
    span: Span,
    #[knuffel(node_name)]
    name: String,
    #[knuffel(children)]
    functional_groups: Vec<FunctionalGroupKdl>,
}

// FIXME: Be sure to write tests to check this errors for missing composition keys but works for nulled ones — also
// check that the `lost` and `gained` keys are the opposite — null isn't allowed, but they can be left out entirely
// NOTE: This forces composition to have a `null` value instead of being a totally optional key
type NullOr<T> = Option<T>;

#[derive(Decode, Debug)]
#[knuffel(span_type=Span)]
struct ResidueKdl {
    #[knuffel(span)]
    span: Span,
    #[knuffel(node_name)]
    residue_type: String,
    #[knuffel(argument)]
    abbr: String,
    #[knuffel(argument)]
    name: String,
    #[knuffel(child, unwrap(argument))]
    composition: NullOr<ChemicalCompositionKdl>,
    #[knuffel(children(name = "functional-group"))]
    functional_groups: Vec<FunctionalGroupKdl>,
}

#[derive(Clone, Decode, Debug)]
#[knuffel(span_type=Span)]
struct FunctionalGroupKdl {
    #[knuffel(span)]
    span: Span,
    #[knuffel(argument)]
    name: String,
    #[knuffel(property(name = "at"))]
    location: String,
}

#[cfg(test)]
mod tests {
    use indoc::indoc;
    use insta::{assert_debug_snapshot, assert_ron_snapshot};
    use miette::Result;
    use once_cell::sync::Lazy;

    use crate::{
        atomic_database::AtomicDatabase, polymer_chemistry::ValidateInto,
        testing_tools::assert_miette_snapshot,
    };

    use super::{
        ChemistryError, Modifications, ModificationsKdl, PolymerChemistry, PolymerChemistryKdl,
        Residues, ResiduesKdl,
    };

    static DB: Lazy<AtomicDatabase> = Lazy::new(|| {
        AtomicDatabase::from_kdl(
            "atomic_database.kdl",
            include_str!("../atomic_database.kdl"),
        )
        .unwrap()
    });

    const KDL: &str = include_str!("../muropeptide_chemistry.kdl");

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
            ".**.isotopes, .**.functional_groups" => insta::sorted_redaction()
        });
    }

    fn parse_residues(kdl: &str) -> Result<Residues, ChemistryError> {
        let residues: ResiduesKdl = knuffel::parse("test", kdl).unwrap();
        residues.validate(&DB).map_err(|e| e.finalize("test", kdl))
    }

    #[test]
    fn parse_residues_types_without_groups() {
        let kdl = indoc! {r#"
            types {
                AminoAcid
            }
            AminoAcid "G" "Glycine" {
                composition "C2H5NO2"
            }
        "#};
        let residues = parse_residues(kdl);
        assert!(residues.is_ok());
        assert_ron_snapshot!(residues.unwrap(), {
            ".**.isotopes" => insta::sorted_redaction()
        });
    }

    #[test]
    fn parse_residues_types_with_merged_groups() {
        let kdl = indoc! {r#"
            types {
                AminoAcid {
                    functional-group "Amino" at="N-terminal"
                    functional-group "Carboxyl" at="C-terminal"
                }
            }
            AminoAcid "K" "Lysine" {
                composition "C6H14N2O2"
                functional-group "Amino" at="Sidechain"
            }
        "#};
        let residues = parse_residues(kdl);
        assert!(residues.is_ok());
        assert_ron_snapshot!(residues.unwrap(), {
            ".**.isotopes, .**.functional_groups" => insta::sorted_redaction()
        });
    }

    #[test]
    fn parse_residues_without_type_section() {
        let kdl = indoc! {r#"
            AminoAcid "G" "Glycine" {
                composition "C2H5NO2"
            }
        "#};
        let residues: Result<ResiduesKdl, _> = knuffel::parse("test", kdl);
        assert_miette_snapshot!(residues);
    }

    #[test]
    fn parse_residues_without_types() {
        let kdl = indoc! {r#"
            types {
            //     AminoAcid
            }
            AminoAcid "G" "Glycine" {
                composition "C2H5NO2"
            }
        "#};
        let residues = parse_residues(kdl);
        assert_miette_snapshot!(residues);
    }

    #[test]
    fn parse_residues_with_duplicate_types() {
        let kdl = indoc! {r#"
            types {
                AminoAcid {
                    functional-group "Amino" at="N-terminal"
                    functional-group "Carboxyl" at="C-terminal"
                }
                Monosaccharide
                AminoAcid
            }
            AminoAcid "G" "Glycine" {
                composition "C2H5NO2"
            }
        "#};
        let residues = parse_residues(kdl);
        assert_miette_snapshot!(residues);
    }

    #[test]
    fn parse_residues_with_no_composition() {
        let kdl = indoc! {r#"
            types {
                AminoAcid
            }
            AminoAcid "K" "Lysine" {
                // composition "C6H14N2O2"
            }
        "#};
        let residues: Result<ResiduesKdl, _> = knuffel::parse("test", kdl);
        assert_miette_snapshot!(residues);
    }

    #[test]
    fn parse_residues_with_null_composition() {
        let kdl = indoc! {r#"
            types {
                AminoAcid
            }
            AminoAcid "K" "Lysine" {
                composition null
            }
        "#};
        let residues = parse_residues(kdl);
        assert!(residues.is_ok());
        assert_ron_snapshot!(residues.unwrap());
    }

    #[test]
    fn parse_residues_with_duplicate_compositions() {
        let kdl = indoc! {r#"
            types {
                AminoAcid
            }
            AminoAcid "K" "Lysine" {
                composition "C6H14N2O2"
                composition "C2H5NO2"
            }
        "#};
        let residues: Result<ResiduesKdl, _> = knuffel::parse("test", kdl);
        assert_miette_snapshot!(residues);
    }

    #[test]
    fn parse_residues_with_invalid_composition() {
        let kdl = indoc! {r#"
            types {
                AminoAcid
            }
            AminoAcid "K" "Lysine" {
                composition "C6H14[100Tc]N2O2"
            }
        "#};
        let residues = parse_residues(kdl);
        assert_miette_snapshot!(residues);
    }

    #[test]
    fn parse_residues_with_duplicate_functional_groups() {
        let kdl = indoc! {r#"
            types {
                AminoAcid
            }
            AminoAcid "K" "Lysine" {
                composition "C6H14N2O2"
                functional-group "Amino" at="N-terminal"
                functional-group "Carboxyl" at="C-terminal"
                functional-group "Amino" at="N-terminal"
            }
        "#};
        let residues = parse_residues(kdl);
        assert_miette_snapshot!(residues);
    }

    #[test]
    fn parse_residues_with_duplicate_merged_functional_groups() {
        let kdl = indoc! {r#"
            types {
                AminoAcid {
                    functional-group "Amino" at="N-terminal"
                }
            }
            AminoAcid "K" "Lysine" {
                composition "C6H14N2O2"
                functional-group "Amino" at="N-terminal"
            }
        "#};
        let residues = parse_residues(kdl);
        assert_miette_snapshot!(residues);
    }

    static RESIDUES: Lazy<Residues> = Lazy::new(|| {
        let residues_kdl = indoc! {r#"
            types {
                Monosaccharide {
                    functional-group "Hydroxyl" at="Reducing End"
                    functional-group "Hydroxyl" at="Nonreducing End"
                    functional-group "Hydroxyl" at="6-Position"
                }
                AminoAcid {
                    functional-group "Amino" at="N-Terminal"
                }
            }
            Monosaccharide "m" "N-Acetylmuramic Acid" {
                composition null
                functional-group "Carboxyl" at="Lactyl Ether"
            }
            AminoAcid "A" "Alanine" {
                composition null
            }
        "#};
        parse_residues(residues_kdl).unwrap()
    });

    fn parse_modifications(kdl: &str) -> Result<Modifications, ChemistryError> {
        let modifications: ModificationsKdl = knuffel::parse("test", kdl).unwrap();
        modifications
            .validate((&DB, &RESIDUES))
            .map_err(|e| e.finalize("test", kdl))
    }

    #[test]
    fn parse_empty_modifications() {
        let kdl = "";
        let modifications = parse_modifications(kdl);
        assert!(modifications.is_ok());
        assert_ron_snapshot!(modifications.unwrap());
    }

    #[test]
    fn parse_untargeted_modification() {
        let kdl = indoc! {r#"
            Ac "O-Acetylation" {
                // targeting "Hydroxyl" at="6-Position"
                // lost "H"
                // gained "C2H3O"
            }
        "#};
        let modifications: Result<ModificationsKdl, _> = knuffel::parse("test", kdl);
        assert_miette_snapshot!(modifications);
    }

    #[test]
    fn parse_massless_modification() {
        let kdl = indoc! {r#"
            Ac "O-Acetylation" {
                targeting "Hydroxyl" at="6-Position"
            }
        "#};
        let modifications = parse_modifications(kdl);
        assert!(modifications.is_ok());
        assert_ron_snapshot!(modifications.unwrap());
    }

    #[ignore]
    #[test]
    fn parse_duplicate_target_modification() {
        let kdl = indoc! {r#"
            Ac "O-Acetylation" {
                targeting "Hydroxyl" at="6-Position"
                targeting "Hydroxyl" at="6-Position"
            }
        "#};
        let modifications = parse_modifications(kdl);
        assert_miette_snapshot!(modifications);
    }
}
