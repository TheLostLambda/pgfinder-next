use std::{
    collections::HashMap,
    iter::{self, zip},
};

use knuffel::{
    span::{Span, Spanned},
    Decode,
};
use miette::{Diagnostic, LabeledSpan, NamedSource, Result};
use thiserror::Error;

use super::chemical_targets::{Target, TargetIndex};
use crate::{atomic_database::AtomicDatabase, ChemicalComposition, FunctionalGroup};

pub type Bonds = HashMap<String, BondDescription>;
pub type Modifications = HashMap<String, ModificationDescription>;
pub type Residues = HashMap<String, ResidueDescription>;

#[derive(Clone, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
pub struct PolymerChemistry {
    pub bonds: Bonds,
    pub modifications: Modifications,
    pub residues: Residues,
}

#[derive(Clone, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
pub struct BondDescription {
    from: Target,
    to: Target,
    lost: ChemicalComposition,
}

#[derive(Clone, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
pub struct ModificationDescription {
    name: String,
    lost: ChemicalComposition,
    gained: ChemicalComposition,
    targets: Vec<Target>,
}

#[derive(Clone, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
pub struct ResidueDescription {
    name: String,
    composition: ChemicalComposition,
    functional_groups: Vec<FunctionalGroup>,
}

// FIXME: Maybe move this `type` definition to `chemical_targets.rs`
type Targets<'a> = Vec<Target<&'a str>>;

impl<'a> From<&'a ResidueDescription> for Targets<'a> {
    fn from(value: &'a ResidueDescription) -> Self {
        value
            .functional_groups
            .iter()
            .map(|group| {
                Target::new(
                    group.name.as_str(),
                    Some(group.location.as_str()),
                    Some(value.name.as_str()),
                )
            })
            .collect()
    }
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
        let target_index: TargetIndex = residues.values().flat_map(Targets::from).collect();
        let ctx = (ctx, &target_index);
        Ok(PolymerChemistry {
            bonds: self.bonds.validate(ctx)?,
            modifications: self.modifications.validate(ctx)?,
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
    #[error("the target {2} overlaps with {} other target specifier{}", .1.len(), if .1.len() > 1 {"s"} else {""})]
    #[diagnostic(help("double-check for typos, or remove the overlapping target specifier"))]
    OverlappingTargets(Span, Vec<Span>, Target),
    #[error("the specifier {1} cannot target any currently defined residues")]
    #[diagnostic(help(
        "double-check for typos, or add new residues / groups that are targeted by this specifier"
    ))]
    NonexistentTarget(Span, Target),
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
            Self::OverlappingTargets(s1, ss, _) => iter::once((*s1, "this target"))
                .chain(zip(ss.clone(), iter::repeat("overlaps with")))
                .collect(),
            Self::NonexistentTarget(s, _) => vec![(*s, "targets nothing")],
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
    type Context = (&'a AtomicDatabase, &'a TargetIndex<'a>);

    fn validate(self, ctx: Self::Context) -> ChemResult<Modifications> {
        // FIXME: This is a bit repetitive with the above...
        // FIXME: This is where I should build an index of functional groups, then pass that down to the validation
        // function for ModificationKdl
        self.modifications
            .into_iter()
            .map(|m| Ok((m.abbr.clone(), m.validate(ctx)?)))
            .collect()
    }
}

impl<'a> ValidateInto<'a, ModificationDescription> for ModificationKdl {
    type Context = (&'a AtomicDatabase, &'a TargetIndex<'a>);

    fn validate(self, ctx: Self::Context) -> ChemResult<ModificationDescription> {
        let targets_and_spans: Vec<_> = self
            .targets
            .into_iter()
            .map(|t| {
                let span = t.span;
                Ok((t.validate(ctx.1)?, span))
            })
            // FIXME: That's also messy...
            .collect::<Result<_, _>>()?;

        // FIXME: Messy!
        let target_index: TargetIndex<Span> =
            targets_and_spans.iter().map(|(t, s)| (t, *s)).collect();
        for (target, span) in &targets_and_spans {
            let overlapping_targets: Vec<_> = target_index
                .get(target)
                .into_iter()
                .filter(|&s| s != span)
                .collect();
            if !overlapping_targets.is_empty() {
                return Err(ChemistryErrorKind::OverlappingTargets(
                    *span,
                    overlapping_targets.into_iter().copied().collect(),
                    target.clone(),
                ));
            }
        }

        // FIXME: Inline?
        let targets = targets_and_spans.into_iter().map(|(t, _)| t).collect();

        Ok(ModificationDescription {
            name: self.name,
            lost: self.lost.validate(ctx.0)?,
            gained: self.gained.validate(ctx.0)?,
            targets,
        })
    }
}

impl<'a> ValidateInto<'a, Bonds> for BondsKdl {
    type Context = (&'a AtomicDatabase, &'a TargetIndex<'a>);

    fn validate(self, ctx: Self::Context) -> ChemResult<Bonds> {
        // FIXME: This is a bit repetitive with the above...
        self.bonds
            .into_iter()
            .map(|b| Ok((b.name.clone(), b.validate(ctx)?)))
            .collect()
    }
}

impl<'a> ValidateInto<'a, BondDescription> for BondKdl {
    type Context = (&'a AtomicDatabase, &'a TargetIndex<'a>);

    fn validate(self, ctx: Self::Context) -> ChemResult<BondDescription> {
        Ok(BondDescription {
            from: self.from.validate(ctx.1)?,
            to: self.to.validate(ctx.1)?,
            lost: self.lost.validate(ctx.0)?,
        })
    }
}

impl<'a> ValidateInto<'a, Target> for TargetKdl {
    type Context = &'a TargetIndex<'a>;

    fn validate(self, ctx: Self::Context) -> ChemResult<Target> {
        let target = Target::new(self.functional_group, self.location, self.residue);
        // FIXME: Okay style?
        if ctx.get(&target).is_empty() {
            Err(ChemistryErrorKind::NonexistentTarget(self.span, target))
        } else {
            Ok(target)
        }
    }
}

// FIXME: Check that field names here line up with those in `lib.rs`!
#[derive(Debug, Decode)]
#[knuffel(span_type=Span)]
struct PolymerChemistryKdl {
    #[knuffel(child)]
    bonds: BondsKdl,
    #[knuffel(child)]
    modifications: ModificationsKdl,
    #[knuffel(child)]
    residues: ResiduesKdl,
}

#[derive(Debug, Decode)]
#[knuffel(span_type=Span)]
struct BondsKdl {
    #[knuffel(children)]
    bonds: Vec<BondKdl>,
}

type ChemicalCompositionKdl = Spanned<String, Span>;

#[derive(Debug, Decode)]
#[knuffel(span_type=Span)]
struct BondKdl {
    #[knuffel(node_name)]
    name: String,
    #[knuffel(child)]
    from: TargetKdl,
    #[knuffel(child)]
    to: TargetKdl,
    #[knuffel(child, unwrap(argument))]
    lost: Option<ChemicalCompositionKdl>,
}

#[derive(Debug, Decode)]
#[knuffel(span_type=Span)]
struct TargetKdl {
    #[knuffel(span)]
    span: Span,
    #[knuffel(argument)]
    functional_group: String,
    #[knuffel(property(name = "at"))]
    location: Option<String>,
    #[knuffel(property(name = "of"))]
    residue: Option<String>,
}

#[derive(Debug, Decode)]
#[knuffel(span_type=Span)]
struct ModificationsKdl {
    #[knuffel(children)]
    modifications: Vec<ModificationKdl>,
}

#[derive(Debug, Decode)]
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

#[derive(Debug, Decode)]
#[knuffel(span_type=Span)]
struct ResiduesKdl {
    #[knuffel(child, unwrap(children))]
    types: Vec<ResidueTypeKdl>,
    #[knuffel(children)]
    residues: Vec<ResidueKdl>,
}

#[derive(Debug, Decode)]
#[knuffel(span_type=Span)]
struct ResidueTypeKdl {
    #[knuffel(span)]
    span: Span,
    #[knuffel(node_name)]
    name: String,
    #[knuffel(children(name = "functional-group"))]
    functional_groups: Vec<FunctionalGroupKdl>,
}

// FIXME: Be sure to write tests to check this errors for missing composition keys but works for nulled ones — also
// check that the `lost` and `gained` keys are the opposite — null isn't allowed, but they can be left out entirely
// NOTE: This forces composition to have a `null` value instead of being a totally optional key
type NullOr<T> = Option<T>;

#[derive(Debug, Decode)]
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
    use miette::{Diagnostic, Report, Result};
    use once_cell::sync::Lazy;
    use thiserror::Error;

    use crate::{
        atomic_database::AtomicDatabase,
        polymerizer::{chemical_targets::TargetIndex, polymer_chemistry::ValidateInto},
        testing_tools::assert_miette_snapshot,
    };

    use super::{
        Bonds, BondsKdl, ChemistryError, Modifications, ModificationsKdl, PolymerChemistry,
        PolymerChemistryKdl, Residues, ResiduesKdl, Targets,
    };

    static DB: Lazy<AtomicDatabase> = Lazy::new(|| {
        AtomicDatabase::from_kdl(
            "atomic_database.kdl",
            include_str!("../../atomic_database.kdl"),
        )
        .unwrap()
    });

    const KDL: &str = include_str!("../../muropeptide_chemistry.kdl");

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

    static RESIDUE_INDEX: Lazy<TargetIndex> =
        Lazy::new(|| RESIDUES.values().flat_map(Targets::from).collect());

    fn parse_modifications(kdl: &str) -> Result<Modifications, ChemistryError> {
        let modifications: ModificationsKdl = knuffel::parse("test", kdl).unwrap();
        modifications
            .validate((&DB, &RESIDUE_INDEX))
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

    #[test]
    fn parse_chemically_invalid_lost_modification() {
        let kdl = indoc! {r#"
            Ac "O-Acetylation" {
                targeting "Hydroxyl" at="6-Position"
                lost "-H"
                gained "C2H3O"
            }
        "#};
        let modifications = parse_modifications(kdl);
        assert_miette_snapshot!(modifications);
    }

    #[test]
    fn parse_chemically_invalid_gained_modification() {
        let kdl = indoc! {r#"
            Ac "O-Acetylation" {
                targeting "Hydroxyl" at="6-Position"
                lost "H"
                gained "C2H-3O"
            }
        "#};
        let modifications = parse_modifications(kdl);
        assert_miette_snapshot!(modifications);
    }

    #[test]
    fn parse_null_lost_modification() {
        let kdl = indoc! {r#"
            Ac "O-Acetylation" {
                targeting "Hydroxyl" at="6-Position"
                lost null
                gained "C2H3O"
            }
        "#};
        let modifications: Result<ModificationsKdl, _> = knuffel::parse("test", kdl);
        assert_miette_snapshot!(modifications);
    }

    #[test]
    fn parse_null_gained_modification() {
        let kdl = indoc! {r#"
            Ac "O-Acetylation" {
                targeting "Hydroxyl" at="6-Position"
                lost "H"
                gained null
            }
        "#};
        let modifications: Result<ModificationsKdl, _> = knuffel::parse("test", kdl);
        assert_miette_snapshot!(modifications);
    }

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

    #[test]
    fn parse_modification_with_overlapping_broadening_targets() {
        let kdl = indoc! {r#"
            Ac "O-Acetylation" {
                targeting "Hydroxyl" at="6-Position" of="N-Acetylmuramic Acid"
                targeting "Hydroxyl" at="Reducing End"
                targeting "Hydroxyl"
            }
        "#};
        let modifications = parse_modifications(kdl);
        assert_miette_snapshot!(modifications);
    }

    #[test]
    fn parse_modification_with_overlapping_narrowing_targets() {
        let kdl = indoc! {r#"
            Ac "O-Acetylation" {
                targeting "Hydroxyl"
                targeting "Hydroxyl" at="Reducing End"
                targeting "Hydroxyl" at="6-Position" of="N-Acetylmuramic Acid"
            }
        "#};
        let modifications = parse_modifications(kdl);
        assert_miette_snapshot!(modifications);
    }

    #[test]
    fn parse_modification_with_overlapping_unordered_targets() {
        let kdl = indoc! {r#"
            Ac "O-Acetylation" {
                targeting "Hydroxyl" at="Reducing End"
                targeting "Hydroxyl"
                targeting "Hydroxyl" at="6-Position" of="N-Acetylmuramic Acid"
            }
        "#};
        let modifications = parse_modifications(kdl);
        assert_miette_snapshot!(modifications);
    }

    #[test]
    fn parse_modification_with_nonexistent_targets() {
        let kdl = indoc! {r#"
            Ac "O-Acetylation" {
                targeting "Hydroxyl" at="6-Position" of="N-Acetylmuramic Acid"
                targeting "Amino" at="Reducing End"
                targeting "Amino"
            }
        "#};
        let modifications = parse_modifications(kdl);
        assert_miette_snapshot!(modifications);
    }

    fn parse_bonds(kdl: &str) -> Result<Bonds, ChemistryError> {
        let bonds: BondsKdl = knuffel::parse("test", kdl).unwrap();
        bonds
            .validate((&DB, &RESIDUE_INDEX))
            .map_err(|e| e.finalize("test", kdl))
    }

    #[test]
    fn parse_muropeptide_bonds() {
        let kdl = indoc! {r#"
            Glycosidic {
                from "Hydroxyl" at="Reducing End"
                to "Hydroxyl" at="Nonreducing End"
            }
            Stem {
                from "Carboxyl" at="Lactyl Ether" of="N-Acetylmuramic Acid"
                to "Amino" at="N-Terminal"
                lost "H2O"
            }
        "#};
        let bonds = parse_bonds(kdl);
        assert!(bonds.is_ok());
        assert_ron_snapshot!(bonds.unwrap(), {
            "." => insta::sorted_redaction(),
            ".**.isotopes" => insta::sorted_redaction()
        });
    }

    #[test]
    fn parse_bond_with_chemically_invalid_lost() {
        let kdl = indoc! {r#"
            Glycosidic {
                from "Hydroxyl" at="Reducing End"
                to "Hydroxyl" at="Nonreducing End"
                lost "2H2O"
            }
        "#};
        let bonds = parse_bonds(kdl);
        assert_miette_snapshot!(bonds);
    }

    #[test]
    fn parse_bond_with_null_lost() {
        let kdl = indoc! {r#"
            Glycosidic {
                from "Hydroxyl" at="Reducing End"
                to "Hydroxyl" at="Nonreducing End"
                lost null
            }
        "#};
        let bonds: Result<BondsKdl, _> = knuffel::parse("test", kdl);
        assert_miette_snapshot!(bonds);
    }

    #[test]
    fn parse_bond_with_missing_to() {
        let kdl = indoc! {r#"
            Glycosidic {
                from "Hydroxyl" at="Reducing End"
            }
        "#};
        let bonds: Result<BondsKdl, _> = knuffel::parse("test", kdl);
        assert_miette_snapshot!(bonds);
    }

    #[test]
    fn parse_bond_with_missing_from() {
        let kdl = indoc! {r#"
            Glycosidic {
                to "Hydroxyl" at="Nonreducing End"
            }
        "#};
        let bonds: Result<BondsKdl, _> = knuffel::parse("test", kdl);
        assert_miette_snapshot!(bonds);
    }

    #[test]
    fn parse_bond_with_nonexistent_from() {
        let kdl = indoc! {r#"
            Glycosidic {
                from "Hydroxyl" at="Sidechain"
                to "Hydroxyl" at="Nonreducing End"
            }
        "#};
        let bonds = parse_bonds(kdl);
        assert_miette_snapshot!(bonds);
    }

    #[test]
    fn parse_bond_with_nonexistent_to() {
        let kdl = indoc! {r#"
            Glycosidic {
                from "Hydroxyl" at="Reducing End"
                to "Hydroxyl" at="Nonreducing End" of="Alanine"
            }
        "#};
        let bonds = parse_bonds(kdl);
        assert_miette_snapshot!(bonds);
    }

    #[test]
    fn parse_bond_with_duplicate_from() {
        let kdl = indoc! {r#"
            Glycosidic {
                from "Hydroxyl" at="Reducing End"
                from "Hydroxyl" at="Nonreducing End"
            }
        "#};
        let bonds: Result<BondsKdl, _> = knuffel::parse("test", kdl);
        assert_miette_snapshot!(bonds);
    }

    #[test]
    fn parse_bond_with_duplicate_to() {
        let kdl = indoc! {r#"
            Glycosidic {
                to "Hydroxyl" at="Reducing End"
                to "Hydroxyl" at="Nonreducing End"
            }
        "#};
        let bonds: Result<BondsKdl, _> = knuffel::parse("test", kdl);
        assert_miette_snapshot!(bonds);
    }

    #[derive(Debug, Error, Diagnostic)]
    #[error("encountered an error during testing")]
    #[diagnostic(transparent)]
    struct WrapErr(Report);

    #[test]
    fn parse_complete_parse_error() {
        let kdl = indoc! {r#"
            bonds {
                Peptide {
                    from "Carboxyl" at="C-Terminal"
                    to "Amino" at="N-Terminal"
                    lost "H2O"
                }
            }

            modifications {
                Am "Amidation" {
                    targeting "Carboxyl" at="Sidechain"
                    lost "OH"
                    gained "NH2"
                }
            }

            residues {
                types {
                    AminoAcid {
                        group "Amino" at="N-Terminal"
                        functional-group "Carboxyl" at="C-Terminal"
                    }
                }
                AminoAcid "E" "Glutamic Acid" {
                    composition "C5H9NO4"
                    functional-group "Carboxyl" at="Sidechain"
                }
            }
        "#};
        let polymer_chemistry = PolymerChemistry::from_kdl(&DB, "test", kdl);
        assert_miette_snapshot!(polymer_chemistry.map_err(WrapErr));
    }

    #[test]
    fn parse_complete_bonds_error() {
        let kdl = indoc! {r#"
            bonds {
                Peptide {
                    from "Carboxyl" at="C-Terminal"
                    to "Amino" at="C-Terminal"
                    lost "H2O"
                }
            }

            modifications {
                Am "Amidation" {
                    targeting "Carboxyl" at="Sidechain"
                    lost "OH"
                    gained "NH2"
                }
            }

            residues {
                types {
                    AminoAcid {
                        functional-group "Amino" at="N-Terminal"
                        functional-group "Carboxyl" at="C-Terminal"
                    }
                }
                AminoAcid "E" "Glutamic Acid" {
                    composition "C5H9NO4"
                    functional-group "Carboxyl" at="Sidechain"
                }
            }
        "#};
        let polymer_chemistry = PolymerChemistry::from_kdl(&DB, "test", kdl);
        assert_miette_snapshot!(polymer_chemistry.map_err(WrapErr));
    }

    #[test]
    fn parse_complete_modifications_error() {
        let kdl = indoc! {r#"
            bonds {
                Peptide {
                    from "Carboxyl" at="C-Terminal"
                    to "Amino" at="N-Terminal"
                    lost "H2O"
                }
            }

            modifications {
                Am "Amidation" {
                    targeting "Carboxyl" at="Sidechain"
                    targeting "Carboxyl"
                    lost "OH"
                    gained "NH2"
                }
            }

            residues {
                types {
                    AminoAcid {
                        functional-group "Amino" at="N-Terminal"
                        functional-group "Carboxyl" at="C-Terminal"
                    }
                }
                AminoAcid "E" "Glutamic Acid" {
                    composition "C5H9NO4"
                    functional-group "Carboxyl" at="Sidechain"
                }
            }
        "#};
        let polymer_chemistry = PolymerChemistry::from_kdl(&DB, "test", kdl);
        assert_miette_snapshot!(polymer_chemistry.map_err(WrapErr));
    }

    #[test]
    fn parse_complete_residues_error() {
        let kdl = indoc! {r#"
            bonds {
                Peptide {
                    from "Carboxyl" at="C-Terminal"
                    to "Amino" at="N-Terminal"
                    lost "H2O"
                }
            }

            modifications {
                Am "Amidation" {
                    targeting "Carboxyl" at="Sidechain"
                    lost "OH"
                    gained "NH2"
                }
            }

            residues {
                types {
                    AminoAcid {
                        functional-group "Amino" at="N-Terminal"
                        functional-group "Carboxyl" at="C-Terminal"
                    }
                }
                SuperAminoAcid "E" "Glutamic Acid" {
                    composition "C5H9NO4"
                    functional-group "Carboxyl" at="Sidechain"
                }
            }
        "#};
        let polymer_chemistry = PolymerChemistry::from_kdl(&DB, "test", kdl);
        assert_miette_snapshot!(polymer_chemistry.map_err(WrapErr));
    }
}
