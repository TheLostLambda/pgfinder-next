// Standard Library Imports
use std::{
    collections::hash_map::Entry,
    iter::{self, zip},
};

// External Crate Imports
use ahash::{HashMap, HashMapExt};
use knuffel::{
    span::{Span, Spanned},
    Decode,
};
use miette::{Diagnostic, LabeledSpan, NamedSource, Result};
use thiserror::Error;

// Local Crate Imports
use super::target::{Index, Target};
use crate::{atoms::atomic_database::AtomicDatabase, errors::PolychemError, ChemicalComposition};

// Public API ==========================================================================================================

#[derive(Clone, Eq, PartialEq, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
pub struct PolymerDatabase<'a> {
    pub bonds: Bonds<'a>,
    pub modifications: Modifications<'a>,
    pub residues: Residues<'a>,
}

impl<'a> PolymerDatabase<'a> {
    pub fn new(
        atomic_db: &'a AtomicDatabase,
        file_name: impl AsRef<str>,
        kdl_text: impl AsRef<str>,
    ) -> Result<Self> {
        let parsed_db: PolymerDatabaseKdl = knuffel::parse(file_name.as_ref(), kdl_text.as_ref())?;
        parsed_db
            .validate(atomic_db)
            .map_err(|e| e.finalize(file_name, kdl_text).into())
    }
}

// Private Types =======================================================================================================

type Bonds<'a> = HashMap<String, BondDescription<'a>>;
type Modifications<'a> = HashMap<String, ModificationDescription<'a>>;
type Residues<'a> = HashMap<String, ResidueDescription<'a>>;

// ---------------------------------------------------------------------------------------------------------------------

#[derive(Clone, Eq, PartialEq, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
pub struct BondDescription<'a> {
    pub name: String,
    pub lost: ChemicalComposition<'a>,
    pub gained: ChemicalComposition<'a>,
    pub from: Target,
    pub to: Target,
}

#[derive(Clone, Eq, PartialEq, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
pub struct ModificationDescription<'a> {
    pub name: String,
    pub lost: ChemicalComposition<'a>,
    pub gained: ChemicalComposition<'a>,
    pub targets: Vec<Target>,
}

#[derive(Clone, Eq, PartialEq, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
pub struct ResidueDescription<'a> {
    pub name: String,
    pub composition: ChemicalComposition<'a>,
    pub functional_groups: Vec<FunctionalGroupDescription>,
}

#[derive(Clone, Eq, PartialEq, Hash, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
pub struct FunctionalGroupDescription {
    pub name: String,
    pub location: String,
}

// KDL File Schema =====================================================================================================

#[derive(Debug, Decode)]
#[knuffel(span_type=Span)]
struct PolymerDatabaseKdl {
    #[knuffel(child)]
    bonds: BondsKdl,
    #[knuffel(child)]
    modifications: ModificationsKdl,
    #[knuffel(child)]
    residues: ResiduesKdl,
}

// ---------------------------------------------------------------------------------------------------------------------

#[derive(Debug, Decode)]
#[knuffel(span_type=Span)]
struct BondsKdl {
    #[knuffel(children)]
    bonds: Vec<BondKdl>,
}

#[derive(Debug, Decode)]
#[knuffel(span_type=Span)]
struct ModificationsKdl {
    #[knuffel(children)]
    modifications: Vec<ModificationKdl>,
}

#[derive(Debug, Decode)]
#[knuffel(span_type=Span)]
struct ResiduesKdl {
    #[knuffel(child, unwrap(children))]
    types: Vec<ResidueTypeKdl>,
    #[knuffel(children)]
    residues: Vec<ResidueKdl>,
}

// ---------------------------------------------------------------------------------------------------------------------

#[derive(Debug, Decode)]
#[knuffel(span_type=Span)]
struct BondKdl {
    #[knuffel(node_name)]
    abbr: String,
    #[knuffel(argument)]
    name: String,
    #[knuffel(child)]
    from: TargetKdl,
    #[knuffel(child)]
    to: TargetKdl,
    #[knuffel(child, unwrap(argument))]
    lost: Option<ChemicalCompositionKdl>,
    #[knuffel(child, unwrap(argument))]
    gained: Option<ChemicalCompositionKdl>,
}

#[derive(Debug, Decode)]
#[knuffel(span_type=Span)]
struct ModificationKdl {
    #[knuffel(node_name)]
    abbr: String,
    #[knuffel(argument)]
    name: String,
    #[knuffel(children(name = "targeting", non_empty))]
    targets: Vec<TargetKdl>,
    #[knuffel(child, unwrap(argument))]
    lost: Option<ChemicalCompositionKdl>,
    #[knuffel(child, unwrap(argument))]
    gained: Option<ChemicalCompositionKdl>,
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

// ---------------------------------------------------------------------------------------------------------------------

#[derive(Debug, Decode)]
#[knuffel(span_type=Span)]
struct TargetKdl {
    #[knuffel(span)]
    span: Span,
    #[knuffel(argument)]
    group: String,
    #[knuffel(property(name = "at"))]
    location: Option<String>,
    #[knuffel(property(name = "of"))]
    residue: Option<String>,
}

type ChemicalCompositionKdl = Spanned<String, Span>;

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

// NOTE: This forces composition to have a `null` value instead of being an entirely optional node
type NullOr<T> = Option<T>;

// Contextual Validation Trait  ========================================================================================

type ChemResult<T> = Result<T, ChemistryErrorKind>;

trait ValidateInto<'c, T> {
    type Context: 'c;

    fn validate(self, ctx: Self::Context) -> ChemResult<T>;
}

// Polymer Database Validation =========================================================================================

impl<'a> ValidateInto<'a, PolymerDatabase<'a>> for PolymerDatabaseKdl {
    type Context = &'a AtomicDatabase;

    fn validate(self, ctx: Self::Context) -> ChemResult<PolymerDatabase<'a>> {
        let residues = self.residues.validate(ctx)?;
        let target_index: Index = residues.values().flat_map(Targets::from).collect();
        let ctx = (ctx, &target_index);
        Ok(PolymerDatabase {
            bonds: self.bonds.validate(ctx)?,
            modifications: self.modifications.validate(ctx)?,
            residues,
        })
    }
}

// Validate Residue Types and Residues =================================================================================

impl<'a> ValidateInto<'a, Residues<'a>> for ResiduesKdl {
    type Context = &'a AtomicDatabase;

    fn validate(self, ctx: Self::Context) -> ChemResult<Residues<'a>> {
        let types = self.types.validate(())?;
        self.residues
            .into_iter()
            .map(|r| r.validate((ctx, &types)))
            .collect()
        // TODO: Report an error if any residues are defined twice / have the same abbr!
    }
}

// ---------------------------------------------------------------------------------------------------------------------

type ResidueTypes = HashMap<String, Vec<FunctionalGroupKdl>>;

// NOTE: The Context = () means this is essentially a TryInto implementation, but Rust's orphan rules mean I can't
// actually implement TryInto, because Vec<_> is not a local type (despite its parameter, ResidueTypeKdl being one)
impl ValidateInto<'_, ResidueTypes> for Vec<ResidueTypeKdl> {
    type Context = ();

    fn validate(self, _ctx: Self::Context) -> ChemResult<ResidueTypes> {
        let mut seen_types = HashMap::new();

        for residue_type in self {
            match seen_types.entry(residue_type.name) {
                Entry::Occupied(e) => {
                    let (residue_name, (first_defined_at, _)) = e.remove_entry();
                    return Err(ChemistryErrorKind::DuplicateResidueType(
                        first_defined_at,
                        residue_type.span,
                        residue_name,
                    ));
                }
                Entry::Vacant(e) => e.insert((residue_type.span, residue_type.functional_groups)),
            };
        }

        Ok(seen_types.into_iter().map(|(k, (_, v))| (k, v)).collect())
    }
}

type ResidueEntry<'a> = (String, ResidueDescription<'a>);

impl<'a: 'r, 'r> ValidateInto<'r, ResidueEntry<'a>> for ResidueKdl {
    type Context = (&'a AtomicDatabase, &'r ResidueTypes);

    fn validate(self, ctx: Self::Context) -> ChemResult<ResidueEntry<'a>> {
        let groups_from_type = ctx
            .1
            .get(&self.residue_type)
            .ok_or_else(|| ChemistryErrorKind::UndefinedResidueType(self.span, self.residue_type))?
            .clone();

        let mut seen_groups = HashMap::new();

        for functional_group_kdl in groups_from_type.into_iter().chain(self.functional_groups) {
            let (functional_group, group_span) = functional_group_kdl.into();

            match seen_groups.entry(functional_group) {
                Entry::Occupied(e) => {
                    let (functional_group, first_defined_at) = e.remove_entry();
                    return Err(ChemistryErrorKind::DuplicateFunctionalGroup(
                        first_defined_at,
                        group_span,
                        functional_group.name,
                        functional_group.location,
                    ));
                }
                Entry::Vacant(e) => e.insert(group_span),
            };
        }

        Ok((
            self.abbr,
            ResidueDescription {
                name: self.name,
                composition: self.composition.validate(ctx.0)?,
                functional_groups: seen_groups.into_keys().collect(),
            },
        ))
    }
}

// ---------------------------------------------------------------------------------------------------------------------

impl<'a> ValidateInto<'a, ChemicalComposition<'a>> for Option<ChemicalCompositionKdl> {
    type Context = &'a AtomicDatabase;

    fn validate(self, ctx: Self::Context) -> ChemResult<ChemicalComposition<'a>> {
        self.map_or_else(|| Ok(ChemicalComposition::default()), |c| c.validate(ctx))
    }
}

// ---------------------------------------------------------------------------------------------------------------------

impl<'a> ValidateInto<'a, ChemicalComposition<'a>> for ChemicalCompositionKdl {
    type Context = &'a AtomicDatabase;

    fn validate(self, ctx: Self::Context) -> ChemResult<ChemicalComposition<'a>> {
        ChemicalComposition::new(ctx, &self)
            .map_err(|e| ChemistryErrorKind::Composition(*self.span(), *e))
    }
}

// Validate Bonds ======================================================================================================

impl<'a: 't, 't> ValidateInto<'t, Bonds<'a>> for BondsKdl {
    type Context = (&'a AtomicDatabase, &'t Index<'t>);

    fn validate(self, ctx: Self::Context) -> ChemResult<Bonds<'a>> {
        self.bonds.into_iter().map(|b| b.validate(ctx)).collect()
        // TODO: Report an error if any bonds are defined twice / have the same abbr!
    }
}

// ---------------------------------------------------------------------------------------------------------------------

type BondEntry<'a> = (String, BondDescription<'a>);

impl<'a: 't, 't> ValidateInto<'t, BondEntry<'a>> for BondKdl {
    type Context = (&'a AtomicDatabase, &'t Index<'t>);

    fn validate(self, ctx: Self::Context) -> ChemResult<BondEntry<'a>> {
        Ok((
            self.abbr,
            BondDescription {
                name: self.name,
                from: self.from.validate(ctx.1)?,
                to: self.to.validate(ctx.1)?,
                lost: self.lost.validate(ctx.0)?,
                gained: self.gained.validate(ctx.0)?,
            },
        ))
    }
}

// ---------------------------------------------------------------------------------------------------------------------

impl<'t> ValidateInto<'t, Target> for TargetKdl {
    type Context = &'t Index<'t>;

    fn validate(self, ctx: Self::Context) -> ChemResult<Target> {
        let target = Target::new(self.group, self.location, self.residue);

        if ctx.contains_target(&target) {
            Ok(target)
        } else {
            Err(ChemistryErrorKind::NonexistentTarget(self.span, target))
        }
    }
}

// Validate Modifications ==============================================================================================

impl<'a: 't, 't> ValidateInto<'t, Modifications<'a>> for ModificationsKdl {
    type Context = (&'a AtomicDatabase, &'t Index<'t>);

    fn validate(self, ctx: Self::Context) -> ChemResult<Modifications<'a>> {
        self.modifications
            .into_iter()
            .map(|m| m.validate(ctx))
            .collect()
        // TODO: Report an error if any modifications are defined twice / have the same abbr!
    }
}

// ---------------------------------------------------------------------------------------------------------------------

type ModificationEntry<'a> = (String, ModificationDescription<'a>);

impl<'a: 't, 't> ValidateInto<'t, ModificationEntry<'a>> for ModificationKdl {
    type Context = (&'a AtomicDatabase, &'t Index<'t>);

    fn validate(self, ctx: Self::Context) -> ChemResult<ModificationEntry<'a>> {
        let targets_and_spans: Vec<_> = self
            .targets
            .into_iter()
            .map(|t| t.validate(ctx.1))
            .collect::<Result<_, _>>()?;

        let target_index: Index<_> = targets_and_spans.iter().map(|(t, s)| (t, *s)).collect();
        for (target, span) in &targets_and_spans {
            let overlapping_targets: Vec<_> = target_index
                .matches(target, |_, _, _, s| s)
                .copied()
                .filter(|s| s != span)
                .collect();

            if !overlapping_targets.is_empty() {
                return Err(ChemistryErrorKind::OverlappingTargets(
                    *span,
                    overlapping_targets,
                    target.clone(),
                ));
            }
        }

        Ok((
            self.abbr,
            ModificationDescription {
                name: self.name,
                lost: self.lost.validate(ctx.0)?,
                gained: self.gained.validate(ctx.0)?,
                targets: targets_and_spans.into_iter().map(|(t, _)| t).collect(),
            },
        ))
    }
}

// ---------------------------------------------------------------------------------------------------------------------

type TargetEntry = (Target, Span);

impl<'t> ValidateInto<'t, TargetEntry> for TargetKdl {
    type Context = &'t Index<'t>;

    fn validate(self, ctx: Self::Context) -> ChemResult<TargetEntry> {
        let span = self.span;
        Ok((self.validate(ctx)?, span))
    }
}

// Infallible Conversions =============================================================================================

type Targets<'a> = Vec<Target<&'a str>>;

impl<'a> From<&'a ResidueDescription<'a>> for Targets<'a> {
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

type FunctionalGroupEntry = (FunctionalGroupDescription, Span);

impl From<FunctionalGroupKdl> for FunctionalGroupEntry {
    fn from(value: FunctionalGroupKdl) -> Self {
        (
            FunctionalGroupDescription {
                name: value.name,
                location: value.location,
            },
            value.span,
        )
    }
}

// Validation Error Types and Trait Implementations  ===================================================================

#[derive(Debug, Error)]
#[error("failed to validate polymer database file")]
struct ChemistryError {
    kdl: NamedSource<String>,
    #[source]
    kind: ChemistryErrorKind,
}

// NOTE: This is manually implemented because the list of labels is dynamic and needs to be extracted from `self.kind`
impl Diagnostic for ChemistryError {
    fn source_code(&self) -> Option<&dyn miette::SourceCode> {
        Some(&self.kdl)
    }

    fn labels(&self) -> Option<Box<dyn Iterator<Item = miette::LabeledSpan> + '_>> {
        Some(Box::new(self.kind.labels().into_iter().map(|(s, l)| {
            LabeledSpan::new_with_span(Some(l.to_owned()), *s)
        })))
    }

    fn diagnostic_source(&self) -> Option<&dyn Diagnostic> {
        Some(&self.kind)
    }
}

#[derive(Clone, Debug, Diagnostic, Error)]
enum ChemistryErrorKind {
    #[error("the residue type {2:?} has already been defined")]
    #[diagnostic(help("consider consolidating duplicate types or picking a new type name"))]
    DuplicateResidueType(Span, Span, String),

    #[error("the residue type {1:?} is undefined")]
    #[diagnostic(help("double-check for typos, or add {1:?} to the types section"))]
    UndefinedResidueType(Span, String),

    #[error("the functional group {2:?} has already been defined at {3:?}")]
    #[diagnostic(help("double-check for typos, or remove the duplicate functional group"))]
    DuplicateFunctionalGroup(Span, Span, String, String),

    #[error("the specifier {1} cannot target any currently defined residues")]
    #[diagnostic(help(
        "double-check for typos, or add new residues / groups that are targeted by this specifier"
    ))]
    NonexistentTarget(Span, Target),

    #[error("the target {2} overlaps with {} other target specifier{}", .1.len(), if .1.len() > 1 {"s"} else {""})]
    #[diagnostic(help("double-check for typos, or remove the overlapping target specifier"))]
    OverlappingTargets(Span, Vec<Span>, Target),

    #[error("polymer database file contained an invalid chemical composition")]
    Composition(
        Span,
        #[source]
        #[diagnostic_source]
        PolychemError,
    ),
}

impl ChemistryErrorKind {
    fn labels(&self) -> Vec<(&Span, &'static str)> {
        match self {
            Self::DuplicateResidueType(s1, s2, _)
            | Self::DuplicateFunctionalGroup(s1, s2, _, _) => {
                vec![(s1, "first defined here"), (s2, "then again here")]
            }
            Self::UndefinedResidueType(s, _) => vec![(s, "undefined residue type")],
            Self::NonexistentTarget(s, _) => vec![(s, "targets nothing")],
            Self::OverlappingTargets(s1, ss, _) => iter::once((s1, "this target"))
                .chain(zip(ss, iter::repeat("overlaps with")))
                .collect(),
            Self::Composition(s, _) => vec![(s, "invalid chemical composition")],
        }
    }

    fn finalize(self, file_name: impl AsRef<str>, kdl: impl AsRef<str>) -> ChemistryError {
        let kdl = NamedSource::new(file_name, kdl.as_ref().to_owned());
        ChemistryError { kdl, kind: self }
    }
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use indoc::indoc;
    use insta::{assert_debug_snapshot, assert_ron_snapshot, with_settings};
    use miette::{Diagnostic, Report};
    use once_cell::sync::Lazy;
    use thiserror::Error;

    use crate::testing_tools::assert_miette_snapshot;

    use super::*;

    static DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);

    const KDL: &str = include_str!("../../tests/data/polymer_database.kdl");

    #[test]
    fn parse_muropeptide_chemistry() {
        let db: PolymerDatabaseKdl = knuffel::parse("test_polymer_database.kdl", KDL).unwrap();
        with_settings!({filters => vec![
            (r"Span\([^)]*\)", "<SPAN>"),
        ]}, {
            assert_debug_snapshot!(db);
        });
    }

    #[test]
    fn build_muropeptide_chemistry() {
        let db = PolymerDatabase::new(&DB, "test_polymer_database.kdl", KDL).unwrap();
        assert_ron_snapshot!(db, {
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

    static RESIDUE_INDEX: Lazy<Index> =
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
            Gly "Glycosidic" {
                from "Hydroxyl" at="Reducing End"
                to "Hydroxyl" at="Nonreducing End"
            }
            Stem "MurNAc -> Stem Peptide" {
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
            Gly "Glycosidic" {
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
            Gly "Glycosidic" {
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
            Gly "Glycosidic" {
                from "Hydroxyl" at="Reducing End"
            }
        "#};
        let bonds: Result<BondsKdl, _> = knuffel::parse("test", kdl);
        assert_miette_snapshot!(bonds);
    }

    #[test]
    fn parse_bond_with_missing_from() {
        let kdl = indoc! {r#"
            Gly "Glycosidic" {
                to "Hydroxyl" at="Nonreducing End"
            }
        "#};
        let bonds: Result<BondsKdl, _> = knuffel::parse("test", kdl);
        assert_miette_snapshot!(bonds);
    }

    #[test]
    fn parse_bond_with_nonexistent_from() {
        let kdl = indoc! {r#"
            Gly "Glycosidic" {
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
            Gly "Glycosidic" {
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
            Gly "Glycosidic" {
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
            Gly "Glycosidic" {
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
                Pep "Peptide" {
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
        let db = PolymerDatabase::new(&DB, "test", kdl);
        assert_miette_snapshot!(db.map_err(WrapErr));
    }

    #[test]
    fn parse_complete_bonds_error() {
        let kdl = indoc! {r#"
            bonds {
                Pep "Peptide" {
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
        let db = PolymerDatabase::new(&DB, "test", kdl);
        assert_miette_snapshot!(db.map_err(WrapErr));
    }

    #[test]
    fn parse_complete_modifications_error() {
        let kdl = indoc! {r#"
            bonds {
                Pep "Peptide" {
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
        let db = PolymerDatabase::new(&DB, "test", kdl);
        assert_miette_snapshot!(db.map_err(WrapErr));
    }

    #[test]
    fn parse_complete_residues_error() {
        let kdl = indoc! {r#"
            bonds {
                Pep "Peptide" {
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
        let db = PolymerDatabase::new(&DB, "test", kdl);
        assert_miette_snapshot!(db.map_err(WrapErr));
    }
}
