use std::{collections::HashMap, convert::Infallible};

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

impl PolymerChemistry {
    pub fn from_kdl(
        db: &AtomicDatabase,
        file_name: impl AsRef<str>,
        text: impl AsRef<str>,
    ) -> Result<Self> {
        let parsed_chemistry: PolymerChemistryKdl =
            knuffel::parse(file_name.as_ref(), text.as_ref())?;
        parsed_chemistry.ctx_try_into(db).map_err(|(span, kind)| {
            let kdl = NamedSource::new(file_name, text.as_ref().to_string());
            ChemistryError { kdl, span, kind }.into()
        })
    }
}

impl<'a> ContextualTryFrom<'a, PolymerChemistryKdl> for PolymerChemistry {
    type Error = ProtoChemistryError;
    type Context = &'a AtomicDatabase;

    fn ctx_try_from(ctx: Self::Context, value: PolymerChemistryKdl) -> Result<Self, Self::Error> {
        Ok(Self {
            bonds: value.bonds.ctx_try_into(ctx)?,
            modifications: value.modifications.ctx_try_into(ctx)?,
            residues: value.residues.ctx_try_into(ctx)?,
        })
    }
}

impl<'a> ContextualTryFrom<'a, ResiduesKdl> for HashMap<String, ResidueDescription> {
    type Error = ProtoChemistryError;
    type Context = &'a AtomicDatabase;

    fn ctx_try_from(ctx: Self::Context, value: ResiduesKdl) -> Result<Self, Self::Error> {
        let types: HashMap<_, _> = value
            .types
            .into_iter()
            .map(ResidueTypeEntry::from)
            .collect();
        let to_residue = |r: ResidueKdl| {
            let mut functional_groups: Vec<_> = r
                .functional_groups
                .into_iter()
                .map(FunctionalGroup::from)
                .collect();
            functional_groups.extend(types[&r.residue_type].clone());
            Ok(ResidueDescription {
                name: r.name,
                composition: r.composition.ctx_try_into(ctx)?,
                functional_groups,
            })
        };
        value
            .residues
            .into_iter()
            // FIXME: Big hmm for the ? operator here...
            .map(|r| Ok((r.abbr.clone(), to_residue(r)?)))
            .collect()
    }
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

type ProtoChemistryError = (SourceSpan, ChemistryErrorKind);

// FIXME: This should definitely live in it's own crate!

trait ContextualTryInto<'a, T>: Sized {
    /// The type returned in the event of a conversion error.
    type Error;
    /// The type of context needed for this conversion
    type Context: 'a;

    // FIXME: Wow, I really hate this ctx naming...
    /// Performs the conversion.
    fn ctx_try_into(self, ctx: Self::Context) -> Result<T, Self::Error>;
}

trait ContextualTryFrom<'a, T>: Sized {
    /// The type returned in the event of a conversion error.
    type Error;
    /// The type of context needed for this conversion
    type Context: 'a;

    /// Performs the conversion.
    fn ctx_try_from(ctx: Self::Context, value: T) -> Result<Self, Self::Error>;
}

// TryFrom implies TryInto
impl<'a, T, U> ContextualTryInto<'a, U> for T
where
    U: ContextualTryFrom<'a, T>,
{
    type Error = U::Error;
    type Context = U::Context;

    #[inline]
    fn ctx_try_into(self, ctx: U::Context) -> Result<U, U::Error> {
        U::ctx_try_from(ctx, self)
    }
}

// FIXME: Check that ContextualTryFrom<T> for T is implemented by this?
// Infallible conversions are semantically equivalent to fallible conversions
// with an uninhabited error type.
impl<T, U> ContextualTryFrom<'static, U> for T
where
    U: Into<T>,
{
    type Error = Infallible;
    type Context = ();

    #[inline]
    fn ctx_try_from(_ctx: Self::Context, value: U) -> Result<Self, Self::Error> {
        Ok(U::into(value))
    }
}

impl<'a> ContextualTryFrom<'a, ChemicalCompositionKdl> for ChemicalComposition {
    type Error = ProtoChemistryError;
    type Context = &'a AtomicDatabase;

    fn ctx_try_from(
        ctx: Self::Context,
        value: ChemicalCompositionKdl,
    ) -> Result<Self, Self::Error> {
        Self::new(ctx, &*value).map_err(|e| (value.span().clone().into(), e.into()))
    }
}

impl<'a> ContextualTryFrom<'a, Option<ChemicalCompositionKdl>> for ChemicalComposition {
    type Error = ProtoChemistryError;
    type Context = &'a AtomicDatabase;

    fn ctx_try_from(
        ctx: Self::Context,
        value: Option<ChemicalCompositionKdl>,
    ) -> Result<Self, Self::Error> {
        // FIXME: Questionable style...
        Ok(if let Some(s) = value {
            s.ctx_try_into(ctx)?
        } else {
            Self::default()
        })
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

impl<'a> ContextualTryFrom<'a, ModificationsKdl> for HashMap<String, ModificationDescription> {
    type Error = ProtoChemistryError;
    type Context = &'a AtomicDatabase;

    fn ctx_try_from(ctx: Self::Context, value: ModificationsKdl) -> Result<Self, Self::Error> {
        // FIXME: This is a bit repetitive with the above...
        value
            .modifications
            .into_iter()
            .map(|m| Ok((m.abbr.clone(), m.ctx_try_into(ctx)?)))
            .collect()
    }
}

impl<'a> ContextualTryFrom<'a, ModificationKdl> for ModificationDescription {
    type Error = ProtoChemistryError;
    type Context = &'a AtomicDatabase;

    fn ctx_try_from(ctx: Self::Context, value: ModificationKdl) -> Result<Self, Self::Error> {
        Ok(Self {
            name: value.name,
            lost: value.lost.ctx_try_into(ctx)?,
            gained: value.gained.ctx_try_into(ctx)?,
            targets: value.targets.into_iter().map(Target::from).collect(),
        })
    }
}

impl<'a> ContextualTryFrom<'a, BondsKdl> for HashMap<String, BondDescription> {
    type Error = ProtoChemistryError;
    type Context = &'a AtomicDatabase;

    fn ctx_try_from(ctx: Self::Context, value: BondsKdl) -> Result<Self, Self::Error> {
        // FIXME: This is a bit repetitive with the above...
        value
            .bonds
            .into_iter()
            .map(|b| Ok((b.name.clone(), b.ctx_try_into(ctx)?)))
            .collect()
    }
}

impl<'a> ContextualTryFrom<'a, BondKdl> for BondDescription {
    type Error = ProtoChemistryError;
    type Context = &'a AtomicDatabase;

    fn ctx_try_from(ctx: Self::Context, value: BondKdl) -> Result<Self, Self::Error> {
        Ok(Self {
            from: value.from.into(),
            to: value.to.into(),
            lost: value.lost.ctx_try_into(ctx)?,
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
type NullOr<T> = Option<T>;

#[derive(Decode, Debug)]
#[knuffel(span_type=Span)]
struct ResidueKdl {
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

#[derive(Decode, Debug)]
struct FunctionalGroupKdl {
    #[knuffel(argument)]
    name: String,
    #[knuffel(property(name = "at"))]
    location: String,
}

#[cfg(test)]
mod tests {
    use insta::{assert_debug_snapshot, assert_ron_snapshot};
    use miette::Result;
    use once_cell::sync::Lazy;

    use crate::atomic_database::AtomicDatabase;

    use super::{PolymerChemistry, PolymerChemistryKdl};

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
    fn parse_muropeptide_chemistry() -> Result<()> {
        let chemistry: PolymerChemistryKdl = knuffel::parse("muropeptide_chemistry.rs", KDL)?;
        assert_debug_snapshot!(chemistry);
        Ok(())
    }

    #[test]
    fn build_muropeptide_chemistry() -> Result<()> {
        let chemistry = PolymerChemistry::from_kdl(&DB, "muropeptide_chemistry.rs", KDL)?;
        assert_ron_snapshot!(chemistry, {
            ".bonds, .modifications, .residues" => insta::sorted_redaction(),
            ".**.isotopes" => insta::sorted_redaction()
        });
        Ok(())
    }
}
