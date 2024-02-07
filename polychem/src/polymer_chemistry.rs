use std::{borrow::Cow, ops::Deref};

use knuffel::{
    ast::{Literal, TypeName},
    decode::{Context, Kind},
    errors::{DecodeError, ExpectedType},
    span::Spanned,
    traits::ErrorSpan,
    Decode, DecodeScalar,
};

use crate::{atomic_database::AtomicDatabase, ChemicalComposition};

// FIXME: Check that field names here line up with those in `lib.rs`!
#[derive(Decode, Debug)]
struct PolymerChemistryKdl {
    #[knuffel(child)]
    bonds: BondsKdl,
    #[knuffel(child)]
    modifications: ModificationsKdl,
    #[knuffel(child)]
    residues: ResiduesKdl,
}

#[derive(Decode, Debug)]
struct BondsKdl {
    #[knuffel(children)]
    bonds: Vec<BondKdl>,
}

#[derive(Decode, Debug)]
struct BondKdl {
    #[knuffel(node_name)]
    name: String,
    #[knuffel(child)]
    from: BondTargetKdl,
    #[knuffel(child)]
    to: BondTargetKdl,
    #[knuffel(child, unwrap(argument))]
    lost: ChemicalCompositionKdl,
}

#[derive(Decode, Debug)]
struct BondTargetKdl {
    #[knuffel(argument)]
    functional_group: String,
    #[knuffel(property(name = "at"))]
    location: Option<String>,
    #[knuffel(property(name = "of"))]
    residue: Option<String>,
}

#[derive(Decode, Debug)]
struct ModificationsKdl {
    #[knuffel(children(name = "targeting"))]
    targets: Vec<ModificationTargetKdl>,
}

// FIXME: Consider using flatten to abstract out common types!
#[derive(Decode, Debug)]
struct ModificationTargetKdl {
    #[knuffel(argument)]
    functional_group: String,
    #[knuffel(property(name = "at"))]
    location: Option<String>,
    #[knuffel(property(name = "of"))]
    residue: Option<String>,
    #[knuffel(children)]
    modifications: Vec<ModificationKdl>,
}

#[derive(Decode, Debug)]
struct ModificationKdl {
    #[knuffel(node_name)]
    abbr: String,
    #[knuffel(argument)]
    name: String,
    #[knuffel(child, unwrap(argument))]
    lost: Option<ChemicalCompositionKdl>,
    #[knuffel(child, unwrap(argument))]
    gained: Option<ChemicalCompositionKdl>,
}

#[derive(Decode, Debug)]
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

#[derive(Decode, Debug)]
struct ResidueKdl {
    #[knuffel(node_name)]
    residue_type: String,
    #[knuffel(argument)]
    abbr: String,
    #[knuffel(argument)]
    name: String,
    #[knuffel(child, unwrap(argument))]
    composition: ChemicalCompositionKdl,
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

#[derive(Debug, Default)]
struct ChemicalCompositionKdl(ChemicalComposition);

impl<S: ErrorSpan> DecodeScalar<S> for ChemicalCompositionKdl {
    fn type_check(type_name: &Option<Spanned<TypeName, S>>, ctx: &mut Context<S>) {
        if let Some(t) = type_name {
            ctx.emit_error(DecodeError::TypeName {
                span: t.span().clone(),
                found: Some(t.deref().clone()),
                expected: ExpectedType::no_type(),
                rust_type: "ChemicalComposition",
            });
        }
    }

    fn raw_decode(
        value: &Spanned<Literal, S>,
        ctx: &mut Context<S>,
    ) -> Result<Self, DecodeError<S>> {
        match &**value {
            Literal::String(s) => {
                let db = ctx.get::<Cow<'static, AtomicDatabase>>().unwrap();
                match ChemicalComposition::new(db, s) {
                    Ok(d) => Ok(Self(d)),
                    Err(e) => {
                        ctx.emit_error(DecodeError::conversion_diagnostic(value, Box::new(e)));
                        Ok(Self::default())
                    }
                }
            }
            Literal::Null => Ok(Self::default()),
            unsupported => {
                ctx.emit_error(DecodeError::unsupported(
                    value,
                    format!(
                        "expected a string or null, found {}",
                        Kind::from(unsupported)
                    ),
                ));
                Ok(Self::default())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use insta::assert_debug_snapshot;
    use knuffel::{decode::Context, span::Span};
    use miette::Result;
    use once_cell::sync::Lazy;

    use crate::atomic_database::AtomicDatabase;

    use super::PolymerChemistryKdl;

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
        let chemistry: PolymerChemistryKdl = knuffel::parse_with_context(
            "muropeptide_chemistry.rs",
            KDL,
            |ctx: &mut Context<Span>| ctx.set::<Cow<'static, AtomicDatabase>>(Cow::Borrowed(&DB)),
        )?;
        insta::with_settings!({filters => vec![
            (r"\d+: Isotope \{[^}]+\}", "[Isotope]"),
        ]}, {
            assert_debug_snapshot!(chemistry);
        });
        Ok(())
    }
}
