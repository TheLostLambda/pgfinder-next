use std::{collections::HashMap, ops::Deref};

use knuffel::{
    ast::{Literal, TypeName},
    decode::{Context, Kind},
    errors::{DecodeError, ExpectedType},
    span::{Span, Spanned},
    traits::ErrorSpan,
    Decode, DecodeScalar,
};
use miette::Result;

use crate::{
    atomic_database::AtomicDatabase, ChemicalComposition, FunctionalGroup, GroupState, Location,
    Residue,
};

// FIXME: Might want to change how this is structured down the line...
#[derive(Clone, PartialEq, Eq, Debug)]
struct PolymerChemistry {
    pub(super) bonds: HashMap<String, BondDescription>,
    pub(super) modifications: HashMap<String, ModificationDescription>,
    pub(super) residues: HashMap<String, Residue>,
}

// impl PolymerChemistry {
//     // FIXME: Consider Cow<'static, AtomicDatabase>
//     pub fn from_kdl(
//         db: AtomicDatabase,
//         file_name: impl AsRef<str>,
//         text: impl AsRef<str>,
//     ) -> Result<Self> {
//         let parsed_chemistry: PolymerChemistryKdl = knuffel::parse_with_context(
//             file_name.as_ref(),
//             text.as_ref(),
//             |ctx: &mut Context<Span>| ctx.set(db),
//         )?;
//         Ok(parsed_chemistry.into())
//     }
// }

// impl From<PolymerChemistryKdl> for PolymerChemistry {
//     fn from(value: PolymerChemistryKdl) -> Self {
//         Self {
//             bonds: value.bonds.into(),
//             modifications: value.modifications.into(),
//             residues: value.residues.into(),
//         }
//     }
// }

// impl From<ResiduesKdl> for HashMap<String, Residue> {
//     fn from(value: ResiduesKdl) -> Self {
//         let types: HashMap<_, _> = value
//             .types
//             .into_iter()
//             .map(ResidueTypeEntry::from)
//             .collect();
//         let to_residue = |r: ResidueKdl| {
//             let mut functional_groups: HashMap<_, _> = r
//                 .functional_groups
//                 .into_iter()
//                 .map(FunctionalGroupEntry::from)
//                 .collect();
//             // FIXME: Convert to TryFrom and use .get() here! The error of an invalid type should be reported!
//             // FIXME: Actually, make the Strings that are part of FunctionalGroupKdl and TargetKdl newtypes — I can
//             // then parse the `residues` first, use ctx.set to add all of the functional groups to a HashMap, then when
//             // I'm parsing the `modifications` and `bonds` sections, I can check that they all refer to real groups
//             // whilst their spans are available.
//             // FIXME: Even newer plan: actually just use knuffel's support for capturing a node's Span, then I can
//             // handle all of the validation in my own TryFrom blocks. Excitingly, if I take that error reporting into
//             // my own hands this way, I might actually be able to move back to mainline knuffel, since all of my
//             // changes have been to do with fixing some buggy error propagation and limiations of trying to do
//             // validation in the context of knuffel's traits
//             functional_groups.extend(types[&r.residue_type].clone());
//             Residue {
//                 id: 0,
//                 abbr: r.abbr,
//                 name: r.name,
//                 composition: r.composition.0,
//                 functional_groups,
//                 offset_modifications: Vec::new(),
//             }
//         };
//         value
//             .residues
//             .into_iter()
//             .map(|r| (r.abbr.clone(), to_residue(r)))
//             .collect()
//     }
// }

// type ResidueTypeEntry = (String, HashMap<Location, FunctionalGroup>);

// impl From<ResidueTypeKdl> for ResidueTypeEntry {
//     fn from(value: ResidueTypeKdl) -> Self {
//         let functional_groups = value
//             .functional_groups
//             .into_iter()
//             .map(FunctionalGroupEntry::from)
//             .collect();
//         (value.name, functional_groups)
//     }
// }

// type FunctionalGroupEntry = (Location, FunctionalGroup);

// impl From<FunctionalGroupKdl> for FunctionalGroupEntry {
//     fn from(value: FunctionalGroupKdl) -> Self {
//         let functional_group = FunctionalGroup {
//             name: value.name,
//             state: GroupState::default(),
//         };
//         (value.location, functional_group)
//     }
// }

// impl From<ModificationsKdl> for HashMap<String, ModificationDescription> {
//     fn from(value: ModificationsKdl) -> Self {
//         // FIXME: This is a bit repetitive with the above...
//         value
//             .modifications
//             .into_iter()
//             .map(|m| (m.abbr.clone(), m.into()))
//             .collect()
//     }
// }

// impl From<ModificationKdl> for ModificationDescription {
//     fn from(value: ModificationKdl) -> Self {
//         Self {
//             abbr: value.abbr,
//             name: value.name,
//             lost: value.lost.map(|l| l.0),
//             gained: value.gained.map(|g| g.0),
//             targets: value.targets.into_iter().map(Target::from).collect(),
//         }
//     }
// }

// impl From<BondsKdl> for HashMap<String, BondDescription> {
//     fn from(value: BondsKdl) -> Self {
//         // FIXME: This is a bit repetitive with the above...
//         value
//             .bonds
//             .into_iter()
//             .map(|b| (b.name.clone(), b.into()))
//             .collect()
//     }
// }

// impl From<BondKdl> for BondDescription {
//     fn from(value: BondKdl) -> Self {
//         Self {
//             name: value.name,
//             from: value.from.into(),
//             to: value.to.into(),
//             lost: value.lost.0,
//         }
//     }
// }

// impl From<TargetKdl> for Target {
//     fn from(value: TargetKdl) -> Self {
//         Self {
//             functional_group: value.functional_group,
//             location: value.location,
//             residue: value.residue,
//         }
//     }
// }

// FIXME: Not 110% sold on using the `...Description` structs and just `Residue`
// FIXME: Is it also silly to duplicate so much of the Kdl structs? Or is this good decoupling of the knuffel forms and
// the internal representation? With a TryFrom to shimmy between?
#[derive(Clone, PartialEq, Eq, Debug)]
struct BondDescription {
    name: String,
    from: Target,
    to: Target,
    lost: ChemicalComposition,
}

#[derive(Clone, PartialEq, Eq, Debug)]
struct ModificationDescription {
    abbr: String,
    name: String,
    lost: Option<ChemicalComposition>,
    gained: Option<ChemicalComposition>,
    targets: Vec<Target>,
}

#[derive(Clone, PartialEq, Eq, Debug)]
struct Target {
    functional_group: String,
    location: Option<String>,
    residue: Option<String>,
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

// FIXME: Naming!
// FIXME: Be sure to write tests to check this errors for missing composition keys but works for nulled ones — also
// check that the `lost` and `gained` keys are the opposite — null isn't allowed, but they can be left out entirely
// NOTE: This forces composition to have a `null` value instead of being a totally optional key
type OrNull<T> = Option<T>;

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
    composition: OrNull<ChemicalCompositionKdl>,
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

// #[derive(Debug, Default)]
// struct ChemicalCompositionKdl(ChemicalComposition);
type ChemicalCompositionKdl = Spanned<String, Span>;

// impl<S: ErrorSpan> DecodeScalar<S> for ChemicalCompositionKdl {
//     fn type_check(type_name: &Option<Spanned<TypeName, S>>, ctx: &mut Context<S>) {
//         if let Some(t) = type_name {
//             ctx.emit_error(DecodeError::TypeName {
//                 span: t.span().clone(),
//                 found: Some(t.deref().clone()),
//                 expected: ExpectedType::no_type(),
//                 rust_type: "ChemicalComposition",
//             });
//         }
//     }

//     fn raw_decode(
//         value: &Spanned<Literal, S>,
//         ctx: &mut Context<S>,
//     ) -> Result<Self, DecodeError<S>> {
//         match &**value {
//             Literal::String(s) => {
//                 let db = ctx.get().unwrap();
//                 match ChemicalComposition::new(db, s) {
//                     Ok(d) => Ok(Self(d)),
//                     Err(e) => {
//                         ctx.emit_error(DecodeError::conversion_diagnostic(value, Box::new(e)));
//                         Ok(Self::default())
//                     }
//                 }
//             }
//             Literal::Null => Ok(Self::default()),
//             unsupported => {
//                 ctx.emit_error(DecodeError::unsupported(
//                     value,
//                     format!(
//                         "expected a string or null, found {}",
//                         Kind::from(unsupported)
//                     ),
//                 ));
//                 Ok(Self::default())
//             }
//         }
//     }
// }

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use knuffel::{decode::Context, span::Span};
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
        let chemistry: PolymerChemistryKdl = knuffel::parse_with_context(
            "muropeptide_chemistry.rs",
            KDL,
            |ctx: &mut Context<Span>| ctx.set(DB.clone()),
        )?;
        insta::with_settings!({filters => vec![
            (r"\d+: Isotope \{[^}]+\}", "[Isotope]"),
        ]}, {
            assert_debug_snapshot!(chemistry);
        });
        Ok(())
    }

    // #[test]
    // fn build_muropeptide_chemistry() -> Result<()> {
    //     let chemistry = PolymerChemistry::from_kdl(DB.clone(), "muropeptide_chemistry.rs", KDL)?;
    //     insta::with_settings!({filters => vec![
    //         (r"\d+: Isotope \{[^}]+\}", "[Isotope]"),
    //     ]}, {
    //         assert_debug_snapshot!(chemistry);
    //     });
    //     Ok(())
    // }
}
