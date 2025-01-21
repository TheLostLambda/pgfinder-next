use std::{iter, sync::LazyLock};

// External Crate Imports
use ahash::HashMap;
use knuffel::Decode;
use miette::{Result, miette};
use regex::Regex;

use crate::{AminoAcid, Muropeptide};

// FIXME: Sort into section

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
enum Position {
    Stem(StemPosition),
    Lateral,
}

// SMILES Database =====================================================================================================

#[derive(Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
pub struct SmilesDatabase {
    // TODO: Add `modifications`
    pub residues: Residues,
}

impl SmilesDatabase {
    pub fn new(file_name: impl AsRef<str>, kdl_text: impl AsRef<str>) -> Result<Self> {
        let parsed_db: SmilesDatabaseKdl = knuffel::parse(file_name.as_ref(), kdl_text.as_ref())?;
        parsed_db
            .try_into()
            .map_err(|()| miette!("failed to validate SMILES KDL file"))
    }
}

#[derive(Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
pub struct ResidueDescription {
    isomers: Isomers,
    isomer_rules: IsomerRules,
}

// FIXME: This feels a bit nasty... Maybe I should have a wrapper type like `Residue` in `polychem`?
impl ResidueDescription {
    /// Fetches the default isomer — used for residues without `Position` dependence
    fn smiles(&self) -> &Smiles {
        &self.isomers[0].smiles
    }

    fn isomer(&self, position: Position) -> &Isomer {
        // NOTE: Don't agree with this nursery lint making things clearer...
        #[allow(clippy::option_if_let_else)]
        if let Some(isomer_name) = self.isomer_rules.get(&position) {
            // SAFETY: The `.unwrap()` should be a-okay since we've already validated that all rules have a
            // corresponding isomer
            self.isomers
                .iter()
                .find(|Isomer { name, .. }| (name == isomer_name))
                .unwrap()
        } else {
            &self.isomers[0]
        }
    }
}

// Private =============================================================================================================

#[derive(Debug)]
struct ResidueTypeDescription {
    templates: Templates,
    isomer_rules: IsomerRules,
}

// FIXME: This is the same as the Kdl schema version... Is it stupid to have this duplicated?
#[derive(Clone, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
struct Isomer {
    name: IsomerName,
    smiles: String,
}

// Private Types =======================================================================================================

type IsomerRules = HashMap<Position, IsomerName>;
type Residues = HashMap<String, ResidueDescription>;
type Templates = Vec<(Placeholders, Isomer)>;
type ResidueTypes = HashMap<String, ResidueTypeDescription>;
type Placeholders = Vec<String>;
type Isomers = Vec<Isomer>;
type Smiles = String;
type IsomerName = Option<String>;
type StemPosition = u8;

// KDL File Schema =====================================================================================================

#[derive(Debug, Decode)]
struct SmilesDatabaseKdl {
    // TODO: Add `modifications`
    #[knuffel(child)]
    residues: ResiduesKdl,
}

// ---------------------------------------------------------------------------------------------------------------------

#[derive(Debug, Decode)]
struct ResiduesKdl {
    #[knuffel(child, unwrap(children))]
    types: Vec<ResidueTypeKdl>,
    #[knuffel(children)]
    residues: Vec<ResidueKdl>,
}

// ---------------------------------------------------------------------------------------------------------------------

#[derive(Debug, Decode)]
struct ResidueTypeKdl {
    #[knuffel(node_name)]
    name: String,
    #[knuffel(children(name = "isomer"))]
    isomers: Vec<IsomerKdl>,
    #[knuffel(children(name = "stem"))]
    stem: Vec<StemKdl>,
    #[knuffel(child)]
    lateral: Option<LateralKdl>,
    #[knuffel(children(name = "bond-site"))]
    bond_site: Vec<BondSiteKdl>,
}

#[derive(Debug, Decode)]
struct ResidueKdl {
    #[knuffel(node_name)]
    kind: String,
    #[knuffel(argument)]
    abbr: String,
    // NOTE: Duplicated from `ResidueTypeKdl` — the limitations of `knuffel::flatten` mean this can't be abstracted...
    #[knuffel(children(name = "isomer"))]
    isomers: Vec<IsomerKdl>,
    #[knuffel(children(name = "stem"))]
    stem: Vec<StemKdl>,
    #[knuffel(child)]
    lateral: Option<LateralKdl>,
    #[knuffel(children(name = "bond-site"))]
    bond_site: Vec<BondSiteKdl>,
    // END DUPLICATION
    // FIXME: What the hell should I call these?
    #[knuffel(children)]
    template_values: Vec<TemplateValueKdl>,
}

// ---------------------------------------------------------------------------------------------------------------------

#[derive(Debug, Decode)]
struct IsomerKdl {
    #[knuffel(argument)]
    name: Option<String>,
    #[knuffel(argument)]
    smiles: String,
}

#[derive(Debug, Decode)]
struct StemKdl {
    #[knuffel(arguments)]
    positions: Vec<StemPosition>,
    #[knuffel(property(name = "use"))]
    use_isomer: IsomerName,
}

#[derive(Debug, Decode)]
struct LateralKdl {
    #[knuffel(property(name = "use"))]
    use_isomer: IsomerName,
}

#[derive(Debug, Decode)]
struct BondSiteKdl {
    #[knuffel(argument)]
    regex: String,
    #[knuffel(property)]
    bond: String,
}

#[derive(Debug, Decode)]
struct TemplateValueKdl {
    #[knuffel(node_name)]
    name: String,
    #[knuffel(argument)]
    smiles: String,
}

// Contextual Validation Trait  ========================================================================================

// FIXME: Obviously this is wrong and should eventually have a real error type put here!
type SmilesResult<T> = Result<T, ()>;

trait ValidateInto<'c, T> {
    type Context: 'c;

    fn validate(self, ctx: Self::Context) -> SmilesResult<T>;
}

// Conversion From Parsed KDL to Internal Representation ===============================================================

impl TryFrom<SmilesDatabaseKdl> for SmilesDatabase {
    // FIXME: Obviously this is wrong and should eventually have a real error type put here!
    type Error = ();

    fn try_from(value: SmilesDatabaseKdl) -> Result<Self, Self::Error> {
        Ok(Self {
            residues: value.residues.try_into()?,
        })
    }
}

// TODO: No context needed currently, though in the future it might be reasonable to use the `PolymerDatabase` to
// validate this (e.g. catch cases where we've given SMILES that doesn't actually have a corresponding residue). If
// that becomes the case, then this should look more like `polychem::moieties::polymer_database::ValidateInto`.
impl TryFrom<ResiduesKdl> for Residues {
    // FIXME: Obviously this is wrong and should eventually have a real error type put here!
    type Error = ();

    fn try_from(ResiduesKdl { types, residues }: ResiduesKdl) -> Result<Self, Self::Error> {
        let types = types.into_iter().map(ResidueTypeEntry::from).collect();
        residues.into_iter().map(|r| r.validate(&types)).collect()
    }
}

type ResidueTypeEntry = (String, ResidueTypeDescription);

impl From<ResidueTypeKdl> for ResidueTypeEntry {
    fn from(value: ResidueTypeKdl) -> Self {
        let templates = value.isomers.into_iter().map(TemplateEntry::from).collect();
        // FIXME: Stupid newtype...
        let isomer_rules = IsomerRulesKdl(value.stem, value.lateral).into();
        (value.name, ResidueTypeDescription {
            templates,
            isomer_rules,
        })
    }
}

type TemplateEntry = (Placeholders, Isomer);

impl From<IsomerKdl> for TemplateEntry {
    fn from(value: IsomerKdl) -> Self {
        static PLACEHOLDER_RE: LazyLock<Regex> = LazyLock::new(|| Regex::new("<([^>]+)>").unwrap());

        let placeholders = PLACEHOLDER_RE
            .captures_iter(&value.smiles)
            .map(|c| c[1].to_owned())
            .collect();

        (placeholders, value.into())
    }
}

// FIXME: This is a `From` because that's obviously correct, but now it's inconsistent with other things that could just
// be `From` instead of `ValidateInto` — even then, this function feels silly...
impl From<IsomerKdl> for Isomer {
    fn from(IsomerKdl { name, smiles }: IsomerKdl) -> Self {
        Self { name, smiles }
    }
}

type Replacer = (Regex, Smiles);
type TemplateValues = HashMap<String, Replacer>;

// FIXME: I should probably have a trait without all of the context and error stuff... Just a `ConvertInto` or something
// NOTE: This is only `ValidateInto` because of orphan rule BS — it should just be a `From` impl...
impl ValidateInto<'_, TemplateValues> for Vec<TemplateValueKdl> {
    type Context = ();

    fn validate(self, _ctx: Self::Context) -> SmilesResult<TemplateValues> {
        Ok(self
            .into_iter()
            // FIXME: Users can inject regex here — replace `.unwrap()` with `?`
            .map(|tv| {
                let re = Regex::new(&format!("<{}>", tv.name)).unwrap();
                (tv.name, (re, tv.smiles))
            })
            .collect())
    }
}

type ResidueEntry = (String, ResidueDescription);

impl<'t> ValidateInto<'t, ResidueEntry> for ResidueKdl {
    type Context = &'t ResidueTypes;

    fn validate(self, ctx: Self::Context) -> SmilesResult<ResidueEntry> {
        let Self {
            kind,
            abbr,
            isomers,
            template_values,
            stem,
            lateral,
            ..
        } = self;

        // FIXME: Placeholder error type — replace with something informative!
        let residue_type = ctx.get(&kind).ok_or(())?;
        let template_values = template_values.validate(())?;

        // Populate `isomers` from templates + explicit isomers
        let isomers: Vec<_> = residue_type
            .templates
            .iter()
            .filter(|(placeholders, _)| {
                placeholders.iter().all(|p| template_values.contains_key(p))
            })
            .cloned()
            .map(|(_, Isomer { name, smiles })| {
                let smiles = template_values
                    .values()
                    .fold(smiles, |smiles, (re, replacement)| {
                        re.replace_all(&smiles, replacement).into_owned()
                    });
                Isomer { name, smiles }
            })
            // TODO: Really, what should be nice is having a KDL schema validator, then just using the raw AST from the
            // `kdl` reference-crate parser. That way I can maintain the order of fields and popluate the `isomers`
            // vector that way — keeping the same order in the configuration file for determining default structures.
            // I've opened an issue about this: https://github.com/kdl-org/kdl-rs/issues/116
            .chain(isomers.into_iter().map(Isomer::from))
            .collect();

        // Merge `isomer_rules` from `residue_type`
        let mut isomer_rules = residue_type.isomer_rules.clone();
        isomer_rules.extend(IsomerRules::from(IsomerRulesKdl(stem, lateral)));

        // And remove any those don't have a corresponding isomer
        isomer_rules.retain(|_, target| isomers.iter().any(|Isomer { name, .. }| name == target));

        Ok((abbr, ResidueDescription {
            isomers,
            isomer_rules,
        }))
    }
}

// FIXME: Newtype hack so I can have a `From` impl for this...
struct IsomerRulesKdl(Vec<StemKdl>, Option<LateralKdl>);
impl From<IsomerRulesKdl> for IsomerRules {
    fn from(value: IsomerRulesKdl) -> Self {
        value
            .0
            .into_iter()
            .flat_map(RuleEntries::from)
            .chain(value.1.into_iter().map(RuleEntry::from))
            .collect()
    }
}

// FIXME: Maybe could do this as an iterator without allocating a `Vec`?
type RuleEntry = (Position, IsomerName);
type RuleEntries = Vec<RuleEntry>;
impl From<StemKdl> for RuleEntries {
    fn from(value: StemKdl) -> Self {
        value
            .positions
            .into_iter()
            .map(|p| (Position::Stem(p), value.use_isomer.clone()))
            .collect()
    }
}

impl From<LateralKdl> for RuleEntry {
    fn from(value: LateralKdl) -> Self {
        (Position::Lateral, value.use_isomer)
    }
}

// SOME REALLY HACKY BULLSHIT ==========================================================================================

impl SmilesDatabase {
    fn muropeptide_to_smiles(&self, muropeptide: &Muropeptide) -> Smiles {
        assert_eq!(muropeptide.monomers.len(), 1);
        let monomer = &muropeptide.monomers[0];
        assert!(monomer.glycan.is_empty());
        let peptide = &monomer.peptide;

        // FIXME: Yikes
        let mut residues = Vec::new();
        // FIXME: Refactor to iterator chain?
        for (
            i,
            AminoAcid {
                residue,
                lateral_chain,
            },
        ) in iter::zip(1.., peptide)
        {
            assert!(lateral_chain.is_none());
            // FIXME: Donno if that `.unwrap()` is safe, lol
            let residue = muropeptide.polymer.residue(*residue).unwrap();
            let smiles = self
                .residues
                .get(residue.abbr())
                .unwrap()
                .isomer(Position::Stem(i))
                .smiles
                .clone();
            residues.push(smiles);
        }
        // FIXME: Yuck
        for i in 0..residues.len() - 1 {
            residues[i] = residues[i][..residues[i].len() - 1].to_string();
        }

        residues.concat()
    }
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use std::iter;

    use insta::{assert_debug_snapshot, assert_ron_snapshot, assert_snapshot};
    use itertools::Itertools;
    use once_cell::sync::Lazy;
    use polychem::{AtomicDatabase, PolymerDatabase, Polymerizer};

    static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);
    static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "polymer_database.kdl",
            include_str!("../tests/data/polymer_database.kdl"),
        )
        .unwrap()
    });

    static POLYMERIZER: Lazy<Polymerizer> = Lazy::new(|| Polymerizer::new(&ATOMIC_DB, &POLYMER_DB));

    use super::*;

    const KDL: &str = include_str!("../data/smiles_database.kdl");

    #[test]
    fn parse_smiles_database() {
        let db: SmilesDatabaseKdl = knuffel::parse("test_smiles_database.kdl", KDL).unwrap();
        assert_debug_snapshot!(db);
    }

    #[test]
    fn build_smiles_database() {
        let db = SmilesDatabase::new("test_smiles_database.kdl", KDL).unwrap();
        assert_ron_snapshot!(db, {
            ".residues" => insta::sorted_redaction(),
            ".**.isomer_rules" => insta::sorted_redaction(),
        });
    }

    #[test]
    fn dump_isomer_smiles() {
        let db = SmilesDatabase::new("test_smiles_database.kdl", KDL).unwrap();
        let dump = db
            .residues
            .iter()
            .sorted_unstable_by(|(a, _), (b, _)| a.cmp(b))
            .flat_map(|(residue, ResidueDescription { isomers, .. })| {
                isomers.iter().map(move |Isomer { name, smiles }| {
                    // NOTE: Not convinced that lint makes things more readable, and `nursery` means I can ignore it
                    #[allow(clippy::option_if_let_else)]
                    if let Some(isomer) = name {
                        format!("{residue}({isomer}): {smiles}")
                    } else {
                        format!("{residue}: {smiles}")
                    }
                })
            })
            .join("\n");
        assert_snapshot!(dump);
    }

    #[test]
    fn dump_isomer_positions() {
        let db = SmilesDatabase::new("test_smiles_database.kdl", KDL).unwrap();
        let dump = db
            .residues
            .iter()
            .sorted_unstable_by(|(a, _), (b, _)| a.cmp(b))
            .flat_map(|(residue, description)| {
                let positions = (1..=5).map(Position::Stem).chain([Position::Lateral]);
                let isomers = positions.map(move |p| {
                    let Isomer { name, smiles } = description.isomer(p);
                    // FIXME: Change `IsomerName` to be just a `String` — in the KDL file, replace `null` with `L/D` for
                    // glycine. It's not worth the complication of having an `Option` everywhere!
                    let isomer = name.as_deref().unwrap_or("L/D");
                    format!("    {p:?} -> {isomer:5} ({smiles})")
                });

                iter::once(format!("{residue}:")).chain(isomers)
            })
            .join("\n");
        assert_snapshot!(dump);
    }

    // TODO: Actually write this properly!
    #[ignore]
    #[test]
    fn muropeptide_to_smiles() {
        let db = SmilesDatabase::new("test_smiles_database.kdl", KDL).unwrap();
        let muropeptide = Muropeptide::new(&POLYMERIZER, "AEJA").unwrap();
        dbg!(db.muropeptide_to_smiles(&muropeptide));
        panic!();
    }
}
