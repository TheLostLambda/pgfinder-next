use std::{iter, sync::LazyLock};

// External Crate Imports
use ahash::HashMap;
use itertools::Itertools;
use knuffel::Decode;
use miette::{Result, miette};
use regex::Regex;

use crate::{GLYCOSIDIC_BOND, Muropeptide, PEPTIDE_BOND, STEM_BOND};

// FIXME: Sort into section

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
enum Position {
    Stem(StemPosition),
    Lateral,
    Other,
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
    // FIXME: By copying all of the shared rules into every single residue, we're wasting a massive amount of memory and
    // time compiling the same `Regex`es over and over again and storing them all independently. Think about just
    // storing references to a master list of rules somewhere?
    bond_sites: BondSites,
}

// FIXME: This is another nasty short-cut... Won't be needed after we add these different isomers to the polymer
// database — this is just a temporary data structure to store residue SMILES after isomer resolution, but before
// bonding
#[derive(Debug)]
struct SmilesResidue {
    // FIXME: Should probably also keep track of the isomer name, but I'm not investing that sort of energy into this
    // hack...
    smiles: Smiles,
    bond_sites: BondSites,
}

impl SmilesResidue {
    // FIXME: `&str` should probably be something like `impl AsRef<str>`...
    fn bond(&mut self, kind: &str, acceptor: &mut Self) {
        // FIXME: This shouldn't be harccoded... I should have some dynamic way to pick an index value
        let index = match kind {
            crate::GLYCOSIDIC_BOND => 2,
            crate::PEPTIDE_BOND | crate::STEM_BOND => 3,
            crate::NTOC_BOND | crate::CTON_BOND => 4,
            crate::CROSSLINK_BOND | crate::LAT_CROSSLINK_BOND => 5,
            _ => panic!(),
        }
        .to_string();

        let replace_capture = |re: &Regex, old, new| {
            // FIXME: The first `.unwrap()` is wrong here, but again should probably be checked for during validation — in
            // the KDL, it should be ensured that the regex given by the user corresponds to exactly one match in each
            // residue it applies to!
            // SAFETY: The second `.unwrap()` is okay, since `Regex`es must have a `static_captures_len()` of `Some(2)`
            let group = re.captures(old).unwrap().get(1).unwrap();
            let prefix = &old[..group.start()];
            let suffix = &old[group.end()..];
            [prefix, new, suffix].concat()
        };

        // FIXME: That indexing can probably fail / needs upstream validation...
        let BondingRegexes {
            donor: Some(donor_re),
            ..
        } = &self.bond_sites[kind]
        else {
            // FIXME: Obviously bad
            panic!()
        };

        let BondingRegexes {
            acceptor: Some(acceptor_re),
            ..
        } = &acceptor.bond_sites[kind]
        else {
            // FIXME: Obviously bad
            panic!()
        };

        let donor_smiles = replace_capture(donor_re, &self.smiles, &index);
        let acceptor_smiles = replace_capture(acceptor_re, &acceptor.smiles, &index);

        self.smiles = donor_smiles;
        acceptor.smiles = acceptor_smiles;
    }
}

// FIXME: This feels a bit nasty... Maybe I should have a wrapper type like `Residue` in `polychem`?
impl ResidueDescription {
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

    fn smiles(&self, position: Position) -> SmilesResidue {
        // FIXME: What a wasteful copy...
        let smiles = self.isomer(position).smiles.clone();
        // FIXME: And even worse
        let bond_sites = self.bond_sites.clone();
        SmilesResidue { smiles, bond_sites }
    }
}

// Private =============================================================================================================

#[derive(Debug)]
struct ResidueTypeDescription {
    templates: Templates,
    isomer_rules: IsomerRules,
    bond_sites: BondSites,
}

// FIXME: This is the same as the Kdl schema version... Is it stupid to have this duplicated?
#[derive(Clone, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
struct Isomer {
    name: IsomerName,
    smiles: String,
}

#[derive(Clone, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
struct BondingRegexes {
    #[cfg_attr(test, serde(serialize_with = "ser_regex"))]
    donor: Option<Regex>,
    #[cfg_attr(test, serde(serialize_with = "ser_regex"))]
    acceptor: Option<Regex>,
}

// FIXME: Nasty and shouldn't live in this part of the file...
// FIXME: Would this be cleaner using the `serde_with` crate?
#[cfg(test)]
// NOTE: Serde expects this type to be `&T`, so `Option<&Regex>` isn't really an option here...
#[allow(clippy::ref_option)]
fn ser_regex<S: serde::Serializer>(value: &Option<Regex>, ser: S) -> Result<S::Ok, S::Error> {
    // Serializes `Regex` objects using their `Display` implementation
    if let Some(re) = value {
        ser.serialize_some(re.as_str())
    } else {
        ser.serialize_none()
    }
}

// Private Types =======================================================================================================

type IsomerRules = HashMap<Position, IsomerName>;
type Residues = HashMap<String, ResidueDescription>;
type Templates = Vec<(Placeholders, Isomer)>;
type ResidueTypes = HashMap<String, ResidueTypeDescription>;
type Placeholders = Vec<String>;
type BondSites = HashMap<BondKind, BondingRegexes>;
type BondKind = String;
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
    #[knuffel(children(name = "bond-sites"))]
    bond_sites: Vec<BondSitesKdl>,
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
    #[knuffel(children(name = "bond-sites"))]
    bond_sites: Vec<BondSitesKdl>,
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
struct BondSitesKdl {
    #[knuffel(argument)]
    bond: String,
    #[knuffel(property)]
    donor: Option<String>,
    #[knuffel(property)]
    acceptor: Option<String>,
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

// FIXME: Taking inspiration from how things like `nom` parsers are written, I should realy think about just turning
// all of these `impl`s into top-level functions... Maybe that will be a bit clearer?
// FIXME: This should *not* be allowed — it's just a hack for now...
#[allow(clippy::fallible_impl_from)]
impl From<ResidueTypeKdl> for ResidueTypeEntry {
    fn from(value: ResidueTypeKdl) -> Self {
        let templates = value.isomers.into_iter().map(TemplateEntry::from).collect();
        // FIXME: Stupid newtype...
        let isomer_rules = IsomerRulesKdl(value.stem, value.lateral).into();
        let bond_sites: BondSites = value
            .bond_sites
            .into_iter()
            .map(BondSitesEntry::try_from)
            // FIXME: `.unwrap()` is obviously wrong here, but I cba to make it work right now
            .collect::<Result<_, _>>()
            .unwrap();
        (value.name, ResidueTypeDescription {
            templates,
            isomer_rules,
            bond_sites,
        })
    }
}

type BondSitesEntry = (BondKind, BondingRegexes);
impl TryFrom<BondSitesKdl> for BondSitesEntry {
    // FIXME: This should probably be wrapped in a more general error type?
    type Error = regex::Error;

    fn try_from(value: BondSitesKdl) -> std::result::Result<Self, Self::Error> {
        let build_and_validate_re = |ref re_str: String| {
            let re = Regex::new(re_str);
            // NOTE: Needs exactly two "captures" — the implicit one representing the whole match, then the single capture
            // group that indicates what substring is replaced during bonding.
            // FIXME: Don't panic — return a helpful error message instead!
            assert_eq!(re.as_ref().unwrap().static_captures_len(), Some(2));

            // FIXME: Horrible and awful
            re.unwrap()
        };

        let bonding_regexes = BondingRegexes {
            donor: value.donor.map(build_and_validate_re),
            acceptor: value.acceptor.map(build_and_validate_re),
        };

        Ok((value.bond, bonding_regexes))
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
            bond_sites,
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

        // FIXME: With this duplication, maybe I should also abstract out this pattern...
        // Merge `bond_sites` from `residue_type`
        // FIXME: Yuck with the clashing names here...
        // FIXME: `.unwrap()` is obviously wrong
        let residue_bond_sites = bond_sites
            .into_iter()
            .map(|b| BondSitesEntry::try_from(b).unwrap());
        let mut bond_sites = residue_type.bond_sites.clone();
        bond_sites.extend(residue_bond_sites);

        // And remove any those don't have a corresponding isomer
        isomer_rules.retain(|_, target| isomers.iter().any(|Isomer { name, .. }| name == target));

        Ok((abbr, ResidueDescription {
            isomers,
            isomer_rules,
            bond_sites,
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
    // FIXME: Split this out into some sub-functions...
    fn muropeptide_to_smiles(&self, muropeptide: &Muropeptide) -> Option<Smiles> {
        // FIXME: This needs to be removed, but I'll need elsewhere to decide when to return `None`!
        if muropeptide.monomers.len() != 1 {
            return None;
        }

        let monomer = &muropeptide.monomers[0];

        let lookup_residue = |id| {
            let abbr = muropeptide.polymer.residue(id).unwrap().abbr();
            self.residues.get(abbr).unwrap()
        };
        // Get SMILES for glycan residues
        let mut glycan: Vec<_> = monomer
            .glycan
            .iter()
            .copied()
            .map(|r| lookup_residue(r).smiles(Position::Other))
            .collect();

        for i in 0..glycan.len().saturating_sub(1) {
            // NOTE: I love Rust, but I fucking hate it sometimes... Yet another failing of the borrow-checker means
            // that I need `split_at_mut()` just to mutate two *different* elements of a vector...
            let (head, tail) = glycan.split_at_mut(i + 1);
            let donor = &mut head[i];
            let acceptor = &mut tail[0];
            donor.bond(GLYCOSIDIC_BOND, acceptor);
        }

        // Pick isomers for all of the residues
        let mut peptide: Vec<_> = iter::zip(1.., &monomer.peptide)
            .map(|(i, aa)| lookup_residue(aa.residue).smiles(Position::Stem(i)))
            .collect();

        // Bond together stem residues
        // FIXME: Would be nice to have `windows_mut` for this sort of thing...
        // FIXME: This needs to be de-duplicated with glycan chain bonding...
        for i in 0..peptide.len().saturating_sub(1) {
            // TODO: Keep and eye on https://github.com/rust-lang/rust/pull/134633 for a better solution!
            let [donor, acceptor] = &mut peptide[i..=i + 1] else {
                unreachable!()
            };
            donor.bond(PEPTIDE_BOND, acceptor);
        }

        // Bond glycan to stem
        if let (Some(donor), Some(acceptor)) = (glycan.last_mut(), peptide.first_mut()) {
            donor.bond(STEM_BOND, acceptor);
        }

        // FIXME: Also nasty and repetitive...
        Some(
            glycan
                .into_iter()
                .map(|sr| sr.smiles)
                .chain(peptide.into_iter().map(|sr| sr.smiles))
                .join("."),
        )
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
            include_str!("../data/polymer_database.kdl"),
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
            ".**.isomer_rules, .**.bond_sites" => insta::sorted_redaction(),
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

    #[ignore]
    #[test]
    fn bond_smiles() {
        let db = SmilesDatabase::new("test_smiles_database.kdl", KDL).unwrap();
        let mut m = db.residues["m"].smiles(Position::Other);
        let mut l_ala = db.residues["A"].smiles(Position::Stem(1));
        m.bond("Stem", &mut l_ala);
        dbg!(format!("{}.{}", m.smiles, l_ala.smiles));
        panic!();
    }

    // TODO: Actually write this properly!
    // #[ignore]
    #[test]
    fn muropeptide_to_smiles() {
        let db = SmilesDatabase::new("test_smiles_database.kdl", KDL).unwrap();
        let muropeptides = ["gm", "gmgm", "AEJA", "gm-AEJA", "gmgm-A", "gm-AQKAA"];
        let smiles = muropeptides
            .map(|structure| {
                let muropeptide = Muropeptide::new(&POLYMERIZER, structure).unwrap();
                db.muropeptide_to_smiles(&muropeptide).unwrap()
            })
            .join("\n");
        assert_snapshot!(smiles);
    }
}
