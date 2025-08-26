use std::{iter, sync::LazyLock};

// External Crate Imports
use ahash::HashMap;
use itertools::Itertools;
use knuffel::Decode;
use miette::{Result, miette};
use polychem::{ModificationInfo, Polymer};
use regex::Regex;

use crate::{
    CROSSLINK_BOND, CTON_BOND, Connection, CrosslinkDescriptor, GLYCOSIDIC_BOND, Monomer,
    Muropeptide, NTOC_BOND, PEPTIDE_BOND, PeptideDirection, STEM_BOND,
};

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

#[derive(Clone, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
pub struct ResidueDescription {
    isomers: Isomers,
    isomer_rules: IsomerRules,
    // FIXME: By copying all of the shared rules into every single residue, we're wasting a massive amount of memory and
    // time compiling the same `Regex`es over and over again and storing them all independently. Think about just
    // storing references to a master list of rules somewhere?
    #[cfg_attr(test, serde(serialize_with = "ser_bond_sites"))]
    bond_sites: BondSites,
    modifications: Modifications,
}

// FIXME: Nasty and shouldn't live in this part of the file...
// FIXME: Would this be cleaner using the `serde_with` crate?
#[cfg(test)]
fn ser_bond_sites<S: serde::Serializer>(value: &BondSites, ser: S) -> Result<S::Ok, S::Error> {
    // Serializes `Regex` objects using their `Display` implementation
    ser.collect_map(value.iter().map(|(k, v)| (k, v.to_string())))
}

// FIXME: This is another nasty short-cut... Won't be needed after we add these different isomers to the polymer
// database — this is just a temporary data structure to store residue SMILES after isomer resolution, but before
// bonding
#[derive(Clone, Debug)]
struct SmilesResidue {
    // FIXME: Should probably also keep track of the isomer name, but I'm not investing that sort of energy into this
    // hack...
    smiles: Smiles,
    bond_sites: BondSites,
}

#[derive(Clone, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
struct Modification {
    #[cfg_attr(test, serde(serialize_with = "ser_regex"))]
    replace: Regex,
    with: String,
}

// FIXME: Nasty and shouldn't live in this part of the file...
// FIXME: Would this be cleaner using the `serde_with` crate?
#[cfg(test)]
fn ser_regex<S: serde::Serializer>(value: &Regex, ser: S) -> Result<S::Ok, S::Error> {
    // Serializes `Regex` objects using their `Display` implementation
    ser.serialize_str(value.as_str())
}

// FIXME: Vile copy-paste just to switch out the residue numbers (`UnbranchedAminoAcid`) for `SmilesResidue`...
#[derive(Debug)]
struct AminoAcid {
    residue: SmilesResidue,
    lateral_chain: Option<LateralChain>,
}

#[derive(Debug)]
struct LateralChain {
    direction: PeptideDirection,
    peptide: Vec<SmilesResidue>,
}

impl SmilesResidue {
    // FIXME: `&str` should probably be something like `impl AsRef<str>`...
    // FIXME: Really really should not have that `index` argument... At least certainly not in the way it's currently
    // implemented...
    // FIXME: Going forward, I'll need delay bonding regex searches until all of the bonds to a residue are settled /
    // determined (which is already the case if I were looking at the functional groups, by the way...), then I can
    // loop through the residues in the order they will be printed to the SMILES string, and, importantly, then also get
    // the order of each functional group in that residue, then I can keep a running tally as I move through each place
    // in the whole string where an index needs to be inserted, so I always know (in string order) what the lowest free
    // index is... Or another way to do this, don't insert these bonds in pairs at all — just go through one residue at
    // a time, check it's functional groups for any bonds, then get the ID of the residue it's bonding to, and assign it
    // the lowest free index value. Later on, when you come across the residue that has that ID, look up the index you
    // assigned it earlier, then free that index back up for use! In even more detail, when you come across something
    // like `Donor(id)`, then immediately check for `Acceptor(id)` in the `HashMap` (and visa versa — looking for
    // `Donor(id)` when you come across `Acceptor(id)`; obviously those ids have to match — just index the map). If that
    // key is present in the `HashMap`, then remove it, insert the stored SMILES ID into the string. If that key isn't
    // already present, then scan the `HashMap::values()` for the lowest free index, and insert that key with the new
    // lowest found index. This works because if we find something like `Donor(id)` on a group, then we *know* that the
    // matching functional group will be `Acceptor(id)`, so we know exactly where to put away that SMILES index for
    // later
    fn bond(&mut self, kind: &str, acceptor: &mut Self, mut index: u8) {
        // FIXME: This shouldn't be hardccoded... I should have some dynamic way to pick an index value
        // FIXME: The number skipping is a hack because of the other code elsewhere that can modify the `index` and
        // we need to avoid index collisions here...
        index += match kind {
            crate::GLYCOSIDIC_BOND => 10,
            crate::STEM_BOND | crate::PEPTIDE_BOND => 20,
            crate::NTOC_BOND | crate::CTON_BOND => 30,
            crate::CROSSLINK_BOND => 40,
            _ => panic!(),
        };

        let index = if index < 10 {
            index.to_string()
        } else if index < 100 {
            format!("%{index}")
        } else {
            // FIXME: Don't panic...
            panic!("index too large");
        };

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

        // FIXME: Really dumb and awful clone (in the form of `to_string()`)...
        let donor_re = &self.bond_sites[&BondSite::Donor(kind.to_string())];
        let acceptor_re = &acceptor.bond_sites[&BondSite::Acceptor(kind.to_string())];

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
        #[expect(clippy::option_if_let_else)]
        if let Some(isomer_name) = self.isomer_rules.get(&position) {
            // SAFETY: The `.unwrap()` should be a-okay since we've already validated that all rules have a
            // corresponding isomer
            self.isomers
                .iter()
                .find(|Isomer { name, .. }| name == isomer_name)
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
    modifications: Modifications,
}

// FIXME: This is the same as the Kdl schema version... Is it stupid to have this duplicated?
#[derive(Clone, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
struct Isomer {
    name: IsomerName,
    smiles: String,
}

#[derive(Clone, Eq, PartialEq, Hash, Debug)]
#[cfg_attr(test, derive(serde::Serialize))]
enum BondSite {
    Acceptor(BondKind),
    Donor(BondKind),
}

// Private Types =======================================================================================================

type IsomerRules = HashMap<Position, IsomerName>;
type Residues = HashMap<String, ResidueDescription>;
type Modifications = HashMap<String, Vec<Modification>>;
type Templates = Vec<(Placeholders, Isomer)>;
type ResidueTypes = HashMap<String, ResidueTypeDescription>;
type Placeholders = Vec<String>;
type BondSites = HashMap<BondSite, Regex>;
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
    #[knuffel(children(name = "modification"))]
    modifications: Vec<ModificationKdl>,
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
    #[knuffel(children(name = "modification"))]
    modifications: Vec<ModificationKdl>,
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
    // FIXME: This is hiding a god-awful hack — on x86_64-pc-windows-msvc, `encoder_unicode` is added as a dependency
    // by `console`, added by `insta`. When that crate is included, since it provides:
    // `impl FromIterator<encode_unicode::utf8_char::Utf8Char> for Vec<u8>`, Rust gets confused and the `Decode` derive
    // macro of `knuffel` is broken. Rust can't choose between that impl and `impl<T> FromIterator<T> for Vec<T>`,
    // which it should be using. That's why this is a u32 instead of a u8 — it dodges this particular instance of the
    // issue without requiring me to make some likely very complex fixes in `knuffel`... Note you can't use u16 either,
    // since Windows has a 16-bit unicode version too...
    #[knuffel(arguments)]
    positions: Vec<u32>,
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
struct ModificationKdl {
    #[knuffel(argument)]
    abbr: String,
    #[knuffel(property)]
    replace: Option<String>,
    #[knuffel(property)]
    with: Option<String>,
    #[knuffel(children(name = "replace"))]
    replacements: Vec<ReplaceKdl>,
}

#[derive(Debug, Decode)]
struct TemplateValueKdl {
    #[knuffel(node_name)]
    name: String,
    #[knuffel(argument)]
    smiles: String,
}

#[derive(Debug, Decode)]
struct ReplaceKdl {
    // FIXME: This is sorta duplicated from `ModificationKdl`?
    #[knuffel(argument)]
    replace: String,
    #[knuffel(property)]
    with: String,
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
#[expect(clippy::fallible_impl_from)]
impl From<ResidueTypeKdl> for ResidueTypeEntry {
    fn from(value: ResidueTypeKdl) -> Self {
        let templates = value.isomers.into_iter().map(TemplateEntry::from).collect();
        // FIXME: Stupid newtype...
        let isomer_rules = IsomerRulesKdl(value.stem, value.lateral).into();
        let bond_sites: BondSites = value
            .bond_sites
            .into_iter()
            .map(BondSitesEntries::try_from)
            // FIXME: `.unwrap()` is obviously wrong here, but I cba to make it work right now
            .collect::<Result<Vec<_>, _>>()
            .unwrap()
            .into_iter()
            .flatten()
            .collect();
        let modifications = value
            .modifications
            .into_iter()
            .map(|m| ModificationEntry::try_from(m).unwrap())
            .collect();
        (
            value.name,
            ResidueTypeDescription {
                templates,
                isomer_rules,
                bond_sites,
                modifications,
            },
        )
    }
}

type BondSitesEntries = Vec<(BondSite, Regex)>;
impl TryFrom<BondSitesKdl> for BondSitesEntries {
    // FIXME: This should probably be wrapped in a more general error type?
    type Error = regex::Error;

    fn try_from(value: BondSitesKdl) -> Result<Self, Self::Error> {
        let build_and_validate_re = |ref re_str: String| {
            let re = Regex::new(re_str);
            // NOTE: Needs exactly two "captures" — the implicit one representing the whole match, then the single capture
            // group that indicates what substring is replaced during bonding.
            // FIXME: Don't panic — return a helpful error message instead!
            assert_eq!(re.as_ref().unwrap().static_captures_len(), Some(2));

            // FIXME: Horrible and awful
            re.unwrap()
        };

        // FIXME: Replace with iterators
        let mut result = Self::new();
        if let Some(re_str) = value.donor {
            result.push((
                BondSite::Donor(value.bond.clone()),
                build_and_validate_re(re_str),
            ));
        }
        if let Some(re_str) = value.acceptor {
            result.push((
                BondSite::Acceptor(value.bond),
                build_and_validate_re(re_str),
            ));
        }
        Ok(result)
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
            modifications,
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
            .flat_map(|b| BondSitesEntries::try_from(b).unwrap());
        let mut bond_sites = residue_type.bond_sites.clone();
        bond_sites.extend(residue_bond_sites);

        // FIXME: With this duplication, maybe I should also abstract out this pattern...
        // Merge `modifications` from `residue_type`
        // FIXME: `.unwrap()` is wrong...
        let residue_modifications = modifications
            .into_iter()
            .map(|m| ModificationEntry::try_from(m).unwrap());
        let mut modifications = residue_type.modifications.clone();
        modifications.extend(residue_modifications);

        // And remove any those don't have a corresponding isomer
        isomer_rules.retain(|_, target| isomers.iter().any(|Isomer { name, .. }| name == target));

        Ok((
            abbr,
            ResidueDescription {
                isomers,
                isomer_rules,
                bond_sites,
                modifications,
            },
        ))
    }
}

type ModificationEntry = (String, Vec<Modification>);
impl TryFrom<ModificationKdl> for ModificationEntry {
    type Error = regex::Error;

    fn try_from(value: ModificationKdl) -> Result<Self, Self::Error> {
        let modifications = if let ModificationKdl {
            replace: Some(replace),
            with: Some(with),
            ..
        } = value
        {
            let replace = Regex::new(&replace)?;
            vec![Modification { replace, with }]
        } else {
            value
                .replacements
                .into_iter()
                .map(|r| {
                    let replace = Regex::new(&r.replace)?;
                    let with = r.with;
                    // FIXME: Probably a better way, but maybe I'll just need to wait for try blocks...
                    Ok(Modification { replace, with })
                })
                .collect::<Result<_, _>>()?
        };
        Ok((value.abbr, modifications))
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
            // FIXME: See the horrific FIXME in `StemKdl`
            .map(|p| {
                (
                    Position::Stem(p.try_into().unwrap()),
                    value.use_isomer.clone(),
                )
            })
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
    pub fn muropeptide_to_smiles(&self, muropeptide: &Muropeptide) -> Option<Smiles> {
        let mut monomers: Vec<_> = muropeptide
            .monomers
            .iter()
            .map(|m| self.build_monomer(&muropeptide.polymer, m))
            .collect::<Option<_>>()?;

        assert_eq!(monomers.len() - 1, muropeptide.connections.len());

        for i in 0..monomers.len() - 1 {
            let connection = &muropeptide.connections[i];
            let [left, right] = &mut monomers[i..=i + 1] else {
                unreachable!();
            };
            match connection {
                Connection::GlycosidicBond => {
                    // FIXME: This is dreadful — that needs to be a struct with named fields...
                    let (left_glycan, _) = left;
                    let (right_glycan, _) = right;

                    // SAFETY: If we made it this far, then the parser will have ensured that glycosidic bonds are only
                    // between two non-empty glycan chains
                    let donor = left_glycan.last_mut().unwrap();
                    let acceptor = right_glycan.first_mut().unwrap();
                    // FIXME: Yikes, bro, that index...
                    donor.bond(GLYCOSIDIC_BOND, acceptor, 0);
                }
                Connection::Crosslink(crosslinks) => {
                    // TODO: The parser doesn't allows multiply-crosslinked monomers yet
                    assert_eq!(crosslinks.len(), 1);
                    if let CrosslinkDescriptor::DonorAcceptor(l, r) = crosslinks[0] {
                        // FIXME: This is dreadful — that needs to be a struct with named fields...
                        let (_, left_stem) = left;
                        let (_, right_stem) = right;

                        // TODO: Currently I don't *think* my parser has any way for the end of a lateral chain to
                        // donate a crosslink bond, so I'll only implement the lateral-chain logic for the acceptor for
                        // now...
                        let donor = &mut left_stem[l as usize - 1].residue;
                        let acceptor = &mut right_stem[r as usize - 1];
                        if let Some(LateralChain { peptide, .. }) = &mut acceptor.lateral_chain {
                            donor.bond(CROSSLINK_BOND, peptide.last_mut().unwrap(), 0);
                        } else {
                            donor.bond(CROSSLINK_BOND, &mut acceptor.residue, 0);
                        }
                    } else {
                        // TODO: The parser doesn't support different bonding directions yet
                        todo!()
                    }
                }
                // TODO: The parser doesn't yet support this, so don't bother trying to serialize it!
                Connection::Both(_) => todo!(),
            }
        }

        // FIXME: Also nasty and repetitive...
        Some(
            monomers
                .into_iter()
                .flat_map(|(glycan, peptide)| {
                    glycan
                        .into_iter()
                        .map(|sr| sr.smiles)
                        .chain(peptide.into_iter().flat_map(|aa| {
                            iter::once(aa.residue.smiles).chain(
                                aa.lateral_chain
                                    .into_iter()
                                    .flat_map(|lc| lc.peptide.into_iter().map(|sr| sr.smiles)),
                            )
                        }))
                })
                .join("."),
        )
    }

    fn build_monomer(
        &self,
        polymer: &Polymer,
        monomer: &Monomer,
    ) -> Option<(Vec<SmilesResidue>, Vec<AminoAcid>)> {
        let lookup_residue = |id| {
            let residue = polymer.residue(id).unwrap();
            if residue.offset_modifications().count() > 0 {
                return None;
            }

            // FIXME: The fact this needs to be mutable is very upsetting
            let mut residue_description = self.residues.get(residue.abbr())?.clone();

            // FIXME: This is awful... I should really really have the isomer picked before applying the modifications,
            // but this will have to do for now (even though it's a bunch of wasteful computation)... And after these
            // modifications have been applied, there is no need to keep that field around in `ResidueDescription`!
            for id in residue.named_modifications() {
                let ModificationInfo::Named(modification, _) = polymer.modification(id).unwrap()
                else {
                    unreachable!()
                };
                let modifications = residue_description
                    .modifications
                    .get_mut(modification.abbr())?;
                for Modification { replace, with } in modifications {
                    // FIXME: This is the worst part of all of this...
                    for Isomer { smiles, .. } in &mut residue_description.isomers {
                        *smiles = replace.replace_all(smiles, with.as_str()).into_owned();
                    }
                }
            }

            Some(residue_description)
        };
        // Get SMILES for glycan residues
        let mut glycan: Vec<_> = monomer
            .glycan
            .iter()
            .copied()
            .map(|r| lookup_residue(r).map(|sr| sr.smiles(Position::Other)))
            .collect::<Option<_>>()?;

        bond_chain(GLYCOSIDIC_BOND, &mut glycan, 0);

        // Pick isomers for all of the residues
        let mut peptide: Vec<_> = iter::zip(1.., &monomer.peptide)
            .map(|(i, aa)| {
                let residue = lookup_residue(aa.residue)?.smiles(Position::Stem(i));
                let lateral_chain = aa.lateral_chain.as_ref().and_then(
                    |&crate::LateralChain {
                         direction,
                         ref peptide,
                     }| {
                        let peptide = peptide
                            .iter()
                            .copied()
                            .map(|r| lookup_residue(r).map(|sr| sr.smiles(Position::Lateral)))
                            .collect::<Option<_>>()?;
                        Some(LateralChain { direction, peptide })
                    },
                );
                Some(AminoAcid {
                    residue,
                    lateral_chain,
                })
            })
            .collect::<Option<_>>()?;

        // Bond together stem residues
        for i in 0..peptide.len().saturating_sub(1) {
            // FIXME: Duplicated again with `bond_chain()`!
            // FIXME: Would be nice to have `windows_mut` for this sort of thing...
            // TODO: Keep and eye on https://github.com/rust-lang/rust/pull/134633 for a better solution!
            let [
                AminoAcid {
                    residue: donor,
                    lateral_chain,
                },
                AminoAcid {
                    residue: acceptor, ..
                },
            ] = &mut peptide[i..=i + 1]
            else {
                unreachable!()
            };
            donor.bond(PEPTIDE_BOND, acceptor, u8::try_from(i % 2).unwrap());
            if let Some(lateral_chain) = lateral_chain {
                match lateral_chain.direction {
                    // FIXME: Really I should have a different enum for post-parsing that doesn't even have this
                    // `Unspecified` value... Obviously I need to parse it, but after that it shouldn't exist!
                    // FIXME: Hey, the smart thing to do there is to make it an `Option<PeptideDirection>` when parsing,
                    // then I can unwrap that after it's been determined?
                    PeptideDirection::Unspecified => panic!("This should be set during parsing!"),
                    PeptideDirection::CToN => {
                        // FIXME: Suspicious .unwrap()
                        lateral_chain
                            .peptide
                            .first_mut()
                            .unwrap()
                            .bond(CTON_BOND, donor, 2);
                        // FIXME: Don't like the mutation here
                        let lateral_chain2 = lateral_chain.peptide.iter_mut().rev();
                        bond_chain(PEPTIDE_BOND, lateral_chain2, 2);
                    }
                    PeptideDirection::NToC => {
                        // FIXME: Awful copy-paste from above!
                        // FIXME: Suspicious .unwrap()
                        donor.bond(NTOC_BOND, lateral_chain.peptide.first_mut().unwrap(), 2);
                        // FIXME: Don't like the mutation here
                        let lateral_chain = lateral_chain.peptide.iter_mut();
                        bond_chain(PEPTIDE_BOND, lateral_chain, 2);
                    }
                }
            }
        }

        // Bond glycan to stem
        if let (
            Some(donor),
            Some(AminoAcid {
                residue: acceptor, ..
            }),
        ) = (glycan.last_mut(), peptide.first_mut())
        {
            donor.bond(STEM_BOND, acceptor, 0);
        }

        // FIXME: Maybe that needs its own type...
        Some((glycan, peptide))
    }
}

// FIXME: Really shouldn't have that `lateral` argument!
fn bond_chain<'a>(
    kind: &str,
    residues: impl IntoIterator<Item = &'a mut SmilesResidue>,
    index: u8,
) {
    let mut residues: Vec<_> = residues.into_iter().collect();
    for i in 0..residues.len().saturating_sub(1) {
        // FIXME: Would be nice to have `windows_mut` for this sort of thing...
        // TODO: Keep and eye on https://github.com/rust-lang/rust/pull/134633 for a better solution!
        let [donor, acceptor] = &mut residues[i..=i + 1] else {
            unreachable!()
        };
        donor.bond(kind, acceptor, index + u8::try_from(i % 2).unwrap());
    }
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use std::iter;

    use insta::{assert_debug_snapshot, assert_ron_snapshot, assert_snapshot};
    use itertools::Itertools;
    use polychem::{AtomicDatabase, PolymerDatabase, Polymerizer};
    use std::sync::LazyLock;

    static ATOMIC_DB: LazyLock<AtomicDatabase> = LazyLock::new(AtomicDatabase::default);
    static POLYMER_DB: LazyLock<PolymerDatabase> = LazyLock::new(|| {
        PolymerDatabase::new(
            &ATOMIC_DB,
            "polymer_database.kdl",
            include_str!("../data/polymer_database.kdl"),
        )
        .unwrap()
    });

    static POLYMERIZER: LazyLock<Polymerizer> =
        LazyLock::new(|| Polymerizer::new(&ATOMIC_DB, &POLYMER_DB));

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
            ".**.isomer_rules, .**.bond_sites, .**.modifications" => insta::sorted_redaction(),
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
                    #[expect(clippy::option_if_let_else)]
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
    #[test]
    fn muropeptide_to_smiles() {
        let db = SmilesDatabase::new("test_smiles_database.kdl", KDL).unwrap();
        let muropeptides = [
            "gm",
            "gmgm",
            "AEJA",
            "gm-AEJA",
            "gmgm-A",
            "gm-AQKAA",
            "gm-AQK[AA]AA",
            "gm-AQ[GG]K[AA]AA",
            "gm-AE[A]J[GG]A",
            "gm-AEJA~gm-AEJ",
            "gm-AEJA=gm-AEJ (4-3)",
            "gm-AEJ=gm-AEJF (3-3)",
            "gm-AQKA=gm-AQK[GGGGG]AA (4-3)",
            "gm-AQK[GGGGG]AA",
            "g(Ac)m(Anh)-AQK[SSA]A=gm-AQK[SSA]AA (4-3)",
            "g",
            "g(Ac)",
            "g(DeAc)",
            "g(Ac, DeAc)",
            "g(Poly)",
            "g(Poly, DeAc)",
            "m(Ac)",
            "m(DeAc)",
            "m(Ac, DeAc)",
            "m(Poly)",
            "m(Poly, DeAc)",
            "m(Glyc)",
            "m(Poly, Glyc)",
            "m(Anh)",
            "m(Anh, DeAc)",
            "D(Am)",
            "E(Am)",
            "J(Am)",
        ];
        let smiles = muropeptides
            .map(|structure| {
                let muropeptide = Muropeptide::new(&POLYMERIZER, structure).unwrap();
                let smiles = db.muropeptide_to_smiles(&muropeptide).unwrap();
                format!("{structure}: {smiles}")
            })
            .join("\n");
        assert_snapshot!(smiles);
    }
}
