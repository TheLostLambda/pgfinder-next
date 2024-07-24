use once_cell::sync::Lazy;
use polychem::{
    moieties::polymer_database::{ModificationDescription, ResidueDescription},
    AtomicDatabase, Massive, PolymerDatabase,
};

static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);
static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
    PolymerDatabase::new(
        &ATOMIC_DB,
        "polymer_database.kdl",
        include_str!("../../crates/muropeptide/data/polymer_database.kdl"),
    )
    .unwrap()
});

fn main() {
    let mut residues: Vec<_> = POLYMER_DB.residues.iter().map(write_residue).collect();
    residues.sort_unstable();
    let residues = residues.join("\n");
    println!("{residues}");

    println!();

    let mut modifications: Vec<_> = POLYMER_DB
        .modifications
        .iter()
        .map(write_modification)
        .collect();
    modifications.sort_unstable();
    let modifications = modifications.join("\n");
    println!("{modifications}");
}

fn write_residue((abbr, description): (&String, &ResidueDescription)) -> String {
    let ResidueDescription {
        name, composition, ..
    } = description;

    let formula = composition.to_string();
    let monoisotopic_mass = composition.monoisotopic_mass();

    format!(r#"{abbr},"{name}",{formula},{monoisotopic_mass:.6}"#)
}

fn write_modification((abbr, description): (&String, &ModificationDescription)) -> String {
    let ModificationDescription {
        name,
        lost,
        gained,
        targets,
    } = description;

    let lost_formula = lost.to_string();
    let gained_formula = gained.to_string();
    let monoisotopic_mass = gained.monoisotopic_mass() - lost.monoisotopic_mass();
    let targets: Vec<_> = targets
        .iter()
        .map(ToString::to_string)
        .map(|s| s.replace('"', "\"\""))
        .collect();
    let targets = targets.join("\n");

    format!(r#"{abbr},"{name}",{lost_formula},{gained_formula},{monoisotopic_mass:.6},"{targets}""#)
}
