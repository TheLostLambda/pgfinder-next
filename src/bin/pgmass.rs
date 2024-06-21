use miette::{Diagnostic, GraphicalReportHandler, GraphicalTheme};
use muropeptide::{Muropeptide, Result};
use once_cell::sync::Lazy;
use polychem::{AtomicDatabase, Charged, ChargedParticle, Massive, PolymerDatabase, Polymerizer};
use rustyline::DefaultEditor;
use std::fmt::Write;

static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);
static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
    PolymerDatabase::new(
        &ATOMIC_DB,
        "polymer_database.kdl",
        include_str!("../../crates/muropeptide/data/polymer_database.kdl"),
    )
    .unwrap()
});

static POLYMERIZER: Lazy<Polymerizer> = Lazy::new(|| Polymerizer::new(&ATOMIC_DB, &POLYMER_DB));

// FIXME: This entire binary is just copy-pasted, but needs to be rewritten for PG!
fn main() {
    let mut rl = DefaultEditor::new().unwrap();
    while let Ok(formula) = rl.readline("Peptidoglycan: ") {
        rl.add_history_entry(&formula).unwrap();
        match pg_info(&formula) {
            Ok(info) => print!("{info}"),
            Err(diagnostic) => render_error(dbg!(diagnostic)),
        }
    }
}

fn pg_info(formula: &str) -> Result<String> {
    let mut buf = String::new();
    let muropeptide = Muropeptide::new(&POLYMERIZER, formula)?;

    let mono_mass = muropeptide.monoisotopic_mass();
    let avg_mass = muropeptide.average_mass();
    let charge = muropeptide.charge();

    writeln!(buf, "Monoisotopic Mass: {mono_mass:.6}").unwrap();
    writeln!(buf, "Average Mass: {avg_mass:.4}").unwrap();
    writeln!(buf, "Charge: {charge}").unwrap();

    if i64::from(charge) != 0 {
        let mono_mz = muropeptide.monoisotopic_mz().unwrap();
        let avg_mz = muropeptide.average_mz().unwrap();
        writeln!(buf, "Monoisotopic m/z: {mono_mz:.6}").unwrap();
        writeln!(buf, "Average m/z: {avg_mz:.4}").unwrap();
    }

    writeln!(buf).unwrap();

    Ok(buf)
}

fn render_error(diagnostic: impl Into<Box<dyn Diagnostic + 'static>>) {
    let mut buf = String::new();
    GraphicalReportHandler::new_themed(GraphicalTheme::unicode())
        .render_report(&mut buf, diagnostic.into().as_ref())
        .unwrap();
    println!("{buf}");
}
