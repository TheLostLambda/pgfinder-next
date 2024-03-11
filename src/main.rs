use miette::{Diagnostic, GraphicalReportHandler, GraphicalTheme};
use once_cell::sync::Lazy;
use polychem::{
    atoms::atomic_database::AtomicDatabase, Charged, ChemicalComposition, Massive, Mz, Result,
};
use rustyline::DefaultEditor;
use std::fmt::Write;

static DB: Lazy<AtomicDatabase> = Lazy::new(|| {
    AtomicDatabase::from_kdl(
        "atomic_database.kdl",
        include_str!("../polychem/atomic_database.kdl"),
    )
    .unwrap()
});

fn main() {
    let mut rl = DefaultEditor::new().unwrap();
    while let Ok(formula) = rl.readline("Molecule: ") {
        rl.add_history_entry(&formula).unwrap();
        match molecule_info(&formula) {
            Ok(info) => print!("{info}"),
            Err(diagnostic) => render_error(*diagnostic),
        }
    }
}

fn molecule_info(formula: &str) -> Result<String> {
    let mut buf = String::new();
    let molecule = ChemicalComposition::new(&DB, formula)?;

    let mono_mass = molecule.monoisotopic_mass().round_dp(6);
    let avg_mass = molecule.average_mass().round_dp(4);
    let charge = molecule.charge();

    writeln!(buf, "Monoisotopic Mass: {mono_mass}").unwrap();
    writeln!(buf, "Average Mass: {avg_mass}").unwrap();
    writeln!(buf, "Charge: {charge}").unwrap();

    if charge != 0 {
        let mono_mz = molecule.monoisotopic_mz().unwrap().round_dp(6);
        let avg_mz = molecule.average_mz().unwrap().round_dp(4);
        writeln!(buf, "Monoisotopic m/z: {mono_mz}").unwrap();
        writeln!(buf, "Average m/z: {avg_mz}").unwrap();
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
