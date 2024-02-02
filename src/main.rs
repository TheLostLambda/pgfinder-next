use miette::{Diagnostic, GraphicalReportHandler, GraphicalTheme};
use once_cell::sync::Lazy;
use polychem::{chemical_database::ChemicalDatabase, ChemicalComposition, PolychemError};
use rustyline::DefaultEditor;
use std::fmt::Write;

static DB: Lazy<ChemicalDatabase> = Lazy::new(|| {
    ChemicalDatabase::from_kdl("chemistry.kdl", include_str!("../polychem/chemistry.kdl")).unwrap()
});

fn main() {
    let mut rl = DefaultEditor::new().unwrap();
    while let Ok(formula) = rl.readline("Molecule: ") {
        rl.add_history_entry(&formula).unwrap();
        match molecule_info(&formula) {
            Ok(info) => print!("{info}"),
            Err(diagnostic) => render_error(diagnostic),
        }
    }
}

fn molecule_info(formula: &str) -> Result<String, PolychemError> {
    let mut buf = String::new();
    let molecule = ChemicalComposition::new(&DB, formula)?;

    let mono_mass = molecule.monoisotopic_mass()?.round_dp(6);
    let avg_mass = molecule.average_mass()?.round_dp(6);
    let charge = molecule.charge();

    writeln!(buf, "Monoisotopic Mass: {mono_mass}").unwrap();
    writeln!(buf, "Average Mass: {avg_mass}").unwrap();
    writeln!(buf, "Charge: {charge}").unwrap();
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
