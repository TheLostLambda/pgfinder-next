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
    loop {
        if let Ok(formula) = rl.readline("Molecule: ") {
            rl.add_history_entry(&formula).unwrap();
            match molecule_info(&formula) {
                Ok(info) => print!("{info}"),
                Err(diagnostic) => render_error(diagnostic),
            }
        } else {
            break;
        }
    }
}

fn molecule_info(formula: &str) -> Result<String, PolychemError> {
    let mut buf = String::new();
    let molecule = ChemicalComposition::new(&DB, formula)?;

    writeln!(buf, "Monoisotopic Mass: {}", molecule.monoisotopic_mass()?).unwrap();
    writeln!(buf, "Average Mass: {}", molecule.average_mass()?).unwrap();
    // TODO: Add charge!
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
