use miette::{Diagnostic, GraphicalReportHandler, GraphicalTheme};
use polychem::{AtomicDatabase, Charged, ChargedParticle, ChemicalComposition, Massive, Result};
use rust_decimal::Decimal;
use rustyline::DefaultEditor;
use std::{fmt::Write, sync::LazyLock};

static DB: LazyLock<AtomicDatabase> = LazyLock::new(AtomicDatabase::default);

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

    let mono_mass = molecule.monoisotopic_mass();
    let avg_mass = molecule.average_mass();
    let charge = molecule.charge();

    writeln!(
        buf,
        "Monoisotopic Mass: {}",
        decimal_round_workaround(mono_mass, 6)
    )
    .unwrap();
    writeln!(
        buf,
        "Average Mass: {}",
        decimal_round_workaround(avg_mass, 4)
    )
    .unwrap();
    writeln!(buf, "Charge: {charge}").unwrap();

    if i64::from(charge) != 0 {
        let mono_mz = molecule.monoisotopic_mz().unwrap();
        let avg_mz = molecule.average_mz().unwrap();
        writeln!(
            buf,
            "Monoisotopic m/z: {}",
            decimal_round_workaround(mono_mz, 6)
        )
        .unwrap();
        writeln!(buf, "Average m/z: {}", decimal_round_workaround(avg_mz, 4)).unwrap();
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

// FIXME: Really this should be fixed in `rust_decimal`...
fn decimal_round_workaround(value: impl Into<Decimal>, decimal_points: u32) -> String {
    let value = value.into().round_dp(decimal_points);
    format!("{value}")
}
