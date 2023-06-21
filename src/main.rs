use std::{
    env,
    error::Error,
    fs::{self, File},
    path::Path,
};

use petgraph::dot::Dot;
use polars::prelude::*;
use rayon::prelude::*;
use rust_decimal_macros::dec;
use smithereens::{fragment, fragments_to_df, Peptidoglycan, MODIFICATIONS};
use filenamify::filenamify;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    let structure = &args[1];
    let output = Path::new(&args[2]);

    let pg = Peptidoglycan::new(structure)?;
    println!("Fragmenting Structure: {}", structure);
    println!("Monoisotopic Mass: {}", (pg.monoisotopic_mass()).round_dp(4));
    println!("1+ Ion: {}", (pg.monoisotopic_mass() + MODIFICATIONS["+"].mass).round_dp(4));
    println!("2+ Ion: {}", ((pg.monoisotopic_mass() + dec!(2) * MODIFICATIONS["+"].mass) / dec!(2)).round_dp(4));
    println!("3+ Ion: {}", ((pg.monoisotopic_mass() + dec!(3) * MODIFICATIONS["+"].mass) / dec!(3)).round_dp(4));

    let fragments = fragment(Peptidoglycan::new(structure).unwrap().into());
    println!("Number of Fragments: {}", fragments.len());

    let output = output.join(format!("{}.csv", filenamify(structure)));
    println!("Output: {}", output.display());
    let mut file = File::create(output).unwrap();
    CsvWriter::new(&mut file)
        .finish(&mut fragments_to_df(&fragments))
        .unwrap();
    Ok(())
}
