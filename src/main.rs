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

fn main() -> Result<(), Box<dyn Error>> {
    // let args: Vec<String> = env::args().collect();
    // let input = &args[1];
    // let output = Path::new(&args[2]);
    // let structures = fs::read_to_string(input)?;
    // structures.par_lines().for_each(|structure| {
    //     let fragments = fragment(Peptidoglycan::new(structure).unwrap().into());
    //     let mut file = File::create(output.join(format!("{}.csv", structure))).unwrap();
    //     CsvWriter::new(&mut file)
    //         .finish(&mut fragments_to_df(&fragments))
    //         .unwrap();
    // });

    // let structure = "g(-Ac)m-AE[GA]K[AKEAG]AA";
    // let structure = "m-AK[AA]AA";
    // let structure = "g(+Ac)m-AQK[GGGGSSG]AA";
    // let structure = "gm-AE[G]K[AKEGA]AA";
    // let structure = "gm-AQK[GGGGG]A";
    // let structure = "g(-Ac)m-AEJA";
    // let structure = "gm(-Ac)-AEJA";
    // let structure = "gm-AEJA";
    let structure = "gm-AEJN";
    let pg = dbg!(Peptidoglycan::new(structure)?);
    dbg!((pg.monoisotopic_mass()).round_dp(4));
    dbg!((pg.monoisotopic_mass() + MODIFICATIONS["+"].mass).round_dp(4));
    dbg!(((pg.monoisotopic_mass() + dec!(2) * MODIFICATIONS["+"].mass) / dec!(2)).round_dp(4));
    dbg!(Dot::new(&pg.graph));

    let fragments = fragment(pg.into());
    dbg!(fragments.len());

    let mut file = File::create(format!("/home/tll/Downloads/{}.csv", structure)).unwrap();
    CsvWriter::new(&mut file).finish(&mut dbg!(fragments_to_df(&fragments)))?;
    Ok(())
}
