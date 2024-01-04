use std::{error::Error, fs::File};

use petgraph::dot::Dot;
use polars::prelude::*;
use rust_decimal_macros::dec;
use smithereens::{fragment, fragments_to_df, Peptidoglycan, MODIFICATIONS};

fn main() -> Result<(), Box<dyn Error>> {
    // let structure = "g(-Ac)m-AE[GA]K[AKEAG]AA";
    // let structure = "m-AK[AA]AA";
    // let structure = "g(+Ac)m-AQK[GGGGSSG]AA";
    // let structure = "gm-AE[G]K[AKEGA]AA";
    let structure = "gm-AQK[GGGGG]A";
    // let structure = "g(-Ac)m-AEJA";
    // let structure = "gm(-Ac)-AEJA";
    // let structure = "gm-AEJA";
    // let structure = "gmgm(-Ac)-AEJA";
    let pg = dbg!(Peptidoglycan::new(structure));
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
