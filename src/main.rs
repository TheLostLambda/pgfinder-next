use polars::io::{csv::CsvWriter, SerWriter};
// use rayon::prelude::*;
use smithereens::{fragment, fragments_to_df, Peptidoglycan};
use std::{
    fs::{self, File},
    path::PathBuf,
};

// NOTE: To use a crate like `clap`, you just need to run `cargo add clap` to
// automatically add and install the dependency!
use anyhow::Result;
use clap::Parser;

// NOTE: I'm stealing at lot of the `clap` example here!
/// This is a docstring that gets included in the CLI about text!
#[derive(Parser, Debug)]
#[command(version, about)]
struct Args {
    /// A file of newline-separated PG structures
    #[arg(short, long)]
    input_structures: PathBuf,
    /// A directory for outputing fragments / masses
    #[arg(short, long)]
    output_dir: PathBuf,
    /// Whether or not to perform fragmentation
    #[arg(long, default_value_t = false)]
    fragment: bool,
}

// NOTE: For flamegraph profiling, compile with:
// `CARGO_PROFILE_RELEASE_DEBUG=true cargo build --release`
// Then run:
// `flamegraph -- ./target/release/smithereens ...`
// This flamegraph profiling can be performed automatically for benchmarks
// using `cargo-flamegraph` — this is just a toy example!
fn main() -> Result<()> {
    // NOTE: Try adding `dbg!()` around the parsing here!
    let args = Args::parse();

    let structures = fs::read_to_string(args.input_structures)?;
    // NOTE: Change `.lines()` to `.lines_par()` for free data-parallelism
    // Or at least when it doesn't overflow the stack ;-;
    let monoisotopic_masses: Vec<_> = structures
        .lines()
        .map(|structure| {
            let pg = Peptidoglycan::new(structure);
            if args.fragment {
                let fragments = fragment(pg.clone().into());
                let output_file = args.output_dir.join(format!("{}.csv", structure));
                let mut file = File::create(output_file).unwrap();
                CsvWriter::new(&mut file)
                    .finish(&mut fragments_to_df(&fragments))
                    .unwrap();
            }
            format!("{},{}", structure, pg.monoisotopic_mass())
        })
        .collect();

    let massfile_path = args.output_dir.join("masses.txt");
    fs::write(massfile_path, monoisotopic_masses.join("\n"))?;

    Ok(())
}
