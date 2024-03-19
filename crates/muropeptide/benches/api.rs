use divan::{black_box, AllocProfiler};
use once_cell::sync::Lazy;
use polychem::{AtomicDatabase, PolymerDatabase};

#[global_allocator]
static ALLOC: AllocProfiler = AllocProfiler::system();

const ATOMIC_KDL: &str = include_str!("../../polychem/data/atomic_database.kdl");
const POLYMER_KDL: &str = include_str!("../data/polymer_database.kdl");
const AMINO_ACIDS: [&str; 3] = ["D(Am)", "E(Am)", "J(Am)"];

static ATOMIC_DB: Lazy<AtomicDatabase> =
    Lazy::new(|| AtomicDatabase::new("atomic_database.kdl", ATOMIC_KDL).unwrap());

static POLYMER_DB: Lazy<PolymerDatabase> =
    Lazy::new(|| PolymerDatabase::new(&ATOMIC_DB, "polymer_database.kdl", POLYMER_KDL).unwrap());

fn main() {
    Lazy::force(&ATOMIC_DB);
    Lazy::force(&POLYMER_DB);
    divan::main();
}

mod polymers {
    use divan::Bencher;
    use muropeptide::parser::unbranched_amino_acid;
    use polychem::polymerizer::Polymerizer;

    use super::*;

    // FIXME: Will probably blow this away...
    #[divan::bench(sample_count = 100, sample_size = 10)]
    fn residues_with_modifications(bencher: Bencher) {
        let amino_acids = AMINO_ACIDS.repeat(1000);
        let mut polymerizer = Polymerizer::new(&ATOMIC_DB, &POLYMER_DB);
        let mut residue_builder = unbranched_amino_acid(&mut polymerizer);
        bencher.bench_local(move || {
            for abbr in &amino_acids {
                black_box(residue_builder(abbr).unwrap());
            }
        });
    }
}
