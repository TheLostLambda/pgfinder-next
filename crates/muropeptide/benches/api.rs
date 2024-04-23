use divan::AllocProfiler;
use once_cell::sync::Lazy;
use polychem::{atoms::atomic_database, AtomicDatabase, PolymerDatabase};

#[global_allocator]
static ALLOC: AllocProfiler = AllocProfiler::system();

const ATOMIC_KDL: &str = atomic_database::DEFAULT_KDL;
const POLYMER_KDL: &str = include_str!("../data/polymer_database.kdl");
const AMINO_ACIDS: [&str; 3] = ["D(Am)", "E(Am)", "J(Am)"];

static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);

static POLYMER_DB: Lazy<PolymerDatabase> =
    Lazy::new(|| PolymerDatabase::new(&ATOMIC_DB, "polymer_database.kdl", POLYMER_KDL).unwrap());

fn main() {
    Lazy::force(&ATOMIC_DB);
    Lazy::force(&POLYMER_DB);
    divan::main();
}

mod polymers {
    // TODO: Add some benchmarks here!
}
