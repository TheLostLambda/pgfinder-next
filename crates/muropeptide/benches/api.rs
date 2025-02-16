use divan::AllocProfiler;
use polychem::{AtomicDatabase, PolymerDatabase, Polymerizer, atoms::atomic_database};
use std::sync::LazyLock;

#[global_allocator]
static ALLOC: AllocProfiler = AllocProfiler::system();

const ATOMIC_KDL: &str = atomic_database::DEFAULT_KDL;
const POLYMER_KDL: &str = include_str!("../data/polymer_database.kdl");
const AMINO_ACIDS: [&str; 3] = ["D(Am)", "E(Am)", "J(Am)"];
const MONOMERS: [&str; 6] = [
    "m",
    "gmgmgmgmgm",
    "A",
    "AEJAAEJAAEJAAEJAAEJA",
    "m-A",
    "gmgmgmgmgm-AEJAAEJAAEJAAEJAAEJA",
];

static ATOMIC_DB: LazyLock<AtomicDatabase> = LazyLock::new(AtomicDatabase::default);
static POLYMER_DB: LazyLock<PolymerDatabase> = LazyLock::new(|| {
    PolymerDatabase::new(&ATOMIC_DB, "polymer_database.kdl", POLYMER_KDL).unwrap()
});

static POLYMERIZER: LazyLock<Polymerizer> =
    LazyLock::new(|| Polymerizer::new(&ATOMIC_DB, &POLYMER_DB));

fn main() {
    LazyLock::force(&ATOMIC_DB);
    LazyLock::force(&POLYMER_DB);
    LazyLock::force(&POLYMERIZER);
    divan::main();
}

mod polymers {
    use muropeptide::Muropeptide;

    // TODO: Add some benchmarks here!
    // FIXME: Very messy and temporary!
    use super::*;

    #[divan::bench]
    fn build_muropeptides<'a, 'p>() -> Vec<Muropeptide<'a, 'p>> {
        MONOMERS
            .iter()
            .map(|structure| Muropeptide::new(&POLYMERIZER, structure).unwrap())
            .collect()
    }
}
