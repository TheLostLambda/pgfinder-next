use divan::{AllocProfiler, black_box};
use polychem::{
    AtomicDatabase, Charged, ChemicalComposition, Massive, PolymerDatabase, atoms::atomic_database,
};
use std::sync::LazyLock;

#[global_allocator]
static ALLOC: AllocProfiler = AllocProfiler::system();

const ATOMIC_KDL: &str = atomic_database::DEFAULT_KDL;
const POLYMER_KDL: &str = include_str!("../tests/data/polymer_database.kdl");
const FORMULAS: [&str; 5] = [
    "C2H5NO2",
    "C5H9NO2",
    "C7H14N2O4",
    "C5H11NO2S",
    "C3H7[15N]O2",
];

static ATOMIC_DB: LazyLock<AtomicDatabase> = LazyLock::new(AtomicDatabase::default);

static POLYMER_DB: LazyLock<PolymerDatabase> = LazyLock::new(|| {
    PolymerDatabase::new(&ATOMIC_DB, "test_polymer_database.kdl", POLYMER_KDL).unwrap()
});

static COMPOSITIONS: LazyLock<Vec<ChemicalComposition>> = LazyLock::new(|| {
    FORMULAS
        .into_iter()
        .map(|formula| ChemicalComposition::new(&ATOMIC_DB, formula).unwrap())
        .collect()
});

fn main() {
    LazyLock::force(&ATOMIC_DB);
    LazyLock::force(&POLYMER_DB);
    LazyLock::force(&COMPOSITIONS);
    divan::main();
}

mod atoms {
    use super::*;

    #[divan::bench]
    fn build_atomic_database() -> AtomicDatabase {
        AtomicDatabase::new("atomic_database.kdl", ATOMIC_KDL).unwrap()
    }

    #[divan::bench]
    fn parse_chemical_compositions() {
        for formula in FORMULAS.into_iter() {
            black_box(ChemicalComposition::new(&ATOMIC_DB, formula).unwrap());
        }
    }

    #[divan::bench]
    fn calculate_monoisotopic_masses() {
        for composition in COMPOSITIONS.iter() {
            black_box(composition.monoisotopic_mass());
        }
    }

    #[divan::bench]
    fn calculate_average_masses() {
        for composition in COMPOSITIONS.iter() {
            black_box(composition.average_mass());
        }
    }

    #[divan::bench]
    fn calculate_charges() {
        for composition in COMPOSITIONS.iter() {
            black_box(composition.charge());
        }
    }
}

mod polymers {
    use super::*;

    #[divan::bench]
    fn build_polymer_database() -> PolymerDatabase<'static> {
        PolymerDatabase::new(&ATOMIC_DB, "test_polymer_database.kdl", POLYMER_KDL).unwrap()
    }

    #[ignore]
    #[divan::bench]
    fn parse_residues() {
        // TODO: Benchmark residue parsing / creation via `Polymer`
        todo!()
    }
}

mod polymerizer {
    // TODO: Benchmark some of the public api!
}
