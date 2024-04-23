use divan::{black_box, AllocProfiler};
use once_cell::sync::Lazy;
use polychem::{
    atoms::atomic_database, AtomicDatabase, Charged, ChemicalComposition, Massive, PolymerDatabase,
};

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

static ATOMIC_DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);

static POLYMER_DB: Lazy<PolymerDatabase> =
    Lazy::new(|| PolymerDatabase::new(&ATOMIC_DB, "polymer_database.kdl", POLYMER_KDL).unwrap());

static COMPOSITIONS: Lazy<Vec<ChemicalComposition>> = Lazy::new(|| {
    FORMULAS
        .into_iter()
        .map(|formula| ChemicalComposition::new(&ATOMIC_DB, formula).unwrap())
        .collect()
});

fn main() {
    Lazy::force(&ATOMIC_DB);
    Lazy::force(&POLYMER_DB);
    Lazy::force(&COMPOSITIONS);
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
        PolymerDatabase::new(&ATOMIC_DB, "polymer_database.kdl", POLYMER_KDL).unwrap()
    }

    #[divan::bench]
    fn parse_residues() {
        // TODO: Benchmark residue parsing / creation via `Polymer`
        todo!()
    }
}

mod polymerizer {
    // TODO: Benchmark some of the public api!
}
