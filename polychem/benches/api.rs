use divan::{black_box, AllocProfiler};
use once_cell::sync::Lazy;
use polychem::{AtomicDatabase, Charged, ChemicalComposition, Massive, PolymerDatabase};

#[global_allocator]
static ALLOC: AllocProfiler = AllocProfiler::system();

const ATOMIC_KDL: &str = include_str!("../atomic_database.kdl");
const POLYMER_KDL: &str = include_str!("../muropeptide_chemistry.kdl");
const FORMULAS: [&str; 5] = [
    "C2H5NO2",
    "C5H9NO2",
    "C7H14N2O4",
    "C5H11NO2S",
    "C3H7[15N]O2",
];

static ATOMIC_DB: Lazy<AtomicDatabase> =
    Lazy::new(|| AtomicDatabase::from_kdl("atomic_database.kdl", ATOMIC_KDL).unwrap());

static POLYMER_DB: Lazy<PolymerDatabase> = Lazy::new(|| {
    PolymerDatabase::from_kdl(&ATOMIC_DB, "muropeptide_chemistry.kdl", POLYMER_KDL).unwrap()
});

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
        AtomicDatabase::from_kdl("atomic_database.kdl", ATOMIC_KDL).unwrap()
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
    use divan::Bencher;
    use polychem::Residue;

    use super::*;

    #[divan::bench]
    fn build_polymer_database() -> PolymerDatabase<'static> {
        PolymerDatabase::from_kdl(&ATOMIC_DB, "muropeptide_chemistry.kdl", POLYMER_KDL).unwrap()
    }

    #[divan::bench]
    fn parse_residues(bencher: Bencher) {
        let abbrs: Vec<_> = ('A'..'Z').map(|c| c.to_string()).collect();
        bencher.bench_local(|| {
            for abbr in &abbrs {
                black_box(Residue::new(&POLYMER_DB, abbr, 0).unwrap());
            }
        });
    }
}
