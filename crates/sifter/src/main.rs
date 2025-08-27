// NOTE: Use `ThermoRawFileParser -i 6ldt_R.raw -f 1 -g`
// NOTE: Run `cargo run --release o+e>| helix`

use std::{cmp::Ordering, fs::File, io};

use mzdata::{prelude::*, spectrum::SpectrumDescription};
use polars::prelude::*;

// REMEMBER: I'm prototyping! Code will be re-written afterwards!

// FIXME: This should really come from `polychem`'s atomic database, and absolutely shouldn't be copied here!
const PROTON_MASS: f64 = 1.007_276_466_621;
// const MZML: &str = "/home/tll/Downloads/MS DATA/6ldt_R.mzML.gz";
const MZML: &str = "tests/data/WT.mzML.gz";
const PGF_OUTPUT: &str = "tests/data/WT.csv";
const MASS_DB: &str = "tests/data/Monomers.csv";

// FIXME: Unify these two structs with a common trait? Or make them perfect hashmaps or something?
struct PGFInputColumns;
impl PGFInputColumns {
    const ALL: [&str; 3] = [Self::STRUCTURE, Self::RT, Self::THEO];
    const STRUCTURE: &str = "Structure";
    const RT: &str = "Consolidated RT (min)";
    const THEO: &str = "Consolidated Theo (Da)";
}

struct DBInputColumns;
impl DBInputColumns {
    const STRUCTURE: &str = "Structure";
    const THEO: &str = "Monoisotopic Mass";
}

struct TemporaryColumns;
impl TemporaryColumns {
    const ION_MZS: &str = "Ion M/Zs";
    const ION_CHARGES: &str = "Ion Charges";
}

fn main() {
    println!("Reading {MZML}...");
    let scans = mz_read_path(MZML).unwrap();
    // println!("Reading {PGF_OUTPUT}...");
    // read_pgfinder_output(PGF_OUTPUT);
    println!("Reading {MASS_DB}...");
    let df = read_db(MASS_DB);
    let found_ms2_scans = dbg!(
        df[DBInputColumns::STRUCTURE]
            .str()
            .unwrap()
            .iter()
            .map(Option::unwrap)
            .zip(
                df[DBInputColumns::THEO]
                    .f64()
                    .unwrap()
                    .iter()
                    .map(Option::unwrap)
            )
            .map(|(name, mass)| (
                name,
                mass,
                ion_series(mass, 150.)
                    .iter()
                    .map(|&mz| (mz, within_ppm(&scans, mz, 10.)))
                    .collect::<Vec<_>>()
            ))
            .collect::<Vec<_>>()
    );

    let (with_data, without_data) = found_ms2_scans
        .iter()
        .map(|&(name, mass, ref ions)| {
            (
                name,
                mass,
                ions.iter().map(|(_, scans)| scans.len()).sum::<usize>(),
            )
        })
        .partition::<Vec<_>, _>(|&(_, _, ions)| ions > 0);

    println!("\n\nStructures with MS2 data ({}):", with_data.len());
    for (name, _, scans) in with_data {
        println!("  - {name} :: {scans}");
    }

    println!("\nStructures without MS2 data ({})", without_data.len());
}

// =====================================================================================================================

// FIXME: Replace &str path with a `impl Read + Seek`?
fn read_pgfinder_output(path: &str) {
    let df = CsvReader::new(File::open(path).unwrap())
        .finish()
        // FIXME: Handle this properly!
        .unwrap();
    dbg!(
        df.lazy()
            .select(PGFInputColumns::ALL.map(col))
            .drop_nulls(None)
            .collect()
            .unwrap()
    );
}

// FIXME: Not DRY!
fn read_db(path: &str) -> DataFrame {
    CsvReader::new(File::open(path).unwrap())
        .finish()
        // FIXME: Handle this properly!
        .unwrap()
    // dbg!(
    //     df.lazy()
    //         .with_columns_seq([col(DBInputColumns::THEO)
    //             .map(
    //                 |column| Ok(Some(
    //                     Series::new(
    //                         TemporaryColumns::ION_MZS.into(),
    //                         column
    //                             .as_materialized_series()
    //                             .f64()
    //                             .unwrap()
    //                             .iter()
    //                             .map(|mass| ion_series(mass.unwrap(), 150.)
    //                                 .iter()
    //                                 .collect::<Series>()
    //                                 .implode()
    //                                 .unwrap())
    //                             .collect::<Vec<_>>()
    //                     )
    //                     .into()
    //                 )),
    //                 GetOutput::from_type(DataType::List(Box::new(DataType::Float64)))
    //             )
    //             .alias(TemporaryColumns::ION_MZS)])
    //         .collect()
    //         .unwrap()
    // );
}

#[derive(Debug, PartialEq)]
struct Ms2Scan {
    index: usize,
    start_time: f64,
    precursor_mz: f64,
}

impl Ms2Scan {
    fn cmp(a: &Self, b: &Self) -> Ordering {
        a.precursor_mz.total_cmp(&b.precursor_mz)
    }
}

// NOTE: It's a private function, and there's less boilerplate when taking by value
#[expect(clippy::needless_pass_by_value)]
fn fetch_ms2_precursor(spectrum: impl SpectrumLike) -> Option<Ms2Scan> {
    let SpectrumDescription {
        index,
        ms_level,
        acquisition,
        precursor,
        ..
    } = spectrum.description();

    (*ms_level == 2).then(|| {
        let precursor_ions = precursor
            .as_ref()
            .expect("MS2 scan was missing precursor ion information")
            .ions
            .as_slice();
        assert_eq!(
            precursor_ions.len(),
            1,
            "multiple precursor ions are not yet supported"
        );

        Ms2Scan {
            index: index + 1,
            start_time: acquisition.start_time(),
            precursor_mz: precursor_ions[0].mz,
        }
    })
}

// FIXME: I should use the writing functionality of this crate to generate some smaller test files!!!
// FIXME: Replace &str path with a `impl Read + Seek`?
fn mz_read_path(path: &str) -> io::Result<Vec<Ms2Scan>> {
    // FIXME: `mz_read!` should probably add that `.as_ref()` itself...
    let mut ms2_scans: Vec<_> = mzdata::mz_read!(
        path.as_ref(),
        reader => {
            reader.filter_map(fetch_ms2_precursor).collect()
    })?;

    ms2_scans.sort_unstable_by(Ms2Scan::cmp);

    // dbg!(within_ppm(&ms2_scans, 941.4077 + PROTON_MASS, 10.));
    // let start_index = 2280;
    // let length = 3;
    // dbg!(&ms2_scans[start_index..start_index + length]);

    Ok(ms2_scans)
}

// FIXME: Should maybe move from f64 to Decimal?
fn ppm_window(theoretical_mz: f64, ppm: f64) -> (f64, f64) {
    let half_window = ppm * theoretical_mz / 1_000_000.;
    (theoretical_mz - half_window, theoretical_mz + half_window)
}

fn within_ppm(scans: &[Ms2Scan], theoretical_mz: f64, ppm: f64) -> &[Ms2Scan] {
    assert!(scans.is_sorted_by(|a, b| Ms2Scan::cmp(a, b) != Ordering::Greater));
    let (min_mz, max_mz) = ppm_window(theoretical_mz, ppm);
    let start = scans.partition_point(|s| s.precursor_mz < min_mz);
    let end = scans.partition_point(|s| s.precursor_mz <= max_mz);

    &scans[start..end]
}

fn ion_series(theoretical_mass: f64, cutoff_mz: f64) -> Vec<f64> {
    #[expect(clippy::maybe_infinite_iter)]
    (1..)
        .map(|charge| {
            let charge = f64::from(charge);
            charge.mul_add(PROTON_MASS, theoretical_mass) / charge
        })
        .take_while(|&mz| mz > cutoff_mz)
        .collect()
}
