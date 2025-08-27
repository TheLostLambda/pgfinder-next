use std::{cmp::Ordering, fs::File, io};

use mzdata::{prelude::*, spectrum::SpectrumDescription};
use polars::prelude::*;

// REMEMBER: I'm prototyping! Code will be re-written afterwards!

// FIXME: This should really come from `polychem`'s atomic database, and absolutely shouldn't be copied here!
const PROTON_MASS: f64 = 1.007_276_466_621;
// const MZML: &str = "/home/tll/Downloads/MS DATA/6ldt_R.mzML.gz";
const MZML: &str = "tests/data/WT.mzML.gz";
const PGF_OUTPUT: &str = "tests/data/WT.csv";

struct InputColumns;
impl InputColumns {
    const ALL: [&str; 3] = [Self::STRUCTURE, Self::RT, Self::THEO];
    const STRUCTURE: &str = "Structure";
    const RT: &str = "Consolidated RT (min)";
    const THEO: &str = "Consolidated Theo (Da)";
}

fn main() {
    println!("Reading {MZML}...");
    mz_read_path(MZML).unwrap();
    println!("Reading {PGF_OUTPUT}...");
    read_pgfinder_output(PGF_OUTPUT);
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
            .select(InputColumns::ALL.map(col))
            .drop_nulls(None)
            .collect()
            .unwrap()
    );
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
fn mz_read_path(path: &str) -> io::Result<()> {
    // FIXME: `mz_read!` should probably add that `.as_ref()` itself...
    let mut ms2_scans: Vec<_> = mzdata::mz_read!(
        path.as_ref(),
        reader => {
            reader.filter_map(fetch_ms2_precursor).collect()
    })?;

    ms2_scans.sort_unstable_by(Ms2Scan::cmp);

    dbg!(within_ppm(&ms2_scans, dbg!(941.4077 + PROTON_MASS), 10.));
    // let start_index = 2280;
    // let length = 3;
    // dbg!(&ms2_scans[start_index..start_index + length]);

    Ok(())
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
