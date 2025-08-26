use std::{fs::File, io};

use mzdata::{prelude::*, spectrum::SpectrumDescription};
use polars::prelude::*;

// REMEMBER: I'm prototyping! Code will be re-written afterwards!

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

#[derive(Debug)]
struct Ms2Scan {
    index: usize,
    start_time: f64,
    precursor_mz: f64,
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
    let ms2_scans: Vec<_> = mzdata::mz_read!(
        path.as_ref(),
        reader => {
            reader.filter_map(fetch_ms2_precursor).collect()
    })?;

    let start_index = 2280;
    let length = 5;
    dbg!(&ms2_scans[start_index..start_index + length]);

    Ok(())
}
