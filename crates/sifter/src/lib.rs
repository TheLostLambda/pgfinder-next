use std::{cmp::Ordering, collections::BTreeSet, io::Cursor, ops::RangeInclusive};

use derive_more::{From, Into};
use mzdata::mzpeaks::Tolerance;
use polars::{error::PolarsResult, frame::DataFrame, io::SerReader, prelude::CsvReader};

fn load_mass_database(csv: &str) -> PolarsResult<DataFrame> {
    let csv_reader = Cursor::new(csv);
    CsvReader::new(csv_reader).finish()
}

// PERF: Try out `f32` and see if that saves enough space to speed things up!
#[derive(Copy, Clone, Debug, From, Into)]
struct Mz(f64);

impl Mz {
    fn ppm_window(mz: f64, ppm: f64) -> RangeInclusive<Self> {
        let (min_mz, max_mz) = Tolerance::PPM(ppm).bounds(mz);
        Self(min_mz)..=Self(max_mz)
    }
}

impl Ord for Mz {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

impl PartialOrd for Mz {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Mz {
    fn eq(&self, other: &Self) -> bool {
        // NOTE: We must use the `.cmp()` from `Ord` here — the results of `.total_cmp()` do *not* match the default
        // implementations of `PartialEq` for floats. The `.total_cmp()` method treats `-0` and `0` as different, but
        // the default `PartialEq` for floats does not!
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Mz {}

// PERF: After getting some benchmarks in place, try cutting the size of this struct in half (f64 → f32) — the 64 bits
// definitely isn't needed precision-wise, so if it helps performance, it's probably worth the small amount of type
// casting! Don't forget about padding! All fields must shrink to the same size.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
struct ScanKey {
    precursor: Mz,
    scan_number: usize,
}

impl ScanKey {
    // MISSING: It's not public API, so `const` is unnecessary here
    #[expect(clippy::missing_const_for_fn)]
    fn new(precursor: f64, scan_number: usize) -> Self {
        Self {
            precursor: Mz(precursor),
            scan_number,
        }
    }

    fn ppm_window(mz: f64, ppm: f64) -> RangeInclusive<Self> {
        let (min_mz, max_mz) = Tolerance::PPM(ppm).bounds(mz);
        Self::new(min_mz, 0)..=Self::new(max_mz, usize::MAX)
    }
}

#[derive(Debug)]
struct Peaks(BTreeSet<Mz>);

impl Peaks {
    fn filter_peaks(&self, mz: f64, ppm: f64) -> impl Iterator<Item = f64> {
        self.0.range(Mz::ppm_window(mz, ppm)).map(|mz| mz.0)
    }
}

#[cfg(test)]
mod tests {
    use std::{
        collections::BTreeSet,
        ops::{Bound, RangeBounds},
    };

    use assert_float_eq::assert_float_absolute_eq;
    use polars::prelude::DataType;

    use super::*;

    const MASS_DATABASE_CSV: &str = include_str!("../tests/data/Exchange.csv");

    #[test]
    fn test_load_mass_database() {
        let mass_database = load_mass_database(MASS_DATABASE_CSV).unwrap();
        assert_eq!(
            mass_database.get_column_names(),
            &["Structure", "Monoisotopic Mass"]
        );
        assert_eq!(
            mass_database.dtypes(),
            &[DataType::String, DataType::Float64]
        );
        assert_eq!(mass_database.height(), 21);
    }

    #[test]
    fn mz_traits() {
        // Test `From` and `Into` impls
        let a: Mz = (941.407_703).into();
        let b = Mz::from(471.711_128);
        assert_float_absolute_eq!(a.into(), 941.407_703);
        assert_float_absolute_eq!(b.into(), 471.711_128);
        // Test `Copy` and `Ord` impl
        assert_eq!(a.cmp(&b), Ordering::Greater);
        assert_eq!(b.cmp(&a), Ordering::Less);
        assert_eq!(a.cmp(&a), Ordering::Equal);
        assert_eq!(b.cmp(&b), Ordering::Equal);
    }

    #[test]
    fn mz_ppm_window() {
        let window = Mz::ppm_window(941.407_703, 7.5);
        let Bound::Included(&Mz(start)) = window.start_bound() else {
            panic!("expected a `.start_bound()` matching `Bound::Included(&Mz(_))`");
        };
        let Bound::Included(&Mz(end)) = window.end_bound() else {
            panic!("expected an `.end_bound()` matching `Bound::Included(&Mz(_))`");
        };
        assert_float_absolute_eq!(start, 941.400_642);
        assert_float_absolute_eq!(end, 941.414_764);
        assert!(window.contains(&Mz(941.407_703)));
    }

    #[test]
    fn scan_key_ppm_window() {
        let window = ScanKey::ppm_window(471.711_128, 10.0);
        let Bound::Included(&ScanKey {
            precursor: Mz(start_mz),
            scan_number: start_scan,
        }) = window.start_bound()
        else {
            panic!("expected a `.start_bound()` matching `Bound::Included(&Mz(_))`");
        };
        let Bound::Included(&ScanKey {
            precursor: Mz(end_mz),
            scan_number: end_scan,
        }) = window.end_bound()
        else {
            panic!("expected an `.end_bound()` matching `Bound::Included(&Mz(_))`");
        };
        assert_float_absolute_eq!(start_mz, 471.706_411);
        assert_eq!(start_scan, 0);
        assert_float_absolute_eq!(end_mz, 471.715_845);
        assert_eq!(end_scan, usize::MAX);
        assert!(window.contains(&ScanKey::new(471.711_128, 42)));
    }

    #[test]
    fn peaks_filter_peaks() {
        let query = 942.412_793_603_516;
        let peaks = Peaks(BTreeSet::from(
            [
                942.413_696_289_063,
                942.4121,
                942.413_452_148_438,
                942.4121,
                942.4119,
                942.413_513_183_594,
            ]
            .map(Mz),
        ));
        let filter_peaks = |ppm| -> Vec<_> { peaks.filter_peaks(query, ppm).collect() };
        assert_eq!(filter_peaks(0.0), vec![]);
        assert_eq!(filter_peaks(0.75), vec![942.4121, 942.413_452_148_438]);
        assert_eq!(
            filter_peaks(10.0),
            vec![
                942.4119,
                942.4121,
                942.413_452_148_438,
                942.413_513_183_594,
                942.413_696_289_063
            ]
        );
    }
}
