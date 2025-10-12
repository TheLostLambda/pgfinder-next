use std::{cmp::Ordering, io::Cursor, ops::RangeBounds};

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
    fn ppm_window(self, ppm: f64) -> impl RangeBounds<Self> {
        let (min_mz, max_mz) = Tolerance::PPM(ppm).bounds(self.0);
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

#[cfg(test)]
mod tests {
    use std::ops::Bound;

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
        let window = Mz::from(941.407_703).ppm_window(7.5);
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
}
