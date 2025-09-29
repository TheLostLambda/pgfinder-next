use std::io::Cursor;

use polars::{error::PolarsResult, frame::DataFrame, io::SerReader, prelude::CsvReader};

fn load_mass_database(csv: &str) -> PolarsResult<DataFrame> {
    let csv_reader = Cursor::new(csv);
    CsvReader::new(csv_reader).finish()
}

fn ppm_window(mz: f64, ppm: f64) -> (f64, f64) {
    let half_window = mz * ppm / 1_000_000.;
    (mz - half_window, mz + half_window)
}

#[cfg(test)]
mod tests {
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
    fn test_ppm_window() {
        let (min, max) = ppm_window(941.407_703, 1.0);
        assert_float_absolute_eq!(min, 941.406_762);
        assert_float_absolute_eq!(max, 941.408_644);

        let (min, max) = ppm_window(941.407_703, 7.5);
        assert_float_absolute_eq!(min, 941.400_642);
        assert_float_absolute_eq!(max, 941.414_764);

        let (min, max) = ppm_window(471.711_128, 10.0);
        assert_float_absolute_eq!(min, 471.706_411);
        assert_float_absolute_eq!(max, 471.715_845);
    }
}
