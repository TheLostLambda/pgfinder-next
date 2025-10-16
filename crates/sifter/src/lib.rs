mod ms2_index;
mod peaks;
mod scan_kv;
mod total_float;

use std::io::Cursor;

use polars::{error::PolarsResult, frame::DataFrame, io::SerReader, prelude::CsvReader};

fn load_mass_database(csv: &str) -> PolarsResult<DataFrame> {
    let csv_reader = Cursor::new(csv);
    CsvReader::new(csv_reader).finish()
}

// FIXME: Use a better error type from `thiserror`
type Result<T> = std::result::Result<T, &'static str>;

#[cfg(test)]
mod tests {
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
}
