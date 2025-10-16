mod peaks;
mod scan_kv;
mod total_float;

use std::{
    borrow::Cow,
    collections::BTreeMap,
    io::{Cursor, Read},
};

use flate2::read::GzDecoder;
use mzdata::{
    MzMLReader,
    io::DetailLevel,
    prelude::{IonProperties, SpectrumLike},
    spectrum::MultiLayerSpectrum,
};
use polars::{error::PolarsResult, frame::DataFrame, io::SerReader, prelude::CsvReader};
use rayon::prelude::*;

use crate::scan_kv::{ScanKey, ScanValue};

fn load_mass_database(csv: &str) -> PolarsResult<DataFrame> {
    let csv_reader = Cursor::new(csv);
    CsvReader::new(csv_reader).finish()
}

// FIXME: Use a better error type from `thiserror`
type Result<T> = std::result::Result<T, &'static str>;

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct Ms2Index(BTreeMap<ScanKey, ScanValue>);

impl Ms2Index {
    pub fn from_bytes(bytes: impl AsRef<[u8]>) -> Result<Self> {
        let bytes = Self::maybe_decode_bytes(bytes.as_ref())?;

        let mzml = MzMLReader::with_buffer_capacity_and_detail_level(
            bytes.as_ref(),
            10_000,
            DetailLevel::Lazy,
        );
        let spectra: Vec<_> = mzml.collect();

        Self::from_spectra(spectra)
    }

    pub fn from_spectra(
        spectra: impl IntoParallelIterator<Item = MultiLayerSpectrum>,
    ) -> Result<Self> {
        spectra
            .into_par_iter()
            .filter(|spectrum| spectrum.ms_level() == 2)
            .map(|mut spectrum| {
                let precursor = spectrum
                    .precursor()
                    .ok_or("MS2 spectrum was missing precursor ion information")?
                    .mz();
                // NOTE: This is assuming that the mzML file we've been given represents a full MS run — if this file
                // is a slice of a larger run, then the scan `.id()` will differ from this!
                let scan_number = spectrum.index() + 1;
                let start_time = spectrum.start_time();

                let peaks = spectrum
                    .try_build_centroids()
                    .map_err(|_| "failed to find centroided peak data")?
                    .into();

                Ok((
                    ScanKey::new(precursor, scan_number),
                    ScanValue::new(start_time, peaks),
                ))
            })
            .collect::<Result<_>>()
            .map(Ms2Index)
    }

    fn maybe_decode_bytes(bytes: &[u8]) -> Result<Cow<'_, [u8]>> {
        let mut gz_decoder = GzDecoder::new(bytes);
        if gz_decoder.header().is_some() {
            // NOTE: Our decompressed data should be *at least* as big as the compressed data, so we can use the size
            // of the compressed data to allocate some initial capacity
            let mut decoded_bytes = Vec::with_capacity(bytes.len());
            gz_decoder
                .read_to_end(&mut decoded_bytes)
                .map_err(|_| "failed to decompress bytes")?;
            Ok(Cow::Owned(decoded_bytes))
        } else {
            Ok(Cow::Borrowed(bytes))
        }
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use mzdata::{MzMLReader, io::DetailLevel};
    use polars::prelude::DataType;

    use super::*;

    const MASS_DATABASE_CSV: &str = include_str!("../tests/data/Exchange.csv");
    const MZML: &[u8] = include_bytes!("../tests/data/WT (6.7–7.3 min).mzML");
    const MZML_GZ: &[u8] = include_bytes!("../tests/data/WT (6.7–7.3 min).mzML.gz");

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
    fn ms2_index_from() {
        let mzml =
            MzMLReader::with_buffer_capacity_and_detail_level(MZML, 10_000, DetailLevel::Lazy);
        let spectra: Vec<_> = mzml.collect();
        let ms2_index = Ms2Index::from_spectra(spectra).unwrap();
        assert_debug_snapshot!(ms2_index);

        let from_bytes = Ms2Index::from_bytes(MZML).unwrap();
        assert_eq!(from_bytes, ms2_index);
        let from_gzipped_bytes = Ms2Index::from_bytes(MZML_GZ).unwrap();
        assert_eq!(from_gzipped_bytes, ms2_index);
    }
}
