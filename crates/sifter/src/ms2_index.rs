// Standard Library Imports
use std::{borrow::Cow, io::Read};

// External Crate Imports
use flate2::read::GzDecoder;
use mzdata::{
    MzMLReader,
    io::DetailLevel,
    prelude::{IonProperties, SpectrumLike},
    spectrum::MultiLayerSpectrum,
};
use rayon::prelude::*;

// Local Crate Imports
use crate::{
    FoundPrecursor, Ms2Index, NamedIon, Result,
    ppm_window::PpmWindow,
    scan_kv::{ScanKey, ScanValue},
};

// Public API ==========================================================================================================

impl Ms2Index {
    pub fn from_bytes(bytes: impl AsRef<[u8]>) -> Result<Self> {
        let bytes = Self::decode_bytes_if_compressed(bytes.as_ref())?;

        // FIXME: I think `mzdata` should have a way to construct an `MzMLReader` *without* a buffer! For now, the
        // buffer capacity is `10,000` bytes, completely arbitrarily
        let mzml = MzMLReader::with_buffer_capacity_and_detail_level(
            bytes.as_ref(),
            10_000,
            DetailLevel::Lazy,
        );

        // TODO: Benchmark if this `.collect()` is faster or slower than `.par_bridge()`
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

    // PERF: Add a `par_find_precursors()` that takes and returns a parallel iterator (if benchmarks justify it)
    pub fn find_precursors<'p, 'n: 'p>(
        &self,
        precursors: impl IntoIterator<Item = &'p NamedIon<'n>>,
        ppm_tolerance: impl Into<f64>,
    ) -> impl Iterator<Item = FoundPrecursor<'p, 'n>> {
        let ppm_tolerance = ppm_tolerance.into();
        precursors.into_iter().flat_map(move |theoretical| {
            self.filter_scans(theoretical.mz(), ppm_tolerance).map(
                |(
                    &ScanKey {
                        precursor,
                        scan_number,
                    },
                    &ScanValue { start_time, .. },
                )| {
                    FoundPrecursor::new(theoretical, precursor, scan_number, start_time)
                },
            )
        })
    }
}

// Private Methods =====================================================================================================

impl Ms2Index {
    fn decode_bytes_if_compressed(bytes: &[u8]) -> Result<Cow<'_, [u8]>> {
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

    fn filter_scans(
        &self,
        precursor_mz: f64,
        ppm_tolerance: f64,
    ) -> impl Iterator<Item = (&ScanKey, &ScanValue)> {
        self.0
            .range(ScanKey::ppm_window(precursor_mz, ppm_tolerance))
    }
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use mzdata::{MzMLReader, io::DetailLevel};

    use super::*;

    const MZML: &[u8] = include_bytes!("../tests/data/WT (6.7–7.3 min).mzML");
    const MZML_GZ: &[u8] = include_bytes!("../tests/data/WT (6.7–7.3 min).mzML.gz");

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

    #[test]
    fn ms2_index_find_precursors() {
        let ms2_index = Ms2Index::from_bytes(MZML).unwrap();
        let monomer_ions = [
            NamedIon::new("g(r)m-AEJA(+p)", 942.414_979),
            NamedIon::new("g(r)m-AEJA(+2p)", 471.711_127_5),
        ];
        let found_precursors: Vec<_> = ms2_index.find_precursors(&monomer_ions, 10).collect();
        assert_debug_snapshot!(found_precursors);
    }
}
