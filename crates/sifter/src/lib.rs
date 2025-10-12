use std::{
    borrow::Cow,
    cmp::Ordering,
    collections::{BTreeMap, BTreeSet},
    io::{Cursor, Read},
    ops::RangeInclusive,
};

use derive_more::{From, Into};
use flate2::read::GzDecoder;
use mzdata::{
    MzMLReader,
    io::DetailLevel,
    mzpeaks::{CentroidPeak, MZPeakSetType, Tolerance},
    prelude::{IonProperties, SpectrumLike},
    spectrum::MultiLayerSpectrum,
};
use polars::{error::PolarsResult, frame::DataFrame, io::SerReader, prelude::CsvReader};
use rayon::prelude::*;

fn load_mass_database(csv: &str) -> PolarsResult<DataFrame> {
    let csv_reader = Cursor::new(csv);
    CsvReader::new(csv_reader).finish()
}

// PERF: Try out `f32` and see if that saves enough space to speed things up!
#[derive(Copy, Clone, Debug, From, Into)]
struct TotalFloat(f64);

impl Ord for TotalFloat {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

impl PartialOrd for TotalFloat {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for TotalFloat {
    fn eq(&self, other: &Self) -> bool {
        // NOTE: We must use the `.cmp()` from `Ord` here — the results of `.total_cmp()` do *not* match the default
        // implementations of `PartialEq` for floats. The `.total_cmp()` method treats `-0` and `0` as different, but
        // the default `PartialEq` for floats does not!
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for TotalFloat {}

type Mz = TotalFloat;

impl Mz {
    fn ppm_window(mz: f64, ppm: f64) -> RangeInclusive<Self> {
        let (min_mz, max_mz) = Tolerance::PPM(ppm).bounds(mz);
        Self(min_mz)..=Self(max_mz)
    }
}

// PERF: After getting some benchmarks in place, try cutting the size of this struct in half (f64 → f32) — the 64 bits
// definitely isn't needed precision-wise, so if it helps performance, it's probably worth the small amount of type
// casting! Don't forget about padding! All fields must shrink to the same size.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
struct ScanKey {
    precursor: Mz,
    scan_number: usize,
}

impl ScanKey {
    fn new(precursor: f64, scan_number: usize) -> Self {
        Self {
            precursor: Mz::from(precursor),
            scan_number,
        }
    }

    fn ppm_window(mz: f64, ppm: f64) -> RangeInclusive<Self> {
        let (min_mz, max_mz) = Tolerance::PPM(ppm).bounds(mz);
        Self::new(min_mz, 0)..=Self::new(max_mz, usize::MAX)
    }
}

#[derive(Clone, Eq, PartialEq, Debug)]
struct Peaks(BTreeSet<Mz>);

impl From<&MZPeakSetType<CentroidPeak>> for Peaks {
    fn from(peak_set: &MZPeakSetType<CentroidPeak>) -> Self {
        Self(peak_set.iter().map(|peak| peak.mz.into()).collect())
    }
}

// FIXME: Use a better error type from `thiserror`
type Result<T> = std::result::Result<T, &'static str>;

impl Peaks {
    fn filter_peaks(&self, mz: f64, ppm: f64) -> impl Iterator<Item = f64> {
        self.0.range(Mz::ppm_window(mz, ppm)).map(|mz| mz.0)
    }
}

type Minutes = TotalFloat;

// PERF: Shrinking this `start_time` to a `f32` won't shrink the size of this struct. If I want to avoid any packing
// anywhere, then I should move `start_time` back to `ScanInfo`, then shrink all of those fields to 4 bytes. But I'll
// need to benchmark if *increasing* the key size, whilst *decreasing* the overall K + V size actually helps
// performance. It's possible that smaller keys are faster anyways!
#[derive(Clone, Eq, PartialEq, Debug)]
struct ScanValue {
    start_time: Minutes,
    peaks: Peaks,
}

impl ScanValue {
    fn new(start_time: f64, peaks: Peaks) -> Self {
        Self {
            start_time: Minutes::from(start_time),
            peaks,
        }
    }
}

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
    use std::{
        collections::BTreeSet,
        ops::{Bound, RangeBounds},
    };

    use assert_float_eq::assert_float_absolute_eq;
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
        let Bound::Included(&TotalFloat(start)) = window.start_bound() else {
            panic!("expected a `.start_bound()` matching `Bound::Included(&Mz(_))`");
        };
        let Bound::Included(&TotalFloat(end)) = window.end_bound() else {
            panic!("expected an `.end_bound()` matching `Bound::Included(&Mz(_))`");
        };
        assert_float_absolute_eq!(start, 941.400_642);
        assert_float_absolute_eq!(end, 941.414_764);
        assert!(window.contains(&Mz::from(941.407_703)));
    }

    #[test]
    fn scan_key_ppm_window() {
        let window = ScanKey::ppm_window(471.711_128, 10.0);
        let Bound::Included(&ScanKey {
            precursor: TotalFloat(start_mz),
            scan_number: start_scan,
        }) = window.start_bound()
        else {
            panic!("expected a `.start_bound()` matching `Bound::Included(&Mz(_))`");
        };
        let Bound::Included(&ScanKey {
            precursor: TotalFloat(end_mz),
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
            .map(Mz::from),
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
