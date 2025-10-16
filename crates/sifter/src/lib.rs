mod found_fragment;
mod found_precursor;
mod ms2_index;
mod named_ion;
mod ordered_floats;
mod peaks;
mod ppm_window;
mod scan_kv;

// Standard Library Imports
use std::{borrow::Cow, collections::BTreeMap};

// External Crate Imports
use derive_more::Constructor;

// Local Crate Imports
use crate::scan_kv::{ScanKey, ScanValue};

// Public API ==========================================================================================================

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Ms2Index(BTreeMap<ScanKey, ScanValue>);

#[derive(Clone, PartialEq, PartialOrd, Debug, Constructor)]
pub struct FoundPrecursor<'p, 'n> {
    theoretical: &'p NamedIon<'n>,
    observed_mz: f64,
    scan_number: usize,
    start_time: f64,
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Constructor)]
pub struct FoundFragment<'p, 'f, 'n> {
    theoretical_precursor: &'p NamedIon<'n>,
    theoretical_fragment: &'f NamedIon<'n>,
    observed_precursor_mz: f64,
    observed_fragment_mz: f64,
    scan_number: usize,
    start_time: f64,
}

#[derive(Clone, PartialEq, PartialOrd, Debug)]
pub struct NamedIon<'n> {
    name: Cow<'n, str>,
    mz: f64,
}

// FIXME: Use a better error type from `thiserror`
type Result<T> = std::result::Result<T, &'static str>;
