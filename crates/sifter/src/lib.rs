mod found_fragment;
mod ms2_index;
mod named_ion;
mod ordered_floats;
mod peaks;
mod ppm_window;
mod scan_kv;

use std::{borrow::Cow, collections::BTreeMap};

use derive_more::Constructor;

use crate::scan_kv::{ScanKey, ScanValue};

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Ms2Index(BTreeMap<ScanKey, ScanValue>);

#[derive(Clone, PartialEq, PartialOrd, Debug, Constructor)]
pub struct FoundFragment<'n> {
    theoretical_precursor: NamedIon<'n>,
    theoretical_fragment: NamedIon<'n>,
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
