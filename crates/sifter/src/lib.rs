mod ms2_index;
mod ordered_floats;
mod peaks;
mod ppm_window;
mod scan_kv;

use std::collections::BTreeMap;

use crate::scan_kv::{ScanKey, ScanValue};

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Ms2Index(BTreeMap<ScanKey, ScanValue>);

// FIXME: Use a better error type from `thiserror`
type Result<T> = std::result::Result<T, &'static str>;
