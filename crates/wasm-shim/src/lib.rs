use consolidation::consolidate_replicates;
use muropeptide::{Muropeptide, POLYMERIZER, SMILES_DB};
use polychem::{ChargedParticle, Massive};
use rust_decimal::Decimal;
use smithereens::Dissociable;
use std::fmt::Write;
use wasm_bindgen::prelude::*;

const VERSION: &str = env!("CARGO_PKG_VERSION");

#[wasm_bindgen]
#[must_use]
pub fn version() -> String {
    VERSION.to_owned()
}

#[wasm_bindgen]
struct Peptidoglycan(Muropeptide<'static, 'static>);

#[wasm_bindgen]
impl Peptidoglycan {
    #[wasm_bindgen(constructor)]
    #[allow(clippy::use_self)]
    pub fn new(structure: &str) -> Result<Peptidoglycan, String> {
        // NOTE: This ensures the panic hook is set before any other shim code can be run!
        console_error_panic_hook::set_once();
        Muropeptide::new(&POLYMERIZER, structure)
            .map(Self)
            .map_err(|e| e.to_string())
    }

    #[must_use]
    pub fn oligomerization_state(&self) -> usize {
        self.0.oligomerization_state()
    }

    #[must_use]
    pub fn monoisotopic_mass(&self) -> String {
        let mass = self.0.monoisotopic_mass();
        decimal_round_workaround(mass, 6)
    }

    #[must_use]
    pub fn smiles(&self) -> String {
        // FIXME: This API needs some work — this should probably be a method on `Muropeptide`?
        SMILES_DB.muropeptide_to_smiles(&self.0).unwrap_or_default()
    }

    #[must_use]
    pub fn fragment(&self) -> String {
        let mut fragments: Vec<_> = self
            .0
            .fragment(None)
            .map(|fragment| (fragment.to_string(), fragment.monoisotopic_mz().unwrap()))
            .collect();
        fragments.sort_unstable_by(|(s1, mz1), (s2, mz2)| {
            mz1.cmp(mz2).then_with(|| s1.cmp(s2)).reverse()
        });
        fragments.dedup();

        let mut csv = String::new();
        writeln!(&mut csv, "Structure,Ion M/Z").unwrap();
        for (structure, mz) in fragments {
            let mz = decimal_round_workaround(mz, 6);
            writeln!(&mut csv, r#""{structure}",{mz}"#).unwrap();
        }
        csv
    }
}

#[wasm_bindgen]
pub struct Replicate {
    number: u32,
    csv: String,
}

#[wasm_bindgen]
#[must_use]
pub fn consolidate(replicates: Vec<Replicate>) -> String {
    // PERF: Could do this without the `collect()` — just pass on the iterator
    let replicates: Vec<_> = replicates.into_iter().map(|r| (r.number, r.csv)).collect();
    consolidate_replicates(&replicates)
}

// FIXME: Really this should be fixed in `rust_decimal`...
fn decimal_round_workaround(value: impl Into<Decimal>, decimal_points: u32) -> String {
    let value = value.into().round_dp(decimal_points);
    format!("{value}")
}
