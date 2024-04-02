use miette::Diagnostic;
use thiserror::Error;

use crate::Count;

// FIXME: All of the errors in `polychem` need some major refactoring: which files / submodules should have their own
// error enums? How should they be named?
// FIXME: This should probably be hidden behind and opaque unit-struct, like `polymerizer::Error` does for
// `PolymerizerError` — then this doesn't need to be public API...
#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
pub enum OffsetMultiplierError {
    #[error("counts must be non-zero")]
    Zero,

    #[error("the count {0} is too large to represent using {} bits", Count::BITS)]
    TooLarge(u64),
}
