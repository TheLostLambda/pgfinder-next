use std::num::TryFromIntError;

use miette::Diagnostic;
use thiserror::Error;

// FIXME: All of the errors in `polychem` need some major refactoring: which files / submodules should have their own
// error enums? How should they be named?
// FIXME: This should probably be hidden behind and opaque unit-struct, like `polymerizer::Error` does for
// `PolymerizerError` — then this doesn't need to be public API...
#[derive(Debug, Diagnostic, Clone, Eq, PartialEq, Error)]
pub enum OffsetMultiplierError {
    #[error(
        "attemped to convert a `SignedCount` into an `OffsetMultiplier`, but the count was zero"
    )]
    Zero,

    #[error("attemped to convert a `SignedCount` into an `OffsetMultiplier`, but the count was too large to represent")]
    TooLarge(#[from] TryFromIntError),
}
