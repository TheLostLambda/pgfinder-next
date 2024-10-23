// NOTE: For Tia — you'll need to import Polars like this for it to compile to WebAssembly
// https://www.reddit.com/r/rust/comments/yb58hk/how_to_use_polars_with_wasm/
// That's because, by default, `polars` will import some code that interacts with the operating
// system / terminal, and those APIs aren't available on the web!
// use polars_core::prelude::*;
// use polars_lazy::prelude::*;

#[must_use]
pub const fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
