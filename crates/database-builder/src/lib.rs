// NOTE: For Tia — you'll need to import Polars like this for it to compile to WebAssembly
// https://www.reddit.com/r/rust/comments/yb58hk/how_to_use_polars_with_wasm/
// That's because, by default, `polars` will import some code that interacts with the operating
// system / terminal, and those APIs aren't available on the web!
// use polars_core::prelude::*;
// use polars_lazy::prelude::*;

#[derive(Clone, Debug, Default)]
struct DatabaseBuilder {
    glycan_residues: Vec<char>,
    stem_residues: Vec<char>,
}

impl DatabaseBuilder {
    fn new() -> Self {
        Self::default()
    }

    fn add_glycan(&mut self, residue: char) {
        self.glycan_residues.push(residue);
    }

    fn add_stem_residues(&mut self, residues: &[char]) {
        self.stem_residues.extend(residues);
    }

    // NOTE: Probably output DataFrame or CSV text?
    fn generate(&self) -> Vec<String> {
        let mut result = Vec::new();
        for glycan in &self.glycan_residues {
            for stem in &self.stem_residues {
                result.push(format!("{glycan}-{stem}"));
            }
        }
        result
    }
}

use muropeptide::{Muropeptide, POLYMERIZER, parser::muropeptide};
use nom_miette::final_parser;

fn hello(name: &mut String) -> String {
    *name += " Duh";
    format!("Hello, {name}!")
}

pub fn parse_muropeptide(
    structure: impl AsRef<str>,
) -> muropeptide::Result<Muropeptide<'static, 'static>> {
    let muropeptide = final_parser(muropeptide(&POLYMERIZER))(structure.as_ref())?;

    Ok(muropeptide)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let mut name = "Tia".to_owned();
        assert_eq!(hello(&mut name), "Hello, Tia Duh!");
        assert_eq!(hello(&mut name), "Hello, Tia Duh Duh!");
        // assert_eq!(hello(name), "Hello, Tia Duh!");
        // assert_eq!(hello(name), "Hello, Tia!");
    }

    #[test]
    fn bark() {
        let mut database_builder = DatabaseBuilder::new();
        assert!(database_builder.glycan_residues.is_empty());
        assert!(database_builder.stem_residues.is_empty());

        database_builder.add_glycan('g');
        assert_eq!(database_builder.glycan_residues, vec!['g']);

        database_builder.add_stem_residues(&['A', 'E', 'J']);
        assert_eq!(database_builder.stem_residues, vec!['A', 'E', 'J']);

        assert_eq!(database_builder.generate(), vec!["g-A", "g-E", "g-J"]);
    }
}
