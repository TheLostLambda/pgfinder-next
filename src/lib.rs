use phf::phf_map;
use rust_decimal::prelude::*;
use rust_decimal_macros::dec;
use nom::{
    bytes::complete::tag,
    character::complete::{char, one_of, space0},
    combinator::map,
    multi::separated_list1,
    sequence::{delimited, tuple},
    IResult,
};

static MOIETIES: phf::Map<&str, Moiety> = phf_map! {
  "Ac" => Moiety::new("Ac"),
};

pub fn monosaccharide(i: &str) -> IResult<&str, char> {
    one_of("gm")(i)
}

pub fn amino_acid(i: &str) -> IResult<&str, char> {
    one_of("AEJ")(i)
}

pub fn residue_modification(i: &str) -> IResult<&str, Vec<Modification>> {
    map(
        delimited(
            char('('),
            separated_list1(char(','), tuple((space0, one_of("+-"), modification, space0))),
            char(')'),
        ),
        |mods| {
            mods.iter()
                .map(|(_, sign, moiety, _)| match sign {
                    '+' => Modification::Add(*moiety),
                    '-' => Modification::Remove(*moiety),
                    _ => unreachable!(),
                })
                .collect()
        },
    )(i)
}

pub fn modification(i: &str) -> IResult<&str, Moiety> {
    map(tag("Ac"), |moiety| MOIETIES[moiety])(i)
}

pub struct Residue {
    id: usize,
    moiety: Moiety,
    mods: Vec<Modification>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Modification {
    Add(Moiety),
    Remove(Moiety),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Moiety {
    abbr: &'static str,
    name: &'static str,
    mass: Decimal,
}

impl Moiety {
    // FIXME: Replace this with a real implementation
    pub const fn new(abbr: &'static str) -> Self {
        Self {
            abbr,
            name: "",
            mass: dec!(0.0),
        }
    }
}

pub fn add(left: usize, right: usize) -> usize {
    left + right
}

#[cfg(test)]
mod tests {
    use std::error::Error;

    use super::*;

    #[test]
    fn it_works() -> Result<(), Box<dyn Error>> {
        assert_eq!(residue_modification("(-Ac, +Ac)")?, Default::default());
        Ok(())
    }
}
