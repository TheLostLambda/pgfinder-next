// FIXME: This should probably be made private at some point!
pub mod parsers;

// static MOIETIES: phf::Map<&str, Moiety> = phf_map! {
//   "Ac" => Moiety::new("Ac"),
// };
//
// pub fn residue_modification(i: &str) -> IResult<&str, Vec<Modification>> {
//     map(
//         delimited(
//             char('('),
//             separated_list1(char(','), tuple((space0, one_of("+-"), modification, space0))),
//             char(')'),
//         ),
//         |mods| {
//             mods.iter()
//                 .map(|(_, sign, moiety, _)| match sign {
//                     '+' => Modification::Add(*moiety),
//                     '-' => Modification::Remove(*moiety),
//                     _ => unreachable!(),
//                 })
//                 .collect()
//         },
//     )(i)
// }
//
// pub fn modification(i: &str) -> IResult<&str, Moiety> {
//     map(tag("Ac"), |moiety| MOIETIES[moiety])(i)
// }
//
// pub struct Residue {
//     id: usize,
//     moiety: Moiety,
//     mods: Vec<Modification>,
// }
//
// #[derive(Clone, Debug, PartialEq)]
// pub enum Modification {
//     Add(Moiety),
//     Remove(Moiety),
// }
//
// #[derive(Clone, Copy, Debug, PartialEq)]
// pub struct Moiety {
//     abbr: &'static str,
//     name: &'static str,
//     mass: Decimal,
// }
//
// impl Moiety {
//     // FIXME: Replace this with a real implementation
//     pub const fn new(abbr: &'static str) -> Self {
//         Self {
//             abbr,
//             name: "",
//             mass: dec!(0.0),
//         }
//     }
// }
//
// #[cfg(test)]
// mod tests {
//     use std::error::Error;
//
//     use super::*;
//
//     #[test]
//     fn it_works() -> Result<(), Box<dyn Error>> {
//         assert_eq!(residue_modification("(-Ac, +Ac)")?, Default::default());
//         Ok(())
//     }
// }
