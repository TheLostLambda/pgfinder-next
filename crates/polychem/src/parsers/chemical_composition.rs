// External Crate Imports
use nom::{
    branch::alt,
    character::complete::char,
    combinator::{cut, map, opt, recognize},
    multi::many1,
    sequence::{delimited, pair},
    Parser,
};
use nom_miette::{expect, into, map_res, wrap_err};

// Local Crate Imports
use super::{
    errors::{ParseResult, PolychemErrorKind, UserErrorKind},
    primitives::{count, lowercase, offset_kind, uppercase},
};
use crate::{
    AtomicDatabase, ChemicalComposition, Count, Element, MassNumber, OffsetKind, Particle,
};

// Public API ==========================================================================================================

/// Chemical Composition
///   = { Atomic Offset }- , [ Offset Kind , Particle Offset ]
///   | Particle Offset
///   ;
pub fn chemical_composition<'a, 's, K: UserErrorKind>(
    db: &'a AtomicDatabase,
) -> impl FnMut(&'s str) -> ParseResult<'s, ChemicalComposition<'a>, K> {
    let chemical_formula = many1(atomic_offset(db));
    let optional_particle_offset = opt(pair(offset_kind, cut(particle_offset(db))));
    let atoms_and_particles = map(
        pair(chemical_formula, optional_particle_offset),
        |(chemical_formula, particle_offset)| {
            let particle_offset = particle_offset.map(|(k, (c, p))| (k, c, p));
            ChemicalComposition {
                chemical_formula,
                particle_offset,
            }
        },
    );

    let just_particles = map(particle_offset(db), |(count, particle)| {
        let particle_offset = Some((OffsetKind::Add, count, particle));
        ChemicalComposition {
            chemical_formula: Vec::new(),
            particle_offset,
        }
    });

    let parser = alt((atoms_and_particles, just_particles));
    into(wrap_err(
        parser,
        PolychemErrorKind::ExpectedChemicalComposition,
    ))
}

// Private Sub-Parsers =================================================================================================

/// Atomic Offset = ( Element | Isotope ) , [ Count ] ;
fn atomic_offset<'a, 's>(
    db: &'a AtomicDatabase,
) -> impl FnMut(&'s str) -> ParseResult<'s, (Element<'a>, Count)> {
    let element_or_isotope = alt((element(db), isotope(db)));
    let optional_count = opt(count).map(Option::unwrap_or_default);
    let parser = pair(element_or_isotope, optional_count);
    wrap_err(parser, PolychemErrorKind::ExpectedAtomicOffset)
}

/// Particle Offset = [ Count ] , Particle ;
fn particle_offset<'a, 's>(
    db: &'a AtomicDatabase,
) -> impl FnMut(&'s str) -> ParseResult<'s, (Count, Particle<'a>)> {
    let optional_count = opt(count).map(Option::unwrap_or_default);
    let parser = pair(optional_count, particle(db));
    wrap_err(parser, PolychemErrorKind::ExpectedParticleOffset)
}

// ---------------------------------------------------------------------------------------------------------------------

/// Element = uppercase , [ lowercase ] ;
fn element<'a, 's>(db: &'a AtomicDatabase) -> impl FnMut(&'s str) -> ParseResult<'s, Element<'a>> {
    map_res(element_symbol, |symbol| Element::new(db, symbol))
}

// NOTE: These are not meant to be links, it's just EBNF
#[allow(clippy::doc_link_with_quotes)]
/// Isotope = "[" , Count , Element , "]" ;
fn isotope<'a, 's>(db: &'a AtomicDatabase) -> impl FnMut(&'s str) -> ParseResult<'s, Element<'a>> {
    map_res(isotope_expr, |(mass_number, symbol)| {
        Element::new_isotope(db, symbol, mass_number)
    })
}

/// Particle = lowercase ;
fn particle<'a, 's>(
    db: &'a AtomicDatabase,
) -> impl FnMut(&'s str) -> ParseResult<'s, Particle<'a>> {
    map_res(particle_symbol, |symbol| Particle::new(db, symbol))
}

// ---------------------------------------------------------------------------------------------------------------------

/// Element = uppercase , [ lowercase ] ;
fn element_symbol(i: &str) -> ParseResult<&str> {
    let parser = recognize(pair(uppercase, opt(lowercase)));
    wrap_err(parser, PolychemErrorKind::ExpectedElementSymbol)(i)
}

// NOTE: These are not meant to be links, it's just EBNF
#[allow(clippy::doc_link_with_quotes)]
/// Isotope = "[" , Count , Element , "]" ;
fn isotope_expr(i: &str) -> ParseResult<(MassNumber, &str)> {
    let opening_bracket = expect(char('['), PolychemErrorKind::ExpectedIsotopeStart);
    let mass_number = map(
        wrap_err(count, PolychemErrorKind::ExpectedMassNumber),
        |c| MassNumber(c.0),
    );
    let closing_bracket = expect(cut(char(']')), PolychemErrorKind::ExpectedIsotopeEnd);
    delimited(
        opening_bracket,
        cut(pair(mass_number, element_symbol)),
        closing_bracket,
    )(i)
}

/// Particle = lowercase ;
fn particle_symbol(i: &str) -> ParseResult<&str> {
    let parser = recognize(lowercase);
    wrap_err(parser, PolychemErrorKind::ExpectedParticleSymbol)(i)
}

// Module Tests ========================================================================================================

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use once_cell::sync::Lazy;

    use super::*;

    static DB: Lazy<AtomicDatabase> = Lazy::new(AtomicDatabase::default);

    #[test]
    fn test_element_symbol() {
        // Valid Element Symbols
        assert_eq!(element_symbol("H"), Ok(("", "H")));
        assert_eq!(element_symbol("He"), Ok(("", "He")));
        // Invalid Element Symbols
        assert!(element_symbol("p").is_err());
        assert!(element_symbol("ep").is_err());
        assert!(element_symbol("1H").is_err());
        assert!(element_symbol("+H").is_err());
        assert!(element_symbol("[H]").is_err());
        // Multiple Element Symbols
        assert_eq!(element_symbol("OH"), Ok(("H", "O")));
        assert_eq!(element_symbol("HeH"), Ok(("H", "He")));
    }

    #[test]
    fn test_element() {
        let mut element = element(&DB);
        macro_rules! assert_element_name {
            ($input:literal, $output:literal, $name:literal) => {
                assert_eq!(
                    element($input).map(|(r, e)| (r, e.name)),
                    Ok(($output, $name))
                );
            };
        }
        // Valid Elements
        assert_element_name!("H", "", "Hydrogen");
        assert_element_name!("He", "", "Helium");
        // Invalid Elements
        assert!(element("p").is_err());
        assert!(element("ep").is_err());
        assert!(element("1H").is_err());
        assert!(element("+H").is_err());
        assert!(element("[H]").is_err());
        // Non-Existent Elements
        assert!(element("X").is_err());
        assert!(element("To").is_err());
        // Multiple Elements
        assert_element_name!("OH", "H", "Oxygen");
        assert_element_name!("HeH", "H", "Helium");
    }

    #[test]
    fn test_isotope_expr() {
        macro_rules! assert_isotope_expr {
            ($input:literal, $output:literal, $mass_number:literal, $symbol:literal) => {
                assert_eq!(
                    isotope_expr($input),
                    Ok(($output, (MassNumber::new($mass_number).unwrap(), $symbol)))
                );
            };
        }
        // Valid Isotope Expressions
        assert_isotope_expr!("[1H]", "", 1, "H");
        assert_isotope_expr!("[18O]", "", 18, "O");
        assert_isotope_expr!("[37Cl]", "", 37, "Cl");
        // Invalid Isotope Expressions
        assert!(isotope_expr("H").is_err());
        assert!(isotope_expr("[H]").is_err());
        assert!(isotope_expr("[H2]").is_err());
        assert!(isotope_expr("[18OH]").is_err());
        assert!(isotope_expr("[18]").is_err());
        assert!(isotope_expr("[[18O]]").is_err());
        assert!(isotope_expr("[-18O]").is_err());
        assert!(isotope_expr("[+18O]").is_err());
        // Multiple Isotope Expressions
        assert_isotope_expr!("[13C]O2", "O2", 13, "C");
        assert_isotope_expr!("[3He]H", "H", 3, "He");
    }

    #[test]
    fn test_isotope() {
        let mut isotope = isotope(&DB);
        macro_rules! assert_isotope_name {
            ($input:literal, $output:literal, $name:literal) => {
                assert_eq!(
                    isotope($input)
                        .map(|(r, e)| (r, format!("{}-{}", e.name, e.mass_number.unwrap()))),
                    Ok(($output, $name.to_owned()))
                );
            };
        }
        // Valid Isotopes
        assert_isotope_name!("[1H]", "", "Hydrogen-1");
        assert_isotope_name!("[18O]", "", "Oxygen-18");
        assert_isotope_name!("[37Cl]", "", "Chlorine-37");
        // Invalid Isotopes
        assert!(isotope("H").is_err());
        assert!(isotope("[H]").is_err());
        assert!(isotope("[H2]").is_err());
        assert!(isotope("[18OH]").is_err());
        assert!(isotope("[18]").is_err());
        assert!(isotope("[[18O]]").is_err());
        assert!(isotope("[-18O]").is_err());
        assert!(isotope("[+18O]").is_err());
        // Non-Existent Isotopes
        assert!(isotope("[42X]").is_err());
        assert!(isotope("[99To]").is_err());
        assert!(isotope("[15C]").is_err());
        assert!(isotope("[100Tc]").is_err());
        // Multiple Isotopes
        assert_isotope_name!("[13C]O2", "O2", "Carbon-13");
        assert_isotope_name!("[3He]H", "H", "Helium-3");
    }

    #[test]
    fn test_particle() {
        let mut particle = particle(&DB);
        macro_rules! assert_particle_name {
            ($input:literal, $output:literal, $name:literal) => {
                assert_eq!(
                    particle($input).map(|(r, e)| (r, e.name)),
                    Ok(($output, $name))
                );
            };
        }
        // Valid Particles
        assert_particle_name!("p", "", "Proton");
        assert_particle_name!("e", "", "Electron");
        // Invalid Particles
        assert!(particle("P").is_err());
        assert!(particle("Ep").is_err());
        assert!(particle("1p").is_err());
        assert!(particle("+e").is_err());
        assert!(particle("[p]").is_err());
        // Non-Existent Particles
        assert!(particle("m").is_err());
        assert!(particle("g").is_err());
        // Multiple Particles
        assert_particle_name!("ep", "p", "Electron");
        assert_particle_name!("pe", "e", "Proton");
    }

    #[test]
    fn test_particle_offset() {
        let mut particle_offset = particle_offset(&DB);
        macro_rules! assert_particle_offset {
            ($input:literal, $output:literal, $count:literal, $name:literal ) => {
                assert_eq!(
                    particle_offset($input).map(|(r, (c, p))| (r, (c, p.name))),
                    Ok(($output, (Count::new($count).unwrap(), $name)))
                );
            };
        }
        // Valid Particle Offsets
        assert_particle_offset!("p", "", 1, "Proton");
        assert_particle_offset!("e", "", 1, "Electron");
        assert_particle_offset!("5p", "", 5, "Proton");
        assert_particle_offset!("5e", "", 5, "Electron");
        assert_particle_offset!("100p", "", 100, "Proton");
        assert_particle_offset!("100e", "", 100, "Electron");
        // Invalid Particle Offsets
        assert!(particle_offset("-p").is_err());
        assert!(particle_offset("-e").is_err());
        assert!(particle_offset("+p").is_err());
        assert!(particle_offset("+e").is_err());
        assert!(particle_offset("P").is_err());
        assert!(particle_offset("Ep").is_err());
        assert!(particle_offset("1+p").is_err());
        assert!(particle_offset("+0p").is_err());
        assert!(particle_offset("-02e").is_err());
        assert!(particle_offset("+-e").is_err());
        assert!(particle_offset("+42").is_err());
        assert!(particle_offset("-[p]").is_err());
        // Non-Existent Particle Offsets
        assert!(particle_offset("+m").is_err());
        assert!(particle_offset("-3g").is_err());
        // Multiple Particle Offsets
        assert_particle_offset!("3ep", "p", 3, "Electron");
        assert_particle_offset!("pe", "e", 1, "Proton");
    }

    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_atomic_offset() {
        let mut atomic_offset = atomic_offset(&DB);
        macro_rules! assert_atomic_offset {
            ($input:literal, $output:literal, $name:literal, $count:literal) => {
                assert_eq!(
                    atomic_offset($input).map(|(r, (e, c))| {
                        let name = if let Some(a) = e.mass_number {
                            format!("{}-{a}", e.name)
                        } else {
                            e.name.to_owned()
                        };
                        (r, (name, c))
                    }),
                    Ok(($output, ($name.to_owned(), Count::new($count).unwrap())))
                );
            };
        }
        // Valid Atomic Offsets
        assert_atomic_offset!("H", "", "Hydrogen", 1);
        assert_atomic_offset!("He", "", "Helium", 1);
        assert_atomic_offset!("H2", "", "Hydrogen", 2);
        assert_atomic_offset!("C18", "", "Carbon", 18);
        assert_atomic_offset!("[1H]", "", "Hydrogen-1", 1);
        assert_atomic_offset!("[18O]", "", "Oxygen-18", 1);
        assert_atomic_offset!("[37Cl]", "", "Chlorine-37", 1);
        assert_atomic_offset!("[1H]2", "", "Hydrogen-1", 2);
        assert_atomic_offset!("[18O]3", "", "Oxygen-18", 3);
        assert_atomic_offset!("[37Cl]5", "", "Chlorine-37", 5);
        // Invalid Atomic Offsets
        assert!(atomic_offset("p").is_err());
        assert!(atomic_offset("-2p").is_err());
        assert!(atomic_offset("ep").is_err());
        assert!(atomic_offset("1H").is_err());
        assert!(atomic_offset("+H").is_err());
        assert!(atomic_offset("2[2H]").is_err());
        assert!(atomic_offset("[H]").is_err());
        assert!(atomic_offset("[H2]").is_err());
        assert!(atomic_offset("[18OH]").is_err());
        assert!(atomic_offset("[18]").is_err());
        assert!(atomic_offset("[[18O]]").is_err());
        assert!(atomic_offset("[-18O]").is_err());
        assert!(atomic_offset("[+18O]").is_err());
        // Non-Existent Atomic Offsets
        assert!(atomic_offset("X2").is_err());
        assert!(atomic_offset("To7").is_err());
        assert!(atomic_offset("[42X]").is_err());
        assert!(atomic_offset("[99To]2").is_err());
        assert!(atomic_offset("[15C]4").is_err());
        assert!(atomic_offset("[100Tc]8").is_err());
        // Multiple Atomic Offsets
        assert_atomic_offset!("OH", "H", "Oxygen", 1);
        assert_atomic_offset!("HeH", "H", "Helium", 1);
        assert_atomic_offset!("H2O", "O", "Hydrogen", 2);
        assert_atomic_offset!("CO2", "O2", "Carbon", 1);
        assert_atomic_offset!("[13C]6O2", "O2", "Carbon-13", 6);
        assert_atomic_offset!("[3He]1H", "H", "Helium-3", 1);
    }

    #[test]
    #[allow(clippy::cognitive_complexity)]
    fn test_chemical_composition() {
        let mut chemical_composition = chemical_composition::<PolychemErrorKind>(&DB);
        macro_rules! check_composition_snapshot {
            ($input:literal, $output:literal) => {
                let (
                    rest,
                    ChemicalComposition {
                        chemical_formula,
                        particle_offset,
                    },
                ) = chemical_composition($input).unwrap();
                assert_eq!(rest, $output);
                let chemical_formula: Vec<_> = chemical_formula
                    .iter()
                    .map(|(e, c)| {
                        let name = if let Some(a) = e.mass_number {
                            format!("{}-{a}", e.name)
                        } else {
                            e.name.to_owned()
                        };
                        let count = c.0.get();
                        (name, count)
                    })
                    .collect();
                let particle_offset = particle_offset.map(|(k, c, p)| (k, c.0.get(), p.name));
                assert_debug_snapshot!((chemical_formula, particle_offset));
            };
        }
        // Valid Chemical Compositions
        check_composition_snapshot!("H2O", "");
        check_composition_snapshot!("C11H12N2O2", "");
        check_composition_snapshot!("OH+e", "");
        check_composition_snapshot!("Na-e", "");
        check_composition_snapshot!("NH3+p", "");
        check_composition_snapshot!("Cr2O7+2e", "");
        check_composition_snapshot!("NH2+2p", "");
        check_composition_snapshot!("D2O", "");
        check_composition_snapshot!("[2H]2O", "");
        check_composition_snapshot!("[37Cl]5-2p", "");
        check_composition_snapshot!("NH2[99Tc]", "");
        // Invalid Chemical Compositions
        assert!(chemical_composition(" ").is_err());
        assert!(chemical_composition("-2p").is_err());
        assert!(chemical_composition("+H").is_err());
        assert!(chemical_composition("2[2H]").is_err());
        assert!(chemical_composition("[H+p]O").is_err());
        assert!(chemical_composition("NH2[100Tc]").is_err());
        // Multiple Chemical Compositions
        check_composition_snapshot!("[37Cl]5-2p10", "10");
        check_composition_snapshot!("[2H]2O*H2O", "*H2O");
        check_composition_snapshot!("NH2{100Tc", "{100Tc");
        check_composition_snapshot!("C11H12N2O2 H2O", " H2O");
    }
}
