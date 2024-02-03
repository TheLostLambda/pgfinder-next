use std::{collections::HashMap, fmt};

use miette::{Diagnostic, LabeledSpan, SourceSpan};
use nom::{
    combinator::{all_consuming, complete, consumed, map_res, success},
    error::{FromExternalError, ParseError},
    Finish, IResult, Parser,
};
use thiserror::Error;

// FIXME: Move all of this error handling and context wrapping to another crate to be shared
#[derive(Debug, Clone, Eq, PartialEq, Error)]
#[error("{error}")]
pub struct LabeledError<E: LabeledErrorKind> {
    full_input: String,
    labels: Vec<LabeledSpan>,
    error: ErrorTree<E>,
}

#[derive(Debug, Clone, Eq, PartialEq, Error)]
pub enum ErrorTree<E: LabeledErrorKind> {
    #[error("{kind}")]
    Node {
        kind: E,
        #[source]
        source: Option<Box<LabeledError<E>>>,
    },
    #[error("attempted {} parse branches unsuccessfully", .0.len())]
    Branch(Vec<LabeledError<E>>),
}

impl<E: LabeledErrorKind> Diagnostic for LabeledError<E> {
    fn source_code(&self) -> Option<&dyn miette::SourceCode> {
        Some(&self.full_input)
    }

    fn help<'a>(&'a self) -> Option<Box<dyn fmt::Display + 'a>> {
        if let ErrorTree::Node { kind, .. } = &self.error {
            kind.help()
        } else {
            None
        }
    }

    fn labels(&self) -> Option<Box<dyn Iterator<Item = miette::LabeledSpan> + '_>> {
        Some(Box::new(self.labels.iter().cloned()))
    }

    fn related<'a>(&'a self) -> Option<Box<dyn Iterator<Item = &'a dyn Diagnostic> + 'a>> {
        if let ErrorTree::Branch(related) = &self.error {
            Some(Box::new(related.iter().map(|e| e as &dyn Diagnostic)))
        } else {
            None
        }
    }

    fn diagnostic_source(&self) -> Option<&dyn Diagnostic> {
        // FIXME: Better way to do this?
        if let ErrorTree::Node { source, .. } = &self.error {
            source.as_ref().map(|e| &**e as &dyn Diagnostic)
        } else {
            None
        }
    }
}

impl<E: LabeledErrorKind> LabeledError<E> {
    // FIXME: Build the whole labelled span here?
    // FIXME: Is mutable best here?
    fn bubble_labels(&mut self) {
        if self.labels.is_empty() {
            match &mut self.error {
                ErrorTree::Node {
                    source: Some(child),
                    ..
                } => {
                    child.bubble_labels();
                    self.labels = child.labels.drain(..).collect();
                }
                ErrorTree::Branch(alternatives) => {
                    let new_labels = alternatives.iter_mut().flat_map(|child| {
                        child.bubble_labels();
                        child.labels.drain(..)
                    });
                    self.labels = Self::merge_labels(new_labels);
                }
                ErrorTree::Node { .. } => (),
            }
        }
    }

    // FIXME: Where in the world does this belong...
    // FIXME: This eventually needs testing — showing that labels with different spans *don't* get merged
    fn merge_labels(labels: impl Iterator<Item = LabeledSpan>) -> Vec<LabeledSpan> {
        let mut span_map: HashMap<SourceSpan, Vec<String>> = HashMap::new();
        for labeled_span in labels {
            let span = labeled_span.inner();
            // FIXME: Gross with the clone() and to_owned() in here...
            let label = labeled_span.label().unwrap().to_owned();
            span_map
                .entry(*span)
                .and_modify(|l| l.push(label.clone()))
                .or_insert_with(|| vec![label]);
        }
        span_map
            .into_iter()
            .map(|(span, labels)| {
                let label = labels.join(" or ");
                LabeledSpan::new_with_span(Some(label), span)
            })
            .collect()
    }
}

// FIXME: Abstract this out into it's own module (or nom-supreme?)
// FIXME: Check that field ordering everywhere matches this!
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct LabeledParseError<'a, E> {
    input: &'a str,
    length: usize,
    // FIXME: Get rid of this too — expect to wipe error stack, then wrap_err for anything above
    reported: bool,
    // FIXME: Get rid of this?
    fatal: bool,
    kind: E,
    // FIXME: Replace this with the ErrorTree!!
    alternatives: Vec<LabeledParseError<'a, E>>,
    source: Option<Box<LabeledParseError<'a, E>>>,
}

pub trait LabeledErrorKind: Diagnostic + Clone + Eq + From<nom::error::ErrorKind> {
    fn label(&self) -> Option<&'static str>;
}

// FIXME: Oh lord, what a mess...
impl<'a, E: LabeledErrorKind> LabeledParseError<'a, E> {
    // FIXME: Add a version that doesn't keep the source?
    // FIXME: Instead of all of this rebuilding with .., just take `mut self` like miette's GraphicalReportHandler
    // FIXME: I don't like this being public!!!
    pub fn kind(self, kind: E) -> Self {
        // FIXME: This erases any non-final kinds... Is that okay?
        let source = if self.reported {
            Some(Box::new(self.clone()))
        } else {
            None
        };
        Self {
            kind,
            reported: true,
            source,
            // FIXME: Boy, we really need a default constructor...
            alternatives: Vec::new(),
            ..self
        }
    }

    // FIXME: I don't like this being public!!!
    pub fn fatal(self, fatal: bool) -> Self {
        Self { fatal, ..self }
    }

    // FIXME: God help with the naming... Arguments too...
    fn into_final_error(self, full_input: &str) -> LabeledError<E> {
        let source = self
            .source
            .clone()
            .map(|e| Box::new(e.into_final_error(full_input)));
        let alternatives: Vec<_> = self
            .alternatives
            .iter()
            .cloned()
            .map(|e| e.into_final_error(full_input))
            .collect();
        // NOTE: The additional space is added so that Diagnostic labels can point to the end of an input
        let input = format!("{full_input} ");
        let span = span_from_input(full_input, self.input, self.length);
        let label = self.kind.label().map(str::to_string);
        let labels = label
            .into_iter()
            .map(|l| LabeledSpan::new_with_span(Some(l), span))
            .collect();
        // FIXME: WET CODE!!!
        if alternatives.is_empty() {
            LabeledError {
                full_input: input,
                labels,
                error: ErrorTree::Node {
                    kind: self.kind,
                    source,
                },
            }
        } else {
            // FIXME: Yuckie!
            let mut other_self = self;
            other_self.alternatives.clear();
            LabeledError {
                full_input: input,
                labels: Vec::new(),
                error: ErrorTree::Branch(
                    [vec![other_self.into_final_error(full_input)], alternatives].concat(),
                ),
            }
        }
    }
}

// FIXME: OMG refactor... Also, can I do better than stealing here?
fn span_from_input(full_input: &str, input: &str, length: usize) -> SourceSpan {
    let base_addr = full_input.as_ptr() as usize;
    let substr_addr = input.as_ptr() as usize;
    // FIXME: Keep this?
    assert!(
        substr_addr >= base_addr,
        "tried to get the span of a non-substring!"
    );
    let start = substr_addr - base_addr;
    let end = start + length;
    SourceSpan::from(start..end)
}

pub fn final_parser<'a, O, P, E>(parser: P) -> impl FnMut(&'a str) -> Result<O, LabeledError<E>>
where
    E: LabeledErrorKind,
    P: Parser<&'a str, O, LabeledParseError<'a, E>>,
{
    // FIXME: I can't inline this because of some borrow-checker closure witchcraft...
    let mut parser = all_consuming(complete(parser));
    move |input| {
        parser.parse(input).finish().map(|(_, c)| c).map_err(|e| {
            // FIXME: Ew
            let mut error = e.into_final_error(input);
            error.bubble_labels();
            error
        })
    }
}

// FIXME: Just rename to map_res? Make always fatal?
pub fn map_res_span<'a, O1, O2, E1, E2, F, G>(
    parser: F,
    f: G,
) -> impl FnMut(&'a str) -> IResult<&'a str, O2, LabeledParseError<'a, E1>>
where
    O1: Clone,
    E1: LabeledErrorKind,
    F: Copy + Parser<&'a str, O1, LabeledParseError<'a, E1>>,
    G: Copy + FnMut(O1) -> Result<O2, E2>,
    LabeledParseError<'a, E1>: FromExternalError<&'a str, E2>,
{
    move |input| {
        let i = input;
        let (input, (consumed, o1)) = consumed(parser)(input)?;
        wrap_err(map_res(success(o1), f), |e| LabeledParseError {
            input: i,
            length: consumed.len(),
            ..e
        })(input)
    }
}

// FIXME: Check if I'm being consistent about using `impl` or generics... I think I should avoid any generics I don't
// use, as long as this is just a sort of "internal" library
// FIXME: See if this signature can be simplified (elide lifetimes?)
// FIXME: Do I really need a closure? Or just a kind?
// FIXME: Standardize the order of all of these generic arguments!
pub fn wrap_err<'a, O, P, F, E>(
    mut parser: P,
    f: F,
) -> impl FnMut(&'a str) -> IResult<&'a str, O, LabeledParseError<'a, E>>
where
    E: LabeledErrorKind,
    P: Parser<&'a str, O, LabeledParseError<'a, E>>,
    F: Clone + FnOnce(LabeledParseError<'a, E>) -> LabeledParseError<'a, E>,
{
    move |i| {
        parser.parse(i).map_err(|e| match e.map(f.clone()) {
            nom::Err::Error(e) if e.fatal => nom::Err::Failure(e),
            rest => rest,
        })
    }
}

// FIXME: I should just replace wrap_err with expect, or I should make expect *not* wrap things
pub fn expect<'a, O, E: LabeledErrorKind, F>(
    parser: F,
    error: E,
) -> impl FnMut(&'a str) -> IResult<&'a str, O, LabeledParseError<'a, E>>
where
    F: Parser<&'a str, O, LabeledParseError<'a, E>>,
{
    wrap_err(parser, |e| e.kind(error))
}

// FIXME: Eventually, I should make everything generic over the input type again... So you'd be able to use this
// library with &'a [u8] like `nom` lets you
impl<'a, E: LabeledErrorKind> ParseError<&'a str> for LabeledParseError<'a, E> {
    fn from_error_kind(input: &'a str, kind: nom::error::ErrorKind) -> Self {
        // FIXME: Trying implementing Default and doing ..Self::default()
        Self {
            input,
            length: 0,
            kind: kind.into(),
            reported: false,
            fatal: false,
            alternatives: Vec::new(),
            source: None,
        }
    }

    fn append(_input: &str, _kind: nom::error::ErrorKind, other: Self) -> Self {
        other
    }

    fn or(self, other: Self) -> Self {
        // FIXME: Oh god...
        let (new_self, alternative) = match (self.reported, other.reported) {
            (true, true) => (self.clone(), vec![other.clone()]),
            (true, false) => (self.clone(), Vec::new()),
            _ => (other.clone(), Vec::new()),
        };
        let alternatives = [self.alternatives, alternative].concat();
        Self {
            alternatives,
            ..new_self
        }
    }
}
