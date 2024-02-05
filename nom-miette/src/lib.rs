use std::{collections::HashMap, fmt};

use miette::{Diagnostic, LabeledSpan, SourceSpan};
use nom::{
    combinator::{all_consuming, complete, consumed},
    error::ParseError,
    Err, Finish, IResult, Parser,
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

// FIXME: Check that field ordering everywhere matches this!
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct LabeledParseError<'a, E> {
    input: &'a str,
    length: usize,
    kind: E,
    // FIXME: Replace this with the ErrorTree!!
    alternatives: Vec<LabeledParseError<'a, E>>,
    source: Option<Box<LabeledParseError<'a, E>>>,
}

pub trait LabeledErrorKind: Diagnostic + Clone + Eq + From<nom::error::ErrorKind> {
    fn label(&self) -> Option<&'static str> {
        None
    }
}

pub trait FromExternalError<'a, E> {
    const FATAL: bool = false;
    fn from_external_error(input: &'a str, error: E) -> Self;
}

// FIXME: Oh lord, what a mess...
impl<'a, E: LabeledErrorKind> LabeledParseError<'a, E> {
    // FIXME: Make sure this is used everywhere it can be!
    pub fn new(input: &'a str, kind: E) -> Self {
        Self::new_with_source(input, kind, None)
    }

    // FIXME: Make sure this is used everywhere it can be!
    pub fn new_with_source(
        input: &'a str,
        kind: E,
        source: Option<LabeledParseError<'a, E>>,
    ) -> Self {
        Self {
            input,
            length: 0,
            kind,
            alternatives: Vec::new(),
            source: source.map(Box::new),
        }
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
        let span = self.span_from_input(full_input);
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

    // FIXME: OMG refactor... Also, can I do better than stealing here?
    fn span_from_input(&self, full_input: &str) -> SourceSpan {
        let base_addr = full_input.as_ptr() as usize;
        let substr_addr = self.input.as_ptr() as usize;
        // FIXME: Keep this?
        assert!(
            substr_addr >= base_addr,
            "tried to get the span of a non-substring!"
        );
        let start = substr_addr - base_addr;
        let end = start + self.length;
        SourceSpan::from(start..end)
    }
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

// FIXME: Why are these generics so much messier than map_res from nom?
pub fn map_res<'a, O1, O2, E1, E2, F, G>(
    parser: F,
    mut f: G,
) -> impl FnMut(&'a str) -> IResult<&'a str, O2, LabeledParseError<'a, E1>>
where
    O1: Clone,
    E1: LabeledErrorKind,
    F: Copy + Parser<&'a str, O1, LabeledParseError<'a, E1>>,
    G: Copy + FnMut(O1) -> Result<O2, E2>,
    LabeledParseError<'a, E1>: FromExternalError<'a, E2>,
{
    move |input| {
        let i = input;
        let (input, (consumed, o1)) = consumed(parser)(input)?;
        match f(o1) {
            Ok(o2) => Ok((input, o2)),
            Err(e) => {
                let e = LabeledParseError {
                    length: consumed.len(),
                    ..LabeledParseError::from_external_error(i, e)
                };
                Err(if LabeledParseError::FATAL {
                    Err::Failure(e)
                } else {
                    Err::Error(e)
                })
            }
        }
    }
}

// FIXME: Check if I'm being consistent about using `impl` or generics... I think I should avoid any generics I don't
// use, as long as this is just a sort of "internal" library
// FIXME: See if this signature can be simplified (elide lifetimes?)
// FIXME: Standardize the order of all of these generic arguments!
pub fn wrap_err<'a, O, P, E>(
    mut parser: P,
    kind: E,
) -> impl FnMut(&'a str) -> IResult<&'a str, O, LabeledParseError<'a, E>>
where
    // FIXME: Eek, this is a where, and below (in expect) is not!
    E: LabeledErrorKind + Clone,
    P: Parser<&'a str, O, LabeledParseError<'a, E>>,
{
    // FIXME: DRY with expect below!
    move |i| {
        parser
            .parse(i)
            // FIXME: Really don't get why this clone is here...
            .map_err(|e| e.map(|e| LabeledParseError::new_with_source(i, kind.clone(), Some(e))))
    }
}

// FIXME: I should just replace wrap_err with expect, or I should make expect *not* wrap things
pub fn expect<'a, O, E: LabeledErrorKind + Clone, F>(
    mut parser: F,
    kind: E,
) -> impl FnMut(&'a str) -> IResult<&'a str, O, LabeledParseError<'a, E>>
where
    F: Parser<&'a str, O, LabeledParseError<'a, E>>,
{
    move |i| {
        parser
            .parse(i)
            // FIXME: Really don't get why this clone is here...
            .map_err(|e| e.map(|_| LabeledParseError::new(i, kind.clone())))
    }
}

// FIXME: Eventually, I should make everything generic over the input type again... So you'd be able to use this
// library with &'a [u8] like `nom` lets you
impl<'a, E: LabeledErrorKind> ParseError<&'a str> for LabeledParseError<'a, E> {
    fn from_error_kind(input: &'a str, kind: nom::error::ErrorKind) -> Self {
        Self::new(input, kind.into())
    }

    fn append(_input: &str, _kind: nom::error::ErrorKind, other: Self) -> Self {
        other
    }

    fn or(self, other: Self) -> Self {
        Self {
            alternatives: [self.alternatives, vec![other]].concat(),
            ..self
        }
    }
}
