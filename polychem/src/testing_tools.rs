use miette::Diagnostic;
use std::fmt::Debug;

pub trait UnwrapDiagnostic<E> {
    fn unwrap_diagnostic(self) -> E;
}

impl<T: Debug, E: Diagnostic> UnwrapDiagnostic<E> for Result<T, E> {
    fn unwrap_diagnostic(self) -> E {
        self.unwrap_err()
    }
}

impl<T: Debug, E: Diagnostic> UnwrapDiagnostic<E> for Result<T, Box<E>> {
    fn unwrap_diagnostic(self) -> E {
        *self.unwrap_err()
    }
}

macro_rules! assert_miette_snapshot {
    ($diag:expr) => {{
        use insta::{with_settings, assert_snapshot};
        use miette::{GraphicalReportHandler, GraphicalTheme};
        use crate::testing_tools::UnwrapDiagnostic;

        let mut out = String::new();
        GraphicalReportHandler::new_themed(GraphicalTheme::unicode_nocolor())
            .with_width(80)
            .render_report(&mut out, &$diag.unwrap_diagnostic())
            .unwrap();
        with_settings!({
            description => stringify!($diag)
        }, {
            assert_snapshot!(out);
        });
    }};
}

pub(crate) use assert_miette_snapshot;
