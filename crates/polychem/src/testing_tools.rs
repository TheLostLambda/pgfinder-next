macro_rules! assert_miette_snapshot {
    ($diag:expr) => {{
        use insta::{with_settings, assert_snapshot};
        use miette::{GraphicalReportHandler, GraphicalTheme};

        let mut out = String::new();
        GraphicalReportHandler::new_themed(GraphicalTheme::unicode_nocolor())
            .with_width(80)
            .render_report(&mut out, $diag.unwrap_err())
            .unwrap();
        with_settings!({
            description => stringify!($diag)
        }, {
            assert_snapshot!(out);
        });
    }};
}

pub(crate) use assert_miette_snapshot;
