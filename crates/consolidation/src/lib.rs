use std::io::Cursor;

use polars::prelude::*;

// Constants ===========================================================================================================

// NOTE: Z-value for 95% confidence interval
const Z_SCORE_FOR_CI_95: f64 = 1.96;
const OLIGO_STATE_REGEX: &str = r"\|(\d+)";

struct InputColumns;
impl InputColumns {
    const ALL: [&str; 5] = [
        Self::STRUCTURE,
        Self::INTENSITY,
        Self::RT,
        Self::THEO,
        Self::PPM,
    ];
    const STRUCTURE: &str = "Structure";
    const INTENSITY: &str = "Consolidated Intensity";
    const RT: &str = "Consolidated RT (min)";
    const THEO: &str = "Consolidated Theo (Da)";
    const PPM: &str = "Consolidated Delta (ppm)";
}

struct TemporaryColumns;
impl TemporaryColumns {
    const REPLICATE: &str = "Replicate";
    const AVG_INTENSITY: &str = "Avg Intensity";
    const CI_INTENSITY: &str = "95% CI Intensity";
}

struct OutputColumns;
impl OutputColumns {
    const STRUCTURE: &str = "Structure";
    const OLIGO_STATE: &str = "Oligo State";
    const AVG_ABUNDANCE: &str = "Average Abundance (%)";
    const CI_ABUNDANCE: &str = "95% CI Abundance (%)";
    const AVG_RT: &str = "Average RT (min)";
    const CI_RT: &str = "95% CI RT (min)";
    const THEO: &str = "Theo (Da)";
    const PPM: &str = "Absolute Delta (ppm)";
    const PRESENT_IN: &str = "Present In";
    const REPLICATE_INTENSITY: &str = "Intensity";
    const REPLICATE_ABUNDANCE: &str = "Abundance";
}

// Public API ==========================================================================================================

pub type ReplicateNumber = u32;
pub type Replicate = (ReplicateNumber, String);

pub fn consolidate_replicates(replicates: &[Replicate]) -> String {
    // SAFETY: Overall this unwrap should be okay, since I also control PGFinder's outputs
    load_replicates(replicates)
        .map(consolidate)
        .and_then(into_csv)
        .unwrap_or_default()
}

// Private Types =======================================================================================================

struct ReplicateFrame(Vec<ReplicateNumber>, LazyFrame);

// Private Functions ===================================================================================================

// TODO: Use a fixed-point Decimal type instead of a float, so that we don't have nasty floating-point errors that show
// up in all of the results...
fn load_replicates(replicates: &[Replicate]) -> PolarsResult<ReplicateFrame> {
    let replicate_numbers: Vec<_> = replicates.iter().map(|&(n, _)| n).collect();
    let lazyframe = concat(
        replicates
            .iter()
            .map(|&(replicate, ref csv)| {
                Ok(CsvReader::new(Cursor::new(csv))
                    .finish()?
                    .lazy()
                    .select(InputColumns::ALL.map(col))
                    .drop_nulls(None)
                    .with_column(lit(replicate).alias(TemporaryColumns::REPLICATE)))
            })
            .collect::<PolarsResult<Vec<_>>>()?,
        UnionArgs::default(),
    );
    Ok(ReplicateFrame(replicate_numbers, lazyframe?))
}

// FIXME: Benchmark if taking in the replicate numbers is even worth it... I'm just trying to not break the chain of
// lazy operations here, so I'm using `group_by` instead of just `pivot`...
fn consolidate(ReplicateFrame(replicates, df): ReplicateFrame) -> LazyFrame {
    let individual_replicate_columns = df
        .clone()
        .group_by([InputColumns::STRUCTURE])
        .agg(wide_formatted_replicate_columns(&replicates))
        .with_columns(replicate_abundance_columns(&replicates));

    df.group_by([InputColumns::STRUCTURE])
        .agg(consolidated_columns())
        .with_column(oligo_state_column(OLIGO_STATE_REGEX))
        .with_columns(average_abundance_columns())
        .join(
            individual_replicate_columns,
            [col(OutputColumns::STRUCTURE)],
            [col(InputColumns::STRUCTURE)],
            JoinArgs::default(),
        )
        .select(formatted_output_columns(&replicates))
        .sort(
            [
                OutputColumns::OLIGO_STATE,
                OutputColumns::AVG_ABUNDANCE,
                OutputColumns::THEO,
            ],
            SortMultipleOptions::new().with_order_descending_multi([false, true, false]),
        )
}

fn into_csv(df: LazyFrame) -> PolarsResult<String> {
    let mut result = Vec::new();
    CsvWriter::new(&mut result).finish(&mut df.collect()?)?;
    // SAFETY: The `CsvWriter` should always return valid UTF-8
    Ok(String::from_utf8(result).unwrap())
}

// ---------------------------------------------------------------------------------------------------------------------

fn wide_formatted_replicate_columns(replicates: &[u32]) -> Vec<Expr> {
    [
        vec![len().alias(OutputColumns::PRESENT_IN)],
        replicates
            .iter()
            .map(|&n| {
                col(InputColumns::INTENSITY)
                    .filter(col(TemporaryColumns::REPLICATE).eq(n))
                    .first()
                    .alias(replicate_column(n, OutputColumns::REPLICATE_INTENSITY))
            })
            .collect(),
    ]
    .concat()
}

fn replicate_abundance_columns(replicates: &[u32]) -> Vec<Expr> {
    replicates
        .iter()
        .map(|&n| {
            let replicate_intensity_column =
                replicate_column(n, OutputColumns::REPLICATE_INTENSITY);
            intensity_to_abundance(&replicate_intensity_column, &replicate_intensity_column)
                .alias(replicate_column(n, OutputColumns::REPLICATE_ABUNDANCE))
        })
        .collect()
}

fn consolidated_columns() -> [Expr; 6] {
    let confidence_interval = |c| lit(Z_SCORE_FOR_CI_95) * col(c).std(1) / len().sqrt();
    [
        mean(InputColumns::INTENSITY).alias(TemporaryColumns::AVG_INTENSITY),
        confidence_interval(InputColumns::INTENSITY).alias(TemporaryColumns::CI_INTENSITY),
        mean(InputColumns::RT).alias(OutputColumns::AVG_RT),
        confidence_interval(InputColumns::RT).alias(OutputColumns::CI_RT),
        col(InputColumns::THEO).first().alias(OutputColumns::THEO),
        col(InputColumns::PPM)
            .abs()
            .mean()
            .alias(OutputColumns::PPM),
    ]
}

fn oligo_state_column(oligo_state_re: &str) -> Expr {
    col(OutputColumns::STRUCTURE)
        .str()
        .extract(lit(oligo_state_re), 1)
        .alias(OutputColumns::OLIGO_STATE)
}

fn average_abundance_columns() -> [Expr; 2] {
    [
        intensity_to_abundance(
            TemporaryColumns::AVG_INTENSITY,
            TemporaryColumns::AVG_INTENSITY,
        )
        .alias(OutputColumns::AVG_ABUNDANCE),
        intensity_to_abundance(
            TemporaryColumns::CI_INTENSITY,
            TemporaryColumns::AVG_INTENSITY,
        )
        .alias(OutputColumns::CI_ABUNDANCE),
    ]
}

fn formatted_output_columns(replicates: &[u32]) -> Vec<Expr> {
    [
        vec![
            col(OutputColumns::STRUCTURE)
                .str()
                .replace_all(lit(OLIGO_STATE_REGEX), lit(""), false),
            col(OutputColumns::OLIGO_STATE),
            col(OutputColumns::AVG_ABUNDANCE).round(3),
            col(OutputColumns::CI_ABUNDANCE).round(3),
            col(OutputColumns::AVG_RT).round(2),
            col(OutputColumns::CI_RT).round(2),
            col(OutputColumns::THEO).round(4),
            col(OutputColumns::PPM).round(1),
            col(OutputColumns::PRESENT_IN),
        ],
        replicates
            .iter()
            .flat_map(|&n| {
                [
                    col(replicate_column(n, OutputColumns::REPLICATE_INTENSITY)),
                    col(replicate_column(n, OutputColumns::REPLICATE_ABUNDANCE)).round(3),
                ]
            })
            .collect(),
    ]
    .concat()
}

// ---------------------------------------------------------------------------------------------------------------------

fn replicate_column(replicate: u32, column: &str) -> String {
    format!("Replicate {replicate} ({column})")
}

fn intensity_to_abundance(column: &str, total_column: &str) -> Expr {
    lit(100) * col(column) / sum(total_column)
}

// Unit Tests ==========================================================================================================

#[cfg(test)]
mod tests {
    use std::sync::LazyLock;

    use insta::assert_snapshot;

    use super::*;

    static REPLICATES: LazyLock<[Replicate; 3]> = LazyLock::new(|| {
        [
            (1, include_str!("../tests/data/E. coli WT_01.csv")),
            (2, include_str!("../tests/data/E. coli WT_02.csv")),
            (3, include_str!("../tests/data/E. coli WT_03.csv")),
        ]
        .map(|(n, c)| (n, c.to_owned()))
    });

    #[test]
    fn test_load_replicates() {
        assert_snapshot!(
            load_replicates(REPLICATES.as_ref())
                .and_then(|ReplicateFrame(_, df)| into_csv(df))
                .unwrap()
        );
    }

    #[test]
    fn test_consolidate() {
        assert_snapshot!(
            load_replicates(REPLICATES.as_ref())
                .map(consolidate)
                .and_then(into_csv)
                .unwrap()
        );
    }

    #[test]
    fn test_consolidate_replicates() {
        assert_snapshot!(consolidate_replicates(REPLICATES.as_ref()));
    }
}
