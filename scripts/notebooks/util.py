"""Common variables and utility functions shared accross noteboooks and scripts."""

from collections import namedtuple
from pathlib import Path

import pandas as pd
import polars as pl
from sklearn import metrics

MATCH_COLUMNS = [
    "query-name",
    "query-length",
    "query-match",
    "strand",
    "target-name",
    "target-length",
    "target-match",
]

MATCH_CHAIN_COLUMS = [
    *MATCH_COLUMNS,
    "query-matches",
    "target-matches",
    "chain-id",
]

PAF_COLUMNS = [
    "query-name",
    "query-length",
    "query-start",
    "query-end",
    "strand",
    "target-name",
    "target-length",
    "target-start",
    "target-end",
    "n-residue-matches",
]


def load_paf_df(path: Path | str) -> pl.DataFrame:
    """Loads overlaps in paf format.

    Args:
        path: A path to overlaps in paf format.

    Returns:
        A polars.DataFrame with PAF_COLUMNS.
    """
    return pl.read_csv(
        path,
        columns=list(
            range(
                len(PAF_COLUMNS),
            )
        ),
        new_columns=PAF_COLUMNS,
        separator="\t",
    )


def load_chains_df(path: Path | str) -> pl.DataFrame:
    """Loads ram chains file.

    Extends loaded DataFrame with chain length info and
    overlap start, end locations.

    Args:
        path: A path to chains file.

    Returns:
        A polars.DataFrame with MATCH_CHAIN_COLUMNS.
    """
    df_chains = pl.read_csv(
        path,
        new_columns=MATCH_CHAIN_COLUMS,
        separator="\t",
    )

    return df_chains.join(
        df_chains.group_by(pl.col("chain-id")).agg(pl.len().alias("chain-length")),
        on="chain-id",
        how="left",
    )


def load_origins_df(path: Path | str) -> pl.DataFrame:
    """Create origins pl.DataFrame.

    This pl.DataFrame contains information about read's
    position on the reverence.

    Args:
        path: A path to csv file containing origin information.
    """
    return pl.read_csv(path).select(
        pl.col("read").alias("query-name"),
        pl.col("strand").alias("origin-strand"),
        pl.col("start").alias("origin-start"),
        pl.col("end").alias("origin-end"),
    )


def create_overlaps_from_chains(df_chains: pl.DataFrame) -> pl.DataFrame:
    """Create overlaps pl.DataFrame from chains.

    Args:
        df_chains: loaded ram chains pl.DataFrame.

    Returns:
        A pl.DataFrame with columns:
            chain-id,
            query-name, query-start, query-end,
            strand,
            target-name, target-start, target-end
    """
    return df_chains.group_by(
        "chain-id",
        maintain_order=True,
    ).agg(
        pl.col("query-name").first(),
        pl.col("query-match").first().alias("query-start"),
        pl.col("query-match").last().alias("query-end"),
        pl.col("strand").first(),
        pl.col("target-name").first(),
        pl.when(pl.col("strand") == "+")
        .then(pl.col("target-match").first())
        .otherwise(pl.col("target-match").last())
        .alias("target-start")
        .first(),
        pl.when(pl.col("strand") == "+")
        .then(pl.col("target-match").last())
        .otherwise(pl.col("target-match").first())
        .alias("target-end")
        .first(),
    )


def create_annotated_ref_overlaps(
    df_ref_overlaps: pl.DataFrame,
    df_origins: pl.DataFrame,
) -> pl.DataFrame:
    """Creates a pl.DataFrame with annoated overlaps.
    
    Args:
        df_ref_overlaps: A pl.DataFrame representing reads to ref overlaps.
        df_origins: A pl.DataFrame with ground truth information
                    about read reference origins.

    Retruns:
        A pl.DataFrame with standard PAF columns and aditional origin info.
    """
    return df_ref_overlaps.join(
        df_origins,
        on="query-name",
        how="inner",
    ).with_columns(
        (
            pl.max_horizontal(
                pl.min_horizontal(
                    pl.col("origin-end"),
                    pl.col("target-end"),
                )
                - pl.max_horizontal(
                    pl.col("origin-start"),
                    pl.col("target-start"),
                ),
                pl.lit(0),
            )
            / pl.max_horizontal(
                pl.col("origin-end") - pl.col("origin-start"),
                pl.col("target-end") - pl.col("target-start"),
            )
        ).alias("ratio")
    )


def create_annotated_ref_overlaps_from_chains(
    df_chains: pl.DataFrame,
    df_origins: pl.DataFrame,
) -> pl.DataFrame:
    """Creates a an overlap pl.DataFrame from chains pl.DataFrame.

    Created pl.DataFrame is joined with read origins info. We also
    calculate the reference overlap between origin information and
    overlaps info.

    Returns:
        A pl.DataFrame with annotated overlaps.
    """
    return create_annotated_ref_overlaps(
        df_ref_overlaps=create_overlaps_from_chains(df_chains),
        df_origins=df_origins,
    )


def expand_ava_with_origin_info(
    df_ava: pl.DataFrame,
    df_origins: pl.DataFrame,
) -> pl.DataFrame:
    """Expand mapping pl.DataFrame with origin information.

    Args:
        df_ava: A pl.DataFrame with `query-name` and `target-name` columns.
        df_origins: A pl.DataFrame with read reference origin information.

    Returns:
        A pl.DataFrame with query and target origin information.
    """
    return df_ava.join(
        df_origins.select(
            pl.col("query-name"),
            pl.col("origin-strand").alias("query-origin-strand"),
            pl.col("origin-start").alias("query-origin-start"),
            pl.col("origin-end").alias("query-origin-end"),
        ),
        on="query-name",
    ).join(
        df_origins.select(
            pl.col("query-name").alias("target-name"),
            pl.col("origin-strand").alias("target-origin-strand"),
            pl.col("origin-start").alias("target-origin-start"),
            pl.col("origin-end").alias("target-origin-end"),
        ),
        on="target-name",
    )


def calc_ava_origin_overlap(df_origin_annotated: pl.DataFrame) -> pl.DataFrame:
    """Given a pl.DataFrame with target and query origin information,
    calculate their overlap and strand match on the origin reference.

    Args:
        pl_origin_annotated: A pl.DataFrame with query and target origin info.

    Returns:
        A pl.DataFrame with overlap ratio on origin reference and strand match.
    """

    return df_origin_annotated.with_columns(
        (
            pl.max_horizontal(
                pl.min_horizontal(
                    pl.col("query-origin-end"),
                    pl.col("target-origin-end"),
                )
                - pl.max_horizontal(
                    pl.col("query-origin-start"), pl.col("target-origin-start")
                ),
                pl.lit(0),
            )
            / pl.max_horizontal(
                pl.col("query-origin-end") - pl.col("query-origin-start"),
                pl.col("target-origin-end") - pl.col("query-origin-start"),
            )
        ).alias("ratio"),
        (
            pl.when(pl.col("strand") == "+")
            .then(pl.col("query-origin-strand") == pl.col("target-origin-strand"))
            .otherwise(pl.col("query-origin-strand") != pl.col("target-origin-strand"))
        ).alias("matching-strands"),
    )


def create_annotated_overlaps_from_ava(
    df_overlaps: pl.DataFrame,
    df_origins: pl.DataFrame,
    min_ratio: float = 0.875,
) -> pl.DataFrame:
    """Given all-vs-all overlaps and read reference origin information
    construct a pl.DataFrame with true/false positive annotations.

    Args:
        df_overlaps: A pl.DataFrame with read to read overlaps.
        df_origins: A pl.DataFrame containing read reference origins.

    Returns:
        A pl.DataFrame with original overlaps along side with true/false
        positive annotation.
    """

    return calc_ava_origin_overlap(
        expand_ava_with_origin_info(
            df_overlaps,
            df_origins,
        )
    ).select(
        pl.col("query-name"),
        pl.col("query-length"),
        pl.col("query-start"),
        pl.col("query-end"),
        pl.when(pl.col("strand") == "+").then(1).otherwise(0).alias("strand"),
        pl.col("target-name"),
        pl.col("target-length"),
        pl.col("target-start"),
        pl.col("target-end"),
        pl.col("n-residue-matches"),
        ((pl.col("ratio") > min_ratio) & pl.col("matching-strands"))
        .cast(pl.Int64)
        .alias("label"),
    )


TrainTestData = namedtuple(
    "TrainTestData",
    [
        "X_train",
        "X_test",
        "y_train",
        "y_test",
    ],
)


def concat_report_dicts_to_df(reports: dict[str, dict]) -> pl.DataFrame:
    """Concat sklearn classification report dicts into a long pl.DataFrame.

    Args:
        reports: A dict mapping models to classification report dictinaries.

    Returns:
        A long pl.DataFrame with columns: model, class, metric and value.
        The pl.DataFrame summarizes model performance.
    """
    df_report = pl.DataFrame(
        schema={"model": str, "class": pl.Int64, "metric": str, "value": pl.Float64}
    )
    for name, report in reports.items():
        df_report = pl.concat(
            [
                df_report,
                pl.DataFrame(
                    pd.melt(
                        pd.DataFrame(report)[["0", "1"]]
                        .transpose()[["precision", "recall", "f1-score"]]
                        .reset_index(),
                        id_vars=["index"],
                        var_name="metric",
                    )[["index", "metric", "value"]]
                ).select(
                    pl.lit(name).alias("model"),
                    pl.col("index").alias("class").cast(pl.Int64),
                    pl.col("metric"),
                    pl.col("value"),
                ),
            ]
        )

    return df_report


def create_df_reports(models: dict, data: TrainTestData) -> pl.DataFrame:
    """Calculates accuracy, precision and recall for models on given data.

    Args:
        models: A dictinary mapping model names to models.
        data: A TrainTestData instance used for evaluating performance.

    Returns:
        A long pl.DataFrame with columns: model, class, metric and value.
        The pl.DataFrame summarizes model performance.
    """
    return concat_report_dicts_to_df(
        {
            name: metrics.classification_report(
                data.y_test, model.predict(data.X_test), output_dict=True
            )
            for name, model in models.items()
        }
    )
