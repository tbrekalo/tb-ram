"""Common variables and utility functions shared accross noteboooks."""

from pathlib import Path

import polars as pl

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


def create_annotated_overlaps_from_chains(
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
    return (
        create_overlaps_from_chains(df_chains)
        .join(
            df_origins,
            on="query-name",
            how="inner",
        )
        .with_columns(
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
    )
