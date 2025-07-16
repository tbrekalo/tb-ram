# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "polars>=1.31.0",
# ]
# ///

from argparse import ArgumentParser
from pathlib import Path
from typing import IO

import polars as pl


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

_DataSource = str | Path | IO[str] | IO[bytes] | bytes


def _load_paf_df(path: _DataSource) -> pl.DataFrame:
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


def _load_origins_df(path: _DataSource) -> pl.DataFrame:
    return pl.read_csv(path).select(
        pl.col("read").alias("query-name"),
        pl.col("strand").alias("origin-strand"),
        pl.col("start").alias("origin-start"),
        pl.col("end").alias("origin-end"),
    )


def _create_annotated_ref_overlaps(
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


def _create_parser() -> ArgumentParser:
    parser = ArgumentParser(
        "eval",
        description="Evaluate mapper precision for mapping simulated reads onto the reference",
    )
    parser.add_argument("overlaps", type=Path)
    parser.add_argument("origins", type=Path)
    parser.add_argument(
        "-n",
        "--name",
        type=str,
        default="-",
        help="Evaluation name.",
    )
    parser.add_argument(
        "-m",
        "--min_ratio",
        type=float,
        default=0.875,
        help="How much read mapping has to overlap with the origin on the reference to be considered a true positive.",
    )
    parser.add_argument(
        "-o",
        "--out",
        type=Path,
        default=None,
        help="Csv output file.",
    )

    return parser


def _main(
    overlaps: _DataSource,
    origins: _DataSource,
    min_ratio: float,
    name: str,
    out: str | Path | None,
) -> None:
    df_overlaps = _load_paf_df(overlaps)
    df_origins = _load_origins_df(origins)

    df_annotated = (
        _create_annotated_ref_overlaps(df_overlaps, df_origins)
        .drop("origin-start", "origin-end", "origin-strand")
        .with_columns(
            pl.when(pl.col("strand") == "+")
            .then(pl.lit(1))
            .otherwise(pl.lit(0))
            .alias("strand"),
            pl.max_horizontal(
                (pl.col("query-end") - pl.col("query-start")),
                (pl.col("target-end") - pl.col("target-start")),
            ).alias("overlap-length"),
            (pl.col("ratio") > min_ratio).alias("positive"),
            pl.lit(name).alias("name"),
        )
    )

    df_stats = (
        df_annotated.group_by("name", "positive")
        .agg(
            pl.len().alias("n-overlaps"),
            pl.struct(["query-name", "target-name"]).n_unique().alias("n-unique-pairs"),
            pl.col("overlap-length").mean().alias("overlap-length-mean"),
            pl.col("overlap-length").std().alias("overlap-length-std"),
            pl.col("overlap-length").quantile(0.25).alias("overlap-length-q25"),
            pl.col("overlap-length").median().alias("overlap-length-median"),
            pl.col("overlap-length").quantile(0.75).alias("overlap-length-q75"),
        )
        .sort(by=["name", "positive"])
    )

    print(df_stats)
    if out is not None:
        df_stats.write_csv(file=out)


if __name__ == "__main__":
    parser = _create_parser()
    _main(**vars(parser.parse_args()))
