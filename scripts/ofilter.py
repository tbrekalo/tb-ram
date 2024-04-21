"""Script for filtering ram format overlaps.
"""

import sys
import logging
from argparse import ArgumentParser
from pathlib import Path

import catboost as cb
import numpy as np
import polars as pl

import util

logging.basicConfig(
    format="time=%(asctime)s loc=%(filename)s:%(lineno)s level=%(levelname)s %(message)s",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

RAM_OVLP_COLS = [
    "query-name",
    "query-length",
    "query-begin",
    "query-end",
    "query-matches",
    "query-gaps",
    "strand",
    "target-name",
    "target-length",
    "target-begin",
    "target-end",
    "target-matches",
    "target-gaps",
    "diff-mean",
    "q75",
    "q90",
    "q95",
    "q98",
]


def create_cli_parser() -> ArgumentParser:
    """Creates a command line argument parser."""
    parser = ArgumentParser(
        "ofilter",
        description="script for filtering ram format overlaps with ML models.",
    )

    parser.add_argument(
        "-m",
        "--model-path",
        required=True,
        help="Path to CatBoostClassifier model.",
    )
    parser.add_argument(
        "-i",
        "--input-path",
        required=True,
        help="Path to ram formatted overlaps.",
    )

    return parser


def main(input_path: str | Path, model_path: str | Path) -> int:
    """Program main entry point.

    Args:
        input_path: A path to ram formatted overlaps.
        model_path: A path to CatBoostClassifier model.

    Returns:
        0 on a successful run, otherwise an error code.
    """
    logger.info(
        "input_path=%s model_path=%s",
        input_path,
        model_path,
    )

    df_input_overlaps = pl.read_csv(
        input_path,
        has_header=False,
        new_columns=RAM_OVLP_COLS,
        separator="\t",
    )
    model = cb.CatBoostClassifier().load_model(model_path)
    predictions = model.predict(
        df_input_overlaps.select(
            pl.col("diff-mean"),
            *[pl.col(f"q{q}") for q in [75, 90, 95, 98]],
            (pl.col("query-end") - pl.col("query-begin")).alias("query-overlap-length"),
            pl.col("query-matches"),
            (pl.col("target-end") - pl.col("target-begin")).alias(
                "target-overlap-length"
            ),
            pl.col("target-matches"),
        ).to_pandas()
    )
    indices = np.where(predictions == 1)[0]
    logging.info("selection_rate=%f", len(indices) / len(predictions))

    df_out = df_input_overlaps[indices].select(
        pl.col("query-name"),
        pl.col("query-length"),
        pl.col("query-begin"),
        pl.col("query-end"),
        pl.col("strand"),
        pl.col("target-name"),
        pl.col("target-length"),
        pl.col("target-begin"),
        pl.col("target-end"),
        pl.min_horizontal(pl.col("query-matches"), pl.col("target-matches")).alias(
            "score"
        ),
        pl.max_horizontal(
            pl.col("query-end") - pl.col("query-begin"),
            pl.col("target-end") - pl.col("target-begin"),
        ).alias("alignment-block-length"),
        pl.lit(255).alias("mapping-quality"),
    )
    df_out.write_csv(file=sys.stdout, include_header=False, separator="\t")


if __name__ == "__main__":
    parser = create_cli_parser()
    try:
        sys.exit(main(**vars(parser.parse_args())))
    except Exception as e:
        logger.exception(e)
        sys.exit(1)
