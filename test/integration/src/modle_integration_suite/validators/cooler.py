# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
from typing import Dict

import hictkpy
import numpy as np
import pandas as pd

from .common import compare_chroms


def _compare_bins(
    expected: pd.DataFrame,
    found: pd.DataFrame,
) -> Dict[str, str]:
    coord_differences = pd.concat([expected, found]).drop_duplicates(
        subset=["chrom", "start", "end"], keep=False, ignore_index=True
    )

    if len(coord_differences) != 0:
        return {"found differences in bin coordinates": f"found {len(coord_differences)} differences"}

    return {}


def _compare_pixels(
    chrom: str,
    expected: pd.DataFrame,
    found: pd.DataFrame,
    rtol: float = 1.0e-5,
) -> Dict[str, str]:
    assert 0 <= rtol <= 1.0

    if (expected.dtypes != found.dtypes).any():
        return {
            f"{chrom}: pixel table has an unexpected data type": f"expected {expected.dtypes}, found {found.dtypes}"
        }

    if len(expected) != len(found):
        return {
            f"{chrom}: pixel table has an unexpected number of records": f"expected {len(expected)} records, found {len(found)}"
        }

    columns = expected.columns.tolist()
    columns.remove("count")

    expected = expected.copy()
    found = found.copy()
    expected["type"] = "expected"
    found["type"] = "found"

    coord_differences = pd.concat([expected, found]).drop_duplicates(subset=columns, keep=False, ignore_index=True)

    if len(coord_differences) != 0:
        return {f"{chrom}: found differences in pixel coordinates": f"found {len(coord_differences)} differences"}

    count_differences = (~np.isclose(expected["count"], found["count"], rtol=rtol, equal_nan=True)).sum()
    if count_differences != 0:
        return {f"{chrom}: found differences in pixel counts": f"found {count_differences} differences"}

    return {}


def compare_coolers(expected: pathlib.Path, found: pathlib.Path) -> Dict[str, str]:
    expected = hictkpy.File(expected)
    found = hictkpy.File(found)

    errors = compare_chroms(expected.chromosomes(), found.chromosomes())
    errors |= _compare_bins(expected.bins().to_df(), found.bins().to_df())

    if len(errors) != 0:
        return errors

    for chrom in expected.chromosomes():
        errors |= _compare_pixels(
            chrom,
            expected.fetch(chrom).to_df(),
            found.fetch(chrom).to_df(),
        )

    return errors
