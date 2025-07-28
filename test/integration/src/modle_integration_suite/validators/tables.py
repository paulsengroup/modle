# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from typing import Dict

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal


def compare_bed6(
    expected: pd.DataFrame,
    found: pd.DataFrame,
    rtol: float = 1.0e-5,
) -> Dict[str, str]:
    assert 0 <= rtol <= 1.0
    if expected.columns.tolist() != found.columns.tolist():
        return {"column name mismatch": f"expected {expected.columns.tolist()}, found {found.columns.tolist()}"}

    if (expected.dtypes != found.dtypes).any():
        return {"column data type mismatch": f"expected {expected.dtypes}, found {found.dtypes}"}

    if len(expected) != len(found):
        return {"record number mismatch": f"expected {len(expected)} records, found {len(found)}"}

    columns = expected.columns.tolist()
    columns.remove("score")

    expected = expected.copy()
    found = found.copy()
    expected["type"] = "expected"
    found["type"] = "found"

    coord_differences = pd.concat([expected, found]).drop_duplicates(
        subset=columns,
        keep=False,
        ignore_index=True,
    )

    if len(coord_differences) != 0:
        return {"found differences in genomic intervals": f"found {len(coord_differences)} differences"}

    score_differences = (~np.isclose(expected["score"], found["score"], rtol=rtol, equal_nan=True)).sum()
    if score_differences != 0:
        return {"found differences in the scores": f"found {score_differences} differences"}

    return {}


def compare_tables(
    expected: pd.DataFrame,
    found: pd.DataFrame,
    rtol: float = 1.0e-5,
) -> Dict[str, str]:
    assert 0 <= rtol <= 1.0
    if expected.columns.tolist() != found.columns.tolist():
        return {"column name mismatch": f"expected {expected.columns.tolist()}, found {found.columns.tolist()}"}

    if (expected.dtypes != found.dtypes).any():
        return {"column data type mismatch": f"expected {expected.dtypes}, found {found.dtypes}"}

    if len(expected) != len(found):
        return {"record number mismatch": f"expected {len(expected)} records, found {len(found)}"}

    try:
        assert_frame_equal(expected, found, rtol=rtol)
    except AssertionError as e:
        return {"tables differ": str(e)}

    return {}
