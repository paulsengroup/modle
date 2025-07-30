# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
from typing import Dict

import numpy as np
import pyBigWig

from .common import compare_chroms


def compare_bigwigs(
    expected: pathlib.Path,
    found: pathlib.Path,
    rtol: float = 1.0e-5,
) -> Dict[str, str]:
    assert expected.is_file()
    if not found.is_file():
        return {"bigwig file is missing": ""}

    with pyBigWig.open(str(expected)) as bw1, pyBigWig.open(str(found)) as bw2:
        errors = compare_chroms(bw1.chroms(), bw2.chroms())
        if len(errors) != 0:
            return errors

        for chrom, size in bw1.chroms().items():
            x1 = bw1.intervals(chrom, 0, size)
            x2 = bw2.intervals(chrom, 0, size)

            if x1 is None:
                x1 = np.empty(0, dtype=np.float32)

            if x2 is None:
                x2 = np.empty(0, dtype=np.float32)

            if len(x1) != len(x2):
                errors |= {f"{chrom}: unexpected number of entries": f"expected {len(x1)}, found {len(x2)}"}
                continue

            differences = (~np.isclose(x1, x2, rtol=rtol, equal_nan=True)).sum()
            if differences != 0:
                errors |= {f"{chrom}: found differences in pixel counts": f"found {differences} differences"}

        return errors
