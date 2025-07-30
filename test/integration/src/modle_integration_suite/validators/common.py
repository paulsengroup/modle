# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from typing import Dict


def compare_chroms(expected: Dict[str, int], found: Dict[str, int]) -> Dict[str, str]:
    if list(expected.keys()) != list(found.keys()):
        return {"chromosome names mismatch": f"expected {list(expected.keys())}, found {list(found.keys())}"}

    if expected != found:
        return {"chromosome sizes mismatch": f"expected {expected}, found {found}"}

    return {}
