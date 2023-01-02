#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import numbers
from typing import Tuple

import numpy as np
from rpy2.robjects import default_converter, numpy2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from scipy.stats import pearsonr, spearmanr


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument("iterations", type=int)
    cli.add_argument("vector-size", type=int)
    cli.add_argument("lower-bound", type=float)
    cli.add_argument("upper-bound", type=float)

    cli.add_argument("type", type=str, choices={"int", "uint", "float"})

    cli.add_argument("--seed", type=int, default=2572287723)
    cli.add_argument("--sep", type=lambda s: "\t" if s == "\\t" else s, default=",")

    return cli


def get_np_type(type_: str):
    if type_ == "int":
        return int

    if type_ == "uint":
        return np.uint

    assert type_ == "float"
    return float


def generate_random_vectors(size: int, lb, ub, dtype) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    assert size > 0
    assert lb < ub

    if isinstance(dtype, str):
        dtype = get_np_type(dtype)

    if dtype == np.uint:
        assert lb >= 0

    lb = dtype(lb)
    ub = dtype(ub)

    weights = np.random.uniform(low=0, high=1, size=size)

    if dtype != float:
        return (
            np.random.randint(low=lb, high=ub, size=size, dtype=int),
            np.random.randint(low=lb, high=ub, size=size, dtype=int),
            weights,
        )

    return (
        np.random.uniform(low=lb, high=ub, size=size),
        np.random.uniform(low=lb, high=ub, size=size),
        weights,
    )


def ndarray_to_tsv_str(v: np.ndarray) -> str:
    if issubclass(v.dtype.type, numbers.Integral):
        return "\t".join((f"{n}" for n in v))

    return "\t".join((f"{n:.18g}" for n in v))


def call_wcorr(v1, v2, w) -> Tuple[float, float, float, float]:
    with localconverter(np_cv_rules) as cv:
        pcc = wcorr.weightedCorr(v1, v2, "Pearson", w)[0]
        rho = wcorr.weightedCorr(v1, v2, "Spearman", w)[0]

        return pcc, np.nan, rho, np.nan


def correlate(v1: np.ndarray, v2: np.ndarray, weights: np.ndarray) -> str:
    buff = list(pearsonr(v1, v2))
    buff.extend(list(spearmanr(v1, v2)))
    buff.extend(call_wcorr(v1, v2, weights))

    return "|".join((f"{n:.18g}" for n in buff))


def main():
    args = vars(make_cli().parse_args())

    assert args["iterations"] > 0
    assert args["vector-size"] > 0

    np.random.seed(args["seed"])

    # Every iteration a line with the following syntax is printed to stdout:
    # |v1|v2|weights|pcc|pcc_pval|rho|rho_pval|pccw|pccw_pval|rhow|rhow_pval|
    # v1, v2 and weights are comma-separated list of numbers (ints or floats)
    for _ in range(args["iterations"]):
        v1, v2, weights = generate_random_vectors(
            args["vector-size"], args["lower-bound"], args["upper-bound"], args["type"]
        )

        cfx = correlate(v1, v2, weights)

        v1 = ndarray_to_tsv_str(v1)
        v2 = ndarray_to_tsv_str(v2)
        weights = ndarray_to_tsv_str(weights)

        print(f"{v1}|{v2}|{weights}|{cfx}")


if __name__ == "__main__":
    np_cv_rules = default_converter + numpy2ri.converter
    wcorr = importr("wCorr")

    main()
