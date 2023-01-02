#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import sys

import numpy as np
from scipy import signal


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument("kernel-radius", type=int)

    cli.add_argument("sigma", type=float)

    cli.add_argument("--sep", type=lambda s: "\t" if s == "\\t" else s, default=",")

    return cli


def compute_2d_kernel(radius: int, std: float) -> np.ndarray:
    """
    Source: https://stackoverflow.com/a/46892763
    """
    assert radius > 0
    assert std >= 0

    size = 2 * radius + 1

    kernel1d = signal.gaussian(size, std=std).reshape(size, 1)
    kernel2d = np.outer(kernel1d, kernel1d)

    return kernel2d / kernel2d.sum()


def main():
    args = vars(make_cli().parse_args())

    sep = args["sep"]
    kernel = compute_2d_kernel(args["kernel-radius"], args["sigma"])

    np.savetxt(sys.stdout, kernel.flatten(), delimiter=sep)


if __name__ == "__main__":
    main()
