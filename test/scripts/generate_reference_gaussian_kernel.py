#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import sys

import numpy as np


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument("kernel-size", type=int)

    cli.add_argument("kernel-sigma", type=float)

    cli.add_argument("--sep", type=lambda s: "\t" if s == "\\t" else s, default=",")

    return cli


def compute_kernel(size: int, sigma: float) -> np.ndarray:
    """
    https://stackoverflow.com/a/43346070
    """
    assert size > 0
    assert sigma > 0

    ax = np.linspace(-(size - 1) / 2.0, (size - 1) / 2.0, size)
    gauss = np.exp(-0.5 * np.square(ax) / np.square(sigma))
    kernel = np.outer(gauss, gauss)
    return kernel / np.sum(kernel)


def main():
    args = vars(make_cli().parse_args())

    sep = args["sep"]
    kernel = compute_kernel(args["kernel-size"], args["kernel-sigma"])

    np.savetxt(sys.stdout, kernel.flatten(), delimiter=sep)


if __name__ == "__main__":
    main()
