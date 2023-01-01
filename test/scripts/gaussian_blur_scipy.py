#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import sys

import numpy as np
from scipy.ndimage import gaussian_filter


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument("sigma", type=float)

    cli.add_argument("--cutoff", type=float, default=3.5)

    cli.add_argument("--sep", type=lambda s: "\t" if s == "\\t" else s, default="\t")

    cli.add_argument("--extension-mode", type=str, default="nearest")

    return cli


def main():
    args = vars(make_cli().parse_args())

    sigma = args["sigma"]
    sep = args["sep"]

    m = np.loadtxt(sys.stdin, dtype=float, delimiter=sep)
    m = gaussian_filter(m, sigma, mode=args["extension_mode"], truncate=args["cutoff"])

    np.savetxt(sys.stdout, m, fmt="%.18g", delimiter=sep)


if __name__ == "__main__":
    main()
