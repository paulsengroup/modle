#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import sys

import cooler
import numpy as np


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument("path-to-cooler", type=str)
    cli.add_argument("coords", type=str)
    return cli


def main():
    args = vars(make_cli().parse_args())

    c = cooler.Cooler(args["path-to-cooler"])

    selector = c.matrix(balance=False, sparse=False)

    m = selector.fetch(args["coords"])

    if np.issubdtype(m.dtype, np.integer):
        fmt = "%d"
    else:
        fmt = "%.18e"

    np.savetxt(sys.stdout, m, delimiter="\t", fmt=fmt)


if __name__ == "__main__":
    main()
