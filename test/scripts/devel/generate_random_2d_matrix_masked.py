#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import sys

import numpy as np


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument("size", type=int)
    cli.add_argument("diagonal-width", type=int)

    cli.add_argument("--seed", type=int, default=3259565483)
    cli.add_argument("--sep", type=lambda s: "\t" if s == "\\t" else s, default="\t")

    return cli


def main():
    args = vars(make_cli().parse_args())

    size = args["size"]
    diagonal_width = args["diagonal-width"]

    np.random.seed(args["seed"])

    assert size > 0
    assert diagonal_width > 0

    if diagonal_width > size:
        print(
            "Matrix size is smaller than the specifiad diagonal width, using size as diagonal-width...", file=sys.stderr
        )
        diagonal_width = size

    ## Generate a symmetric matrix of random integers
    print(
        f"Generating a matrix with the following shape: {size}x{size} ({((size ** 2) * 8) / 1e9} GB)", file=sys.stderr
    )
    m = np.random.randint(0, 2**32, (size, size))

    ## Mask elements that are far from the diagonal
    print(f"Making elements that are more than {diagonal_width / 2} cells away from the diagonal...", file=sys.stderr)
    for i in range(0, len(m)):
        start = max(0, i - (diagonal_width / 2))
        end = min(len(m[i]), i + (diagonal_width / 2))
        for j in range(0, len(m[i])):
            if j < start or j > end:
                m[i][j] = 0

    print(f"Making matrix symmetric...", file=sys.stderr)
    m = np.tril(m) + np.tril(m, -1).T

    # Output result as TSV
    print(f"Writing result to stdout...", file=sys.stderr)
    np.savetxt(sys.stdout, m, delimiter=args["sep"], fmt="%d")


if __name__ == "__main__":
    main()
