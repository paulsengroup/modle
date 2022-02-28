#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import numpy as np
import sys

if len(sys.argv) < 3:
    print(f"Usage:   {sys.argv[0]} size diagonal width")
    print(f"Example: {sys.argv[0]} 1000 20")
    sys.exit(1)

size = int(sys.argv[1])
diagonal_width = int(sys.argv[2])

## Generate a symmetrix matrix of random integers
print(f"Generating a matrix with the following shape: {size}x{size} ({((size**2) * 8) / 1e9} GB)", file=sys.stderr)
m1 = np.random.randint(0, 2**32, (size, size))

## Mask elements that are far from the diagonal
print(f"Making elements that are more than {diagonal_width / 2} cells away from the diagonal...", file=sys.stderr)
for i in range(0, len(m1)):
    start = max(0, i - (diagonal_width / 2))
    end = min(len(m1[i]), i + (diagonal_width / 2))
    for j in range(0, len(m1[i])):
        if j < start or j > end:
            m1[i][j] = 0
    if (i % (len(m1) // 100) == 0):
        print(f"{(i / len(m1) * 100):.2f}% ...", file=sys.stderr)

print(f"Making matrix symmetric...", file=sys.stderr)
m2 = np.tril(m1) + np.tril(m1, -1).T
m1 = None

# Output result as TSV
print(f"Writing result to stdout...", file=sys.stderr)
for row in m2:
    print("\t".join((str(n) for n in row)))

