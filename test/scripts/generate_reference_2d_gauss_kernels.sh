#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


set -e
set -o pipefail
set -u

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

wd="$tmpdir/reference_gaussian_kernels"
mkdir "$wd"

git_root="$(git rev-parse --show-toplevel)"
PATH="$git_root/test/scripts:$PATH"
dest="$git_root/test/data/unit_tests/reference_gaussian_kernels.tar.xz"

for size in $(seq 1 15); do
  for sigma in 0.5 1.0 1.5 2.5 6.3 10.0; do
    generate_2d_gauss_kernel.py "$size" "$sigma" |
      tee "$wd/gaussian_kernel_${size}_${sigma}.csv" > /dev/null
  done
done

mkdir -p "$(dirname "$dest")"

tar -C "$tmpdir" -cf - "$(basename "$wd")" |
  xz -9 --extreme | tee "$dest" > /dev/null

2>&1 echo "$dest"
