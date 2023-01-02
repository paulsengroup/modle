#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


set -e
set -o pipefail
set -u

git_root="$(git rev-parse --show-toplevel)"
PATH="$git_root/test/scripts:$PATH"
dest_dir="$git_root/test/data/unit_tests/reference_gaussian_kernels/"

mkdir -p "$dest_dir"

for size in $(seq 1 15); do
  for sigma in 0.5 1.0 1.5 2.5 6.3 10.0; do
    generate_2d_gauss_kernel.py "$size" "$sigma" |
      xz -9 --extreme |
      tee "$dest_dir/gaussian_kernel_${size}_${sigma}.tsv.xz" > /dev/null &
  done
done

wait
