#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


set -e
set -o pipefail
set -u


git_root="$(git rev-parse --show-toplevel)"

git_root="$(git rev-parse --show-toplevel)"
PATH="$git_root/test/scripts/devel:$PATH"

dest_dir="$git_root/test/data/unit_tests/correlation/"
mkdir -p "$dest_dir"

correlation_scipy.py 250 1000     0 15000  uint | xz -9 --extreme -T0 > "$dest_dir/correlation_uint.txt.xz" &
correlation_scipy.py 250 1000 -7250  7250   int | xz -9 --extreme -T0 > "$dest_dir/correlation_int.txt.xz" &
correlation_scipy.py 250 1000 -7250  7250 float | xz -9 --extreme -T0 > "$dest_dir/correlation_float.txt.xz" &

correlation_scipy.py 5 100000     0 15000  uint | xz -9 --extreme -T0 > "$dest_dir/correlation_uint_long_vect.txt.xz" &
correlation_scipy.py 5 100000 -7250  7250   int | xz -9 --extreme -T0 > "$dest_dir/correlation_int_long_vect.txt.xz" &
correlation_scipy.py 5 100000 -7250  7250 float | xz -9 --extreme -T0 > "$dest_dir/correlation_float_long_vect.txt.xz" &

wait
