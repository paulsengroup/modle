#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


set -e
set -o pipefail
set -u


git_root="$(git rev-parse --show-toplevel)"

img_url='https://upload.wikimedia.org/wikipedia/commons/thumb/2/29/Japanese_Squirrel_edited_version.jpg/640px-Japanese_Squirrel_edited_version.jpg'
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

(cd "$tmpdir" && curl -LO "$img_url")

img="$tmpdir/$(basename "$img_url")"

git_root="$(git rev-parse --show-toplevel)"
PATH="$git_root/test/scripts:$PATH"

dest_dir="$git_root/test/data/unit_tests/contact_matrices/"
mkdir -p "$dest_dir"


# Generate fp test matrix in cool and TSV format
2>&1 echo "Generating test matrix (FP)..."
generate_test_cooler.py --add-noise "$img" --img-source="$img_url" "$dest_dir/contact_matrix_dense_fp_001.cool"
cooler_to_dense_matrix.py "$dest_dir/contact_matrix_dense_fp_001.cool" test |
  xz -9 --extreme -T0 > "$dest_dir/contact_matrix_dense_fp_001.tsv.xz" &

# Generate int test matrix in cool and TSV format
2>&1 echo "Generating test matrix (INT)..." && \
generate_test_cooler.py "$img" --img-source="$img_url" "$dest_dir/contact_matrix_dense_int_001.cool" && \
cooler_to_dense_matrix.py "$dest_dir/contact_matrix_dense_int_001.cool" test |
  xz -9 --extreme -T0 > "$dest_dir/contact_matrix_dense_int_001.tsv.xz"


# Generate blurred matrices
mkdir -p "$dest_dir/blurred/"
for sigma in 0.01 0.50 1.00 1.50 5.00; do
    src="$dest_dir/contact_matrix_dense_int_001.tsv.xz"
    dest="$dest_dir/blurred/contact_matrix_dense_int_001_blurred_${sigma}.tsv.xz"

    2>&1 echo "Generating blurred matrix with sigma=$sigma..." && \
    xz -dc "$src" | gaussian_blur_scipy.py "$sigma" | xz -9 --extreme -T0 > "$dest" &
done

wait
