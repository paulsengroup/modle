#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

status=0

if [ $# -ne 1 ]; then
  2>&1 echo "Usage: $0 path_to_modle_tools_bin"
  status=1
fi

if ! command -v cooler &> /dev/null; then
  2>&1 echo "Unable to find cooler in your PATH"
  status=1
fi

# Try to detect the error outlined below as early as possible:
# https://github.com/open2c/cooler/pull/298
cooler --help > /dev/null

if [ $status -ne 0 ]; then
  exit $status
fi

modle_tools_bin="$1"

data_dir="$(readlink -f "$(dirname "$0")/../data/integration_tests")"

matrix="$data_dir/4DNFI9GMP2J8_chr20_25kbp.cool"
matrix_blurred="$data_dir/4DNFI9GMP2J8_chr20_25kbp_blurred.cool"
matrix_dog="$data_dir/4DNFI9GMP2J8_chr20_25kbp_dog.cool"

if [ ! -f "$matrix_blurred" ]; then
  2>&1 echo "Unable to find test file \"$matrix_blurred\""
  status=1
fi

if [ ! -f "$matrix_dog" ]; then
  2>&1 echo "Unable to find test file \"$matrix_dog\""
  status=1
fi

if [ ! "$status" -eq 0 ]; then
  exit 1
fi

outdir="$(mktemp -d -t modle-tools-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

"$modle_tools_bin" transform -m gaussian_blur \
                             -w 3mbp          \
                             -i "$matrix"     \
                             -o "$outdir/out_blurred.cool"

"$modle_tools_bin" transform -m difference_of_gaussians \
                             -w 3mbp                    \
                             -i "$matrix"               \
                             -o "$outdir/out_dog.cool"

# mkdir -p /tmp/test/data/integration_tests/
# cp "$outdir/out_blurred.cool" /tmp/test/data/integration_tests/4DNFI9GMP2J8_chr20_25kbp_blurred.cool
# cp "$outdir/out_dog.cool" /tmp/test/data/integration_tests/4DNFI9GMP2J8_chr20_25kbp_dog.cool

function compare_coolers {
  set -o pipefail
  set -e

  2>&1 echo "Comparing $1 with $2..."
  if diff <(cooler dump -t chroms "$1") \
          <(cooler dump -t chroms "$2") \
     && \
     diff <(cooler dump --join "$1") \
          <(cooler dump --join "$2");
  then
    2>&1 echo "Files are identical"
    return 0
  else
    2>&1 echo "Files differ"
    return 1
  fi
}

status=0

if ! compare_coolers "$outdir/out_blurred.cool" "$data_dir/4DNFI9GMP2J8_chr20_25kbp_blurred.cool"; then
  status=1
fi

if ! compare_coolers "$outdir/out_dog.cool" "$data_dir/4DNFI9GMP2J8_chr20_25kbp_dog.cool"; then
  status=1
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
