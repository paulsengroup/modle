#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

ok=true

if [ $# -ne 1 ]; then
  2>&1 echo "Usage: $0 path_to_modle_tools_bin"
  ok=false
fi

if ! command -v h5diff &> /dev/null; then
  2>&1 echo "Unable to find h5diff in your PATH"
  ok=false
fi

if ! $ok; then
  exit 1
fi

modle_tools_bin="$1"

data_dir="$(readlink -f "$(dirname "$0")/../data/integration_tests")"

matrix="$data_dir/4DNFI9GMP2J8_chr20_25kbp.cool"
matrix_blurred="$data_dir/4DNFI9GMP2J8_chr20_25kbp_blurred.cool"
matrix_dog="$data_dir/4DNFI9GMP2J8_chr20_25kbp_dog.cool"

status=0
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

# cp "$outdir/out_blurred.cool" /tmp/test/data/integration_tests/4DNFI9GMP2J8_chr20_25kbp_blurred.cool
# cp "$outdir/out_dog.cool" /tmp/test/data/integration_tests/4DNFI9GMP2J8_chr20_25kbp_dog.cool

function compare_group () {
  ref="$1"
  tgt="$2"
  grp="$3"
  # This is a workaround to make the test fail if two groups/datasets are not comparable
  msg="$(h5diff -c "$ref" "$tgt" "$grp" 2>&1)"

  if [ -n "$msg" ]; then
    >&2 echo "$msg"
    return 1
  fi

  return 0
}

status=0

echo "Comparing $outdir/out_blurred.cool with $data_dir/4DNFI9GMP2J8_chr20_25kbp_blurred.cool..."
# We test each group individually because older versions of h5diff do not have
# the --exclude-attribute flag
if ! compare_group "$outdir/out_blurred.cool" "$data_dir/4DNFI9GMP2J8_chr20_25kbp_blurred.cool" chroms; then status=1; fi
if ! compare_group "$outdir/out_blurred.cool" "$data_dir/4DNFI9GMP2J8_chr20_25kbp_blurred.cool" bins; then status=1; fi
if ! compare_group "$outdir/out_blurred.cool" "$data_dir/4DNFI9GMP2J8_chr20_25kbp_blurred.cool" pixels; then status=1; fi
if ! compare_group "$outdir/out_blurred.cool" "$data_dir/4DNFI9GMP2J8_chr20_25kbp_blurred.cool" indexes; then status=1; fi


echo "Comparing $outdir/out_dog.cool with $data_dir/4DNFI9GMP2J8_chr20_25kbp_dog.cool..."
if ! compare_group "$outdir/out_dog.cool" "$data_dir/4DNFI9GMP2J8_chr20_25kbp_dog.cool" chroms; then status=1; fi
if ! compare_group "$outdir/out_dog.cool" "$data_dir/4DNFI9GMP2J8_chr20_25kbp_dog.cool" bins; then status=1; fi
if ! compare_group "$outdir/out_dog.cool" "$data_dir/4DNFI9GMP2J8_chr20_25kbp_dog.cool" pixels; then status=1; fi
if ! compare_group "$outdir/out_dog.cool" "$data_dir/4DNFI9GMP2J8_chr20_25kbp_dog.cool" indexes; then status=1; fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
  exit "$status"
fi
