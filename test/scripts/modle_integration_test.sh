#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

# readlink -f is not available on macos...
function readlink_py {
  set -eu
  python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$1"
}

function check_files_exist {
  set -eu
  status=0
  for f in "$@"; do
    if [ ! -f "$f" ]; then
      2>&1 echo "Unable to find test file \"$f\""
      status=1
    fi
  done

  return "$status"
}

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

function compare_bwigs {
  set -o pipefail
  set -e

  2>&1 echo "Comparing $1 with $2..."
  if cmp "$1" "$2"; then
    2>&1 echo "Files are identical"
    return 0
  else
    2>&1 echo "Files differ"
    return 1
  fi
}

export function readlink_py

status=0

if [ $# -ne 1 ]; then
  2>&1 echo "Usage: $0 path_to_modle_bin"
  status=1
fi

modle_bin="$1"

data_dir="$(readlink_py "$(dirname "$0")/../data/integration_tests")"
script_dir="$(readlink_py "$(dirname "$0")")"

chrom_sizes="$data_dir/grch38.chrom.sizes"
regions_bed="$data_dir/grch38_regions_of_interest.bed"
extr_barriers="$data_dir/grch38_h1_extrusion_barriers.bed.xz"
ref_cooler="$data_dir/reference_001.cool"
ref_bw="$data_dir/reference_001.bw"

export PATH="$PATH:$script_dir"

if ! command -v cooler &> /dev/null; then
  2>&1 echo "Unable to find cooler in your PATH"
  status=1
fi

# Try to detect the error outlined below as early as possible:
# https://github.com/open2c/cooler/pull/298
cooler --help > /dev/null

if ! command -v compare_bwigs.py &> /dev/null; then
  2>&1 echo "Unable to find compare_bwigs.py"
  status=1
fi

if ! command -v xz &> /dev/null; then
  2>&1 echo "Unable to find xz in your PATH"
  status=1
fi

if [ $status -ne 0 ]; then
  exit $status
fi

if ! check_files_exist "$chrom_sizes" "$extr_barriers" "$ref_cooler" "$ref_bw" "$regions_bed"; then
  exit 1
fi

outdir="$(mktemp -d -t modle-ci-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

"$modle_bin" sim -c "$chrom_sizes" \
                 -g "$regions_bed" \
                 -b "$extr_barriers" \
                 -o "$outdir/out" \
                 -r 20kb \
                 --verbose \
                 --target-contact-density 20 \
                 --ncells 2 \
                 --track-1d-lef-position \
                 --max-burnin-epochs 5000

# mkdir -p /tmp/test/data/integration_tests/
# cp "$outdir/out_lef_1d_occupancy.bw" /tmp/test/data/integration_tests/reference_001.bw
# cp "$outdir/out.cool" /tmp/test/data/integration_tests/reference_001.cool

if ! compare_coolers "$outdir/out.cool" "$ref_cooler"; then
  status=1
fi

if ! compare_bwigs.py "$ref_bw" "$outdir/out_lef_1d_occupancy.bw"; then
  status=1
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
