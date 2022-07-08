#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

ok=true

if [ $# -ne 1 ]; then
  2>&1 echo "Usage: $0 path_to_modle_bin"
  ok=false
fi

if ! command -v h5diff &> /dev/null; then
  2>&1 echo "Unable to find h5diff in your PATH"
  ok=false
fi

if ! command -v xz &> /dev/null; then
  2>&1 echo "Unable to find xz in your PATH"
  ok=false
fi

if ! $ok; then
  exit 1
fi

modle_bin="$1"

data_dir="$(readlink -f "$(dirname "$0")/../data/integration_tests")"

chrom_sizes="$data_dir/grch38.chrom.sizes"
extr_barriers="$data_dir/grch38_h1_extrusion_barriers.bed.xz"

status=0
if [ ! -f "$chrom_sizes" ]; then
  2>&1 echo "Unable to find test file \"$chrom_sizes\""
  status=1
fi

if [ ! -f "$extr_barriers" ]; then
  2>&1 echo "Unable to find test file \"$extr_barriers\""
  status=1
fi

if [ ! "$status" -eq 0 ]; then
  exit 1
fi

outdir="$(mktemp -d -t modle-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

# Only include chr2, chr20, chr21 and chr22
"$modle_bin" sim -c <(grep '^chr2' "$chrom_sizes") \
                 -b <(xz -dc "$extr_barriers" | grep '^chr2') \
                 -o "$outdir/out" \
                 -r 20kb \
                 --target-contact-density 20 \
                 --ncells 4 \
                 --max-burnin-epochs 5000

# cp "$outdir/out.cool" /tmp/test/data/integration_tests/reference_001.cool
echo "Comparing $outdir/out.cool with $data_dir/reference_001.cool..."

function compare_group () {
  # This is a workaround to make the test fail if two groups/datasets are not comparable
  msg="$(h5diff -c "$outdir/out.cool"             \
                   "$data_dir/reference_001.cool" \
                   "$1" 2>&1)"

  if [ -n "$msg" ]; then
    >&2 echo "$msg"
    return 1
  fi

  return 0
}

status=0

# We test each group individually because older versions of h5diff do not have
# the --exclude-attribute flag
if ! compare_group chroms; then status=1; fi
if ! compare_group bins; then status=1; fi
if ! compare_group pixels; then status=1; fi
if ! compare_group indexes; then status=1; fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
