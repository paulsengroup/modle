#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

status=0

if [ $# -ne 1 ]; then
  2>&1 echo "Usage: $0 path_to_modle_bin"
  status=1
fi

if ! command -v cooler &> /dev/null; then
  2>&1 echo "Unable to find cooler in your PATH"
  status=1
fi

if ! command -v shasum &> /dev/null; then
  2>&1 echo "Unable to find shasum in your PATH"
  status=1
fi

if ! command -v xz &> /dev/null; then
  2>&1 echo "Unable to find xz in your PATH"
  status=1
fi

if [ $status -ne 0 ]; then
  exit $status
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
                 --track-1d-lef-position \
                 --max-burnin-epochs 5000

# cp "$outdir/out_lef_1d_occupancy.bw" /tmp/test/data/integration_tests/reference_001.bw
# cp "$outdir/out.cool" /tmp/test/data/integration_tests/reference_001.cool

function compare_coolers {
  set -o pipefail
  2>&1 echo "Comparing $1 with $2..."
  diff --report-identical-files   \
       --brief                    \
       <(cooler dump --join "$1") \
       <(cooler dump --join "$2")
}

if ! compare_coolers "$outdir/out.cool" "$data_dir/reference_001.cool"; then
  status=1
fi

bw_checksum="54eb48e520176d0b80a9ee6df66a9d239f5d08e0237499320f4d6b8d1d70b972"
if ! shasum -c <(echo "$bw_checksum  $outdir/out_lef_1d_occupancy.bw"); then
  status=1
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
