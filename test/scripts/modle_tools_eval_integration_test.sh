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

if ! command -v shasum &> /dev/null; then
  2>&1 echo "Unable to find shasum in your PATH"
  ok=false
fi

if ! $ok; then
  exit 1
fi

modle_tools_bin="$1"

data_dir="$(readlink -f "$(dirname "$0")/../data/integration_tests")"

matrix1="$data_dir/4DNFI9GMP2J8_chr20_25kbp_mt_eval.cool"
matrix2="$data_dir/4DNFIFJH2524_chr20_25kbp_mt_eval.cool"

status=0
if [ ! -f "$matrix1" ]; then
  2>&1 echo "Unable to find test file \"$matrix1\""
  status=1
fi

if [ ! -f "$matrix2" ]; then
  2>&1 echo "Unable to find test file \"$matrix2\""
  status=1
fi

if [ ! "$status" -eq 0 ]; then
  exit 1
fi

outdir="$(mktemp -d -t modle-tools-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

"$modle_tools_bin" eval -i "$matrix1" \
                        -r "$matrix2" \
                        -m custom     \
                        -w 3mbp       \
                        -o "$outdir/out_custom_score"

# cp "$outdir/"*horizontal.bw" /tmp/test/data/integration_tests/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_horizontal.bw
# cp "$outdir/"*vertical.bw" /tmp/test/data/integration_tests/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_vertical.bw
# cp "$outdir/"*horizontal.tsv.gz" /tmp/test/data/integration_tests/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_horizontal.tsv.gz
# cp "$outdir/"*vertical.tsv.gz" /tmp/test/data/integration_tests/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_vertical.tsv.gz

(cd "$data_dir" && shasum -c checksums.sha256)

status="$?"

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
  exit "$status"
fi
