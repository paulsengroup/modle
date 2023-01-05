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
  python -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$1"
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

function compare_files {
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
  2>&1 echo "Usage: $0 path_to_modle_tools_bin"
  exit 1
fi

modle_tools_bin="$1"

data_dir="$(readlink -f "$(dirname "$0")/../data/integration_tests")"

matrix1="$data_dir/4DNFI9GMP2J8_chr20_25kbp_mt_eval.cool"
matrix2="$data_dir/4DNFIFJH2524_chr20_25kbp_mt_eval.cool"

ref=("$data_dir/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_"*.{bw,tsv.gz})


if ! check_files_exist "$matrix1" "$matrix2" "${ref[@]}"; then
  exit 1
fi

outdir="$(mktemp -d -t modle-tools-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

"$modle_tools_bin" eval -i "$matrix1" \
                        -r "$matrix2" \
                        -m custom     \
                        -w 3mbp       \
                        -o "$outdir/out_custom_score"

# mkdir -p /tmp/test/data/integration_tests/
# cp "$outdir/"*horizontal.bw" /tmp/test/data/integration_tests/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_horizontal.bw
# cp "$outdir/"*vertical.bw" /tmp/test/data/integration_tests/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_vertical.bw
# cp "$outdir/"*horizontal.tsv.gz" /tmp/test/data/integration_tests/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_horizontal.tsv.gz
# cp "$outdir/"*vertical.tsv.gz" /tmp/test/data/integration_tests/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_vertical.tsv.gz

tgt=("$outdir/out_custom_score_custom_metric_"*.{bw,tsv.gz})

if [ ${#ref[@]} -ne ${#tgt[@]} ]; then
  2>&1 echo "Expected ${#ref[@]} files, found ${#tgt[@]}!"
  2>&1 echo "Please file a bug at https://github.com/paulsengroup/modle"
  exit 1
fi

for i in "${!ref[@]}"; do
  if ! compare_files "${ref[$i]}" "${tgt[$i]}"; then
    status=1
  fi
done


if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
