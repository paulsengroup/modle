#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
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

function compare_files {
  set -o pipefail
  set -e

  if [[ "$1" == *.bw ]]; then
    compare_bwigs.py "$1" "$2"
    return "$?"
  fi

  2>&1 echo "Comparing $1 with $2..."
  if cmp <(xz -dc "$1") <(xz -dc "$2"); then
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

data_dir="$(readlink_py "$(dirname "$0")/../data/integration_tests")"

bwig="$data_dir/ENCSR942XQI_fc.bw"
bed="$data_dir/ENCSR942XQI_candidate_barriers.bed.xz"

ref="$data_dir/mt_annotate_barriers_reference_001.bed.xz"

if ! check_files_exist "$bwig" "$bed" "${ref[@]}"; then
  exit 1
fi

outdir="$(mktemp -d -t modle-tools-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT


outname="$outdir/out.bed.xz"
"$modle_tools_bin" annotate-barriers \
                   "$bwig" \
                   "$bed" |
                   xz -0 > "$outname"

if ! compare_files "$outname" "$ref"; then
    status=1
fi


if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
