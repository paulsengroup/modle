#!/usr/bin/env bash

set -e
set -o pipefail
set -u

if [ $# -ne 1 ]; then
  2>&1 echo "Usage: $0 path_to_modle_bin"
  exit 1
fi

modle_bin="$1"

data_dir="$(readlink -f "$(dirname "$0")/../data/integration_tests")"

chrom_sizes="$data_dir/GRCh38.chrom.sizes"
extr_barriers="$data_dir/GRCh38_H1_RAD21_occupancy.bed.xz"

outdir="$(mktemp -d -t ci-XXXXXXXXXX)"

trap 'rm -rf -- "$outdir"' EXIT

"$modle_bin" sim -c <(grep '^chr22' "$chrom_sizes") \
                 --extrusion-barrier-file <(xz -dc "$extr_barriers" | grep '^chr22') \
                 -o "$outdir/out" \
                 --bin-size 15000 \
                 --ncells "$(nproc)" \
                 --max-burnin-epochs 1000
