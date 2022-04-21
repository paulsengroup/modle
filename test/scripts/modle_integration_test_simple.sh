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

chrom_sizes="$data_dir/grch38.chrom.sizes"
extr_barriers="$data_dir/grch38_h1_extrusion_barriers.bed.xz"

outdir="$(mktemp -d -t ci-XXXXXXXXXX)"

trap 'rm -rf -- "$outdir"' EXIT

# Only include chr2, chr20, chr21 and chr22
"$modle_bin" sim -c <(grep '^chr2' "$chrom_sizes") \
                 -b <(xz -dc "$extr_barriers" | grep '^chr2') \
                 -o "$outdir/out" \
                 -r 20000 \
                 --ncells 4 \
                 --max-burnin-epochs 5000

echo "Comparing $outdir/out.cool with $data_dir/reference_001.cool..."
h5diff --exclude-attribute / \
       -c "$outdir/out.cool" \
          "$data_dir/reference_001.cool"
