#!/usr/bin/env bash

set -e
set -o pipefail
set -u

ok=true

if [ $# -ne 1 ]; then
  2>&1 echo "Usage: $0 path_to_modle_bin"
  ok=false
fi

if ! command -v h5diff; then
  2>&1 echo "Unable to find h5diff in your PATH"
  ok=false
fi

if ! command -v xz 2; then
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

outdir="$(mktemp -d -t ci-XXXXXXXXXX)"

trap 'rm -rf -- "$outdir"' EXIT

# Only include chr2, chr20, chr21 and chr22
"$modle_bin" sim -c <(grep '^chr2' "$chrom_sizes") \
                 -b <(xz -dc "$extr_barriers" | grep '^chr2') \
                 -o "$outdir/out" \
                 -r 20kb \
                 --rev-extrusion-speed 2.5kb \
                 --fwd-extrusion-speed 2.5kb \
                 --ncells 4 \
                 --max-burnin-epochs 5000

echo "Comparing $outdir/out.cool with $data_dir/reference_001.cool..."
h5diff --exclude-attribute / \
       -c "$outdir/out.cool" \
          "$data_dir/reference_001.cool"
