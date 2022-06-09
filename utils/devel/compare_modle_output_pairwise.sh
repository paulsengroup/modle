#!/usr/bin/env bash

# Copyright (c) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -x
set -o pipefail

# This helper script is used to visually compare the output of two versions of MoDLE.
# This can be useful to assess whether a change produces (unwanted) systematic differences in the output.

# Script dependencies:
# - Cooler
# - HiCExplorer
# - ImageMagick

argc="$#"

if [ $argc -lt 5 ]; then
  echo "Usage: $0 modle_bin_v1 modle_bin_v2 output_dir output_basename region_of_interest modle_param1 modle_param2 ..."
  echo "Example: $0 v1/modle v2/modle /tmp/modle_comparison my_comparison chr10:68000000-70000000 -c GRCh38.chrom.sizes --ncells 64 --extrusion-barrier-file barriers.bed --force"
  exit 1
fi

if ! command -v cooler --help &> /dev/null; then
  echo "cooler could not be found!"
  exit 1
fi

if ! command -v hicCompareMatrices --help &> /dev/null; then
  echo "hicCompareMatrices (HiCExplorer) suite could not be found!"
  exit 1
fi

if ! command -v montage &> /dev/null; then
  echo "montage command (ImageMagick) could not be found!"
  exit 1
fi

if ! command -v convert &> /dev/null; then
  echo "convert command (ImageMagick) could not be found!"
  exit 1
fi

modle_v1="$1"
modle_v2="$2"
outdir="$3"
bname="$4"
region_of_interest="$5"

# shellcheck disable=SC2206
modle_params=(${@:6})
cpus=$(($(nproc) / 2))
plot_suffix="$(printf '%s' "$region_of_interest" | tr -cs '[:alnum:]' _)"

mkdir -p "$outdir/v1" "$outdir/v2" "$outdir/diff" "$outdir/collage"

"$modle_v1" sim "${modle_params[@]}" -o "$outdir/v1/$bname" &> /dev/null \
&& cooler zoomify -p $cpus "$outdir/v1/$bname.cool" &> /dev/null \
&& cooler show --dpi 600 -o "$outdir/v1/${bname}_$plot_suffix.png" \
    "$outdir/v1/$bname.cool" "$region_of_interest" \
     &> /dev/null &

"$modle_v2" sim "${modle_params[@]}" -o "$outdir/v2/$bname" &> /dev/null \
&& cooler zoomify -p "$cpus" "$outdir/v2/$bname.cool" &> /dev/null  \
&& cooler show --dpi 600 -o "$outdir/v2/${bname}_$plot_suffix.png" \
    "$outdir/v2/$bname.cool" "$region_of_interest" \
    &> /dev/null &

wait

hicCompareMatrices -m "$outdir/v1/$bname.cool" "$outdir/v2/$bname.cool" \
                   --outFileName "$outdir/diff/$bname.cool" &> /dev/null

cooler show --dpi 600 -o "$outdir/diff/${bname}_$plot_suffix.png" \
     --cmap bwr --scale linear --zmin='-10' --zmax=10 \
    "$outdir/diff/$bname.cool" "$region_of_interest" \
     &> /dev/null

rm "$outdir"/v?/"$bname.cool"

wd="$(mktemp -d)"

function cleanup {
    rm -r "$wd"
}

trap cleanup EXIT

function trim_and_add_caption {
    input="$1"
    output="$2"
    title="$3"

    tmp1="$wd/$(mktemp XXXXXXXXXX).png"

    convert "$input" -trim "$tmp1"
    montage -label "$title" "$tmp1" -pointsize 150 -geometry +0+0 -gravity South "$output"
}


trim_and_add_caption "$outdir/v1/${bname}_$plot_suffix.png" "$wd/fig1.png" "v1" &
trim_and_add_caption "$outdir/v2/${bname}_$plot_suffix.png" "$wd/fig2.png" "v2" &
trim_and_add_caption "$outdir/diff/${bname}_$plot_suffix.png" "$wd/fig3.png" "DIFF (v1 - v2)" &
wait

out="$outdir/collage/${bname}_$plot_suffix.png"

montage "$wd"/fig1.png "$wd"/fig2.png "$wd"/fig3.png -title "$(basename "$out" .png)" -pointsize 200 -tile 3x1 -geometry +0+0 "$wd/out.png"
convert "$wd/out.png" -scale '25%' "$out"
