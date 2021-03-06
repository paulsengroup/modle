#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

if [ $# -ne 5 ]; then
    echo "Usage:    $0 file.bed chr chr_start chr_end prob_of_barr_block"
    echo "Example:  $0 file.bed chr7 78000000 87000000 0.8"
    exit 1
fi

awk -v chr="$2"   \
    -v start="$3" \
    -v end="$4"   \
    -v prob_of_block="$5" \
'
BEGIN {
    # Set field separator to use for output
    OFS = "\t";
}
{
    if ($1 == chr && $2 > start && $3 < end) {
        # Make coordinates zero-based
        $2 -= start;
        $3 -= start;
        if (NF >= 8) {
          $7 -= start;
          $8 -= start;
        }
        $5 = prob_of_block;
        # In this particular file, strand information is stored in the name field.
        # The following copies the name column to the strand column
        $6 = tolower($4);
        #$4 = "";
        # Replace forward/reverse with + and - and anything else with .
        sub(/forward/, "+", $6);
        sub(/reverse/, "-", $6);
        sub(/[^+-]+/, ".", $6);
        gsub(/\t+/, "\t", $0);
        print $0;
    }
}
' \
"$1"
