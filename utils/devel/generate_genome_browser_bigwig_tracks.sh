#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e

for f in $@; do
    url="$(dropbox sharelink "${f/.*Dropbox.*\///}")"
    echo "track type=bigWig name=$(basename "$f") bigDataUrl=${url/dl=0$/dl=1}"
done;
