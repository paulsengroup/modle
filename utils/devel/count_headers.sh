#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

wd="$(dirname "$0")"

if [ $# -gt 2 ]; then
    echo >&2 "Usage:   $0 path_relative_to_repo_root/"
    echo >&2 "Default: $0 src/"
fi


path='src/'
if [ $# -eq 2 ]; then
    path="$1"
fi

cat $(find "$wd/../$1" -type f -name "*.hpp" | tr '\n' ' ') | grep -oE '#include .[[:alnum:]_/\.\-]+.' | sort | uniq -c | sort -r
