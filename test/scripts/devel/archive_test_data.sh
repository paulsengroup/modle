#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


set -e
set -o pipefail
set -u


git_root="$(git rev-parse --show-toplevel)"

(cd "$git_root/test/data" &&
  find . -type f -exec sha256sum {} + |
      grep -vF checksums.sha256 |
      grep -vF README.md |
      grep -vF modle_test_data.tar.xz |
      sort -k2,2V |
      tee checksums.sha256 > /dev/null
)

tar -C "$git_root" -cf - test/data/{integration,unit}_tests/ |
  xz -9 --extreme -T0 |
  tee "$git_root/test/data/modle_test_data.tar.xz" > /dev/null

sha512sum "$git_root/test/data/modle_test_data.tar.xz"
