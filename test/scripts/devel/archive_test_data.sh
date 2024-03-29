#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


set -e
set -o pipefail
set -u


git_root="$(git rev-parse --show-toplevel)"

(cd "$git_root/test/data" &&
  find . -type f \
         \( -path './integration_tests/*' \
            -or \
            -path './unit_tests/*' \) \
         -exec sha256sum {} + |
      sort -k2,2V |
      tee checksums.sha256 > /dev/null
)

tar -C "$git_root" -cf - test/data/{integration,unit}_tests/ test/data/checksums.sha256 |
  xz -9 --extreme -T0 |
  tee "$git_root/test/data/modle_test_data.tar.xz" > /dev/null

tar -tf "$git_root/test/data/modle_test_data.tar.xz"
sha512sum "$git_root/test/data/modle_test_data.tar.xz"
