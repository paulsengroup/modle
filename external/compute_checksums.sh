#!/usr/bin/env sh

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

shasum -a512 *.xz | tee checksums.sha512
