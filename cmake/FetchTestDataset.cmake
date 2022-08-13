# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/6989635/files/modle_test_data.tar.gz?download=1
  EXPECTED_HASH SHA512=d59c4768173a7e10516b2e60db1248660ee863421b3d698d015dc84de881ab075f0003208a8a77a3224f98f5b58d083c1bfc474eec81f3cfac753c66b278d8ff
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
