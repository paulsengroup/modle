# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/6638906/files/modle_test_data.tar.gz?download=1
  EXPECTED_HASH SHA512=65808e251056f125fa0f7aeb523281a173f057cce184931af5d73d81682de50f8553a9f7ee58d894008b95d0fbe084f744ee8ee4ebf8c9d34454efdf53c3d07e
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
