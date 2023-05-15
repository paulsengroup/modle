# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/7937395/files/modle_test_data.tar.xz?download=1
  EXPECTED_HASH SHA512=7c74479d0570b2d59419f360054255efd1656e1401f442e8dcf2754d84bcb797011374ef3c18b2aa23f06aeb3c9fc64888e1c43a1bccaee6f5cf460c1b9d6e35
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.xz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.xz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
