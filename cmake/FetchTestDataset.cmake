# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/7916701/files/modle_test_data.tar.xz?download=1
  EXPECTED_HASH SHA512=06a034a5a07b0ecef41c719ca56413c1175bfadc30029250539cc16600db3591af0918b88553bb28e0801e2a3e10fcf7c3d1c6e266c5f1dd4e4052bdba68846b
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.xz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.xz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
