# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/7250085/files/modle_test_data.tar.gz?download=1
  EXPECTED_HASH SHA512=b30d8d7bc63fa3f4a8be418a29f16ed9ba55ca116d26dffc8ab6d0ced7d9e82b9e3bb6bcf96355f1cba74b870f25a38fb344f366cd84fc15abe5886d5f7a2b82
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
