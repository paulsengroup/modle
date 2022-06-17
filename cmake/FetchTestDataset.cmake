# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/6656879/files/modle_test_data.tar.gz?download=1
  EXPECTED_HASH SHA512=358a33408a8e307c96a24d32c8a30ad21ee2e5871c02b9f27db43c9b0cd04056761c47847327834785d453db92c18c71af2ed11a8b887e8c4bd7244129987c89
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
