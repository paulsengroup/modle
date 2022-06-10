# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/6631983/files/modle_test_data.tar.gz?download=1
  EXPECTED_HASH SHA512=135712297e460088333162199afdd44a66ba82bde42b093f842fc776f883aaf31137b4f3fa23218eb1695f759d5f074c43afc5e39d4d2de7404d3ef14148db21
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
