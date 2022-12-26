# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/7254741/files/modle_test_data.tar.gz?download=1
  EXPECTED_HASH SHA512=5e24c4df6def21ba626f564efef161d8408fa864e811c2b904edb573e6001ca7535a8f715bfb0931fb01380767b2cb7ba301f7000bb65eeda53a054126256b79
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
