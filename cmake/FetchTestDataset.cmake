# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/7506960/files/modle_test_data.tar.xz?download=1
  EXPECTED_HASH SHA512=c7b6b5fa8e016ce8ed69412439601751d7b6315426ff80f77d15786f583e6baaacf39432616d228bdbd888646e9ffa3a3e9c0cb02e804bc37d2aa350c40f1100
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.xz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.xz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
