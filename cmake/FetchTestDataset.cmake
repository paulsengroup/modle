# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/7455247/files/modle_test_data.tar.gz?download=1
  EXPECTED_HASH SHA512=a3d22f1981e9c96012188e559ad8684919fcdb18b3405bd581b9f8c80239b1566674de8eb9cb992e809449e619c4e9e9e229fca36c14e4f5671b896fb6ee3dfb
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
