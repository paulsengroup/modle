# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/6638790/files/modle_test_data.tar.gz?download=1
  EXPECTED_HASH SHA512=2f82e008429e4c79074e1acecbb9497af15b826e0ac604e2f0cb72c016bc9dee7aaea1c29330a599c0035f85048119590ffd347a2ccff603c413eb9edb23d857
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
