# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/6628728/files/modle_test_data.tar.gz?download=1
  EXPECTED_HASH SHA512=687e915df37d60e8cfae02d1c041f5fe33ecd3262a159ca243697705f4890aecfd6588d0ffb292d86dae79b9d9a37c0ac20bbed7b9beade4223ffd569474de87
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
