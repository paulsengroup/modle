# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/6991251/files/modle_test_data.tar.gz?download=1
  EXPECTED_HASH SHA512=de867c8750325c3ed713c7de92978abfe6b467bbd0f7329c3ca4ddc90ddb4d29ccf1d7981aef6426bd50379e2125b51ec7a2d04562a5e3b32fb87fe4f49f0771
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
