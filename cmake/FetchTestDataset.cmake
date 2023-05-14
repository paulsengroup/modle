# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/7934763/files/modle_test_data.tar.xz?download=1
  EXPECTED_HASH SHA512=ff15fa17b1a418654827f46149d248224963c86e2b54ff2a289a2db6ccc71ff05464e02d479972df864698748ba09e6f334d716b6b3d91c8cf1cf08a1d785eaa
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.xz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.xz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
