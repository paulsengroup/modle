# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/7499258/files/modle_test_data.tar.xz?download=1
  EXPECTED_HASH SHA512=794c7f8dd0cf0725bd6ec3a44b1da494b76cb899dc48816ab5668c4d43390a8e77b71adb1fd27432ca12486f4056bdc4376f21e8111b13cf3fd5ae56eef2f86d
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.xz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.xz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
