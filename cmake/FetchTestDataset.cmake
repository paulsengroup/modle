# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://www.dropbox.com/s/avqohi7xzj0voi2/modle_test_data.tar.xz?dl=1
  EXPECTED_HASH 65f287835560917382097acb22ace9f086338d988cce0827188acc846e18b7de8b8c380a4145e51d81ebd2bf5d0282af12dd7b0dd955f32bebcdfd2348d5f148
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.xz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.xz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
