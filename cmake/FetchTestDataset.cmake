# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/7499993/files/modle_test_data.tar.xz?download=1
  EXPECTED_HASH SHA512=90d7899f9f86d86c87bd79de1a88b237e8de6265510ecdb8ca046204c26fd7752d76066b1bb17b55c3b03c876258a2ac25fa533a12f280ddbafbbb2c8a15c2bf
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.xz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.xz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
