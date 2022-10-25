# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/7249712/files/modle_test_data.tar.gz?download=1
  EXPECTED_HASH SHA512=f710cb1cfacf6fccbe2d2826e1f6c731734b457906870cbf172ab205b01e774ab406f3e8372262d4f7f2dccf6c2eddce37c90f695a23b18d40373f8fd61e7b72
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/modle_test_data.tar.gz
  DESTINATION
  ${CMAKE_CURRENT_SOURCE_DIR})
