# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

if(NOT WIN32)
  file(LOCK "${PROJECT_SOURCE_DIR}/test/data/" DIRECTORY GUARD FILE)
endif()

set(TEST_DATASET_TAR "${PROJECT_SOURCE_DIR}/test/data/modle_test_data.tar.zst")

message(STATUS "Fetching the test dataset")

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/records/10790566/files/modle_test_data.tar.zst?download=1
  EXPECTED_HASH SHA256=9db43cc25c0ecd8f6928a93f42f31814cd70ab45b7edcc3e80ca28955f174237
  "${TEST_DATASET_TAR}"
)
# cmake-format: on

message(STATUS "Fetching the test dataset - done")

message(STATUS "Extracting the test dataset")

file(ARCHIVE_EXTRACT INPUT "${TEST_DATASET_TAR}" DESTINATION "${PROJECT_SOURCE_DIR}")

message(STATUS "Extracting the test dataset - done")
message(STATUS "Test datasets can be found under \"${PROJECT_SOURCE_DIR}/test/data/\"")
