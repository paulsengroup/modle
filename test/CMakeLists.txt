# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Boost REQUIRED QUIET)
find_package(Catch2 REQUIRED QUIET)
find_package(fmt REQUIRED QUIET)

include(CTest)
include(Catch)

add_executable(test_main)

target_sources(
  test_main
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_main.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/contact_matrix/contact_matrix_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/stats/common.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/stats/correlation_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/stats/correlation_utils_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/stats/descriptive_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/stats/tests_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/interval_tree/interval_tree_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/libmodle_io/bed_parser_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/libmodle_io/bed_tree_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/libmodle_io/compressed_io_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/libmodle_io/cooler_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/libmodle_io/hdf5_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/simulation_cpu/common.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/simulation_cpu/simulation_complex_unit_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/simulation_cpu/simulation_simple_unit_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/simulation_internal/extrusion_barriers_test.cpp)

target_link_libraries(
  test_main
  PRIVATE project_warnings
          project_options
          Boost::headers
          Catch2::Catch2
          fmt::fmt
          Modle::cmatrix
          Modle::common
          Modle::interval_tree
          Modle::libmodle
          Modle::stats)
# cmake-format: off
          # Not explicitly linking to the following targets saves quite a bit of time when compiling
          # Modle::libmodle_internal
          # Modle::libmodle_io)
# cmake-format: on

file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Testing/)

# automatically discover tests that are defined in catch based test files you can modify the unittests. TEST_PREFIX to
# whatever you want, or use different for different binaries
catch_discover_tests(
  test_main
  TEST_SPEC
  "[short]"
  TEST_SUFFIX
  " - SHORT"
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/../
  REPORTER
  xml
  OUTPUT_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/Testing/
  OUTPUT_SUFFIX
  .xml
  EXTRA_ARGS
  --success)

catch_discover_tests(
  test_main
  TEST_SPEC
  "[medium]"
  TEST_SUFFIX
  " - MEDIUM"
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/../
  REPORTER
  xml
  OUTPUT_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/Testing/
  OUTPUT_SUFFIX
  .xml
  EXTRA_ARGS
  --success)

catch_discover_tests(
  test_main
  TEST_SPEC
  "[long]"
  TEST_SUFFIX
  " - LONG"
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/../
  REPORTER
  compact # This dramatically speeds up the test suit
  EXTRA_ARGS
  --success)
