# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED)
find_package(Catch2 CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(HDF5 CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)

include(CTest)
include(Catch)

add_executable(test_main)

target_sources(
  test_main
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/units/common/cli_utils_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/common/const_map_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/common/dna_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/contact_matrix/contact_matrix_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/interval_tree/interval_tree_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/libmodle_io/bed_parser_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/libmodle_io/bed_tree_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/libmodle_io/compressed_io_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/libmodle_io/cooler_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/libmodle_io/hdf5_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/simulation_cpu/collision_encoding_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/simulation_cpu/common.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/simulation_cpu/simulation_complex_unit_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/simulation_cpu/simulation_simple_unit_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/simulation_internal/extrusion_barriers_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/stats/correlation_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/stats/correlation_utils_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/stats/descriptive_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/stats/misc_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/stats/tests_test.cpp)

target_include_directories(test_main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/units/include)

target_link_libraries(
  test_main
  PRIVATE project_warnings
          project_options
          Modle::cmatrix
          Modle::common
          Modle::interval_tree
          Modle::io_bed
          Modle::io_bigwig
          Modle::io_chrom_sizes
          Modle::io_compressed
          Modle::io_cooler
          Modle::io_hdf5
          Modle::libmodle_cpu
          Modle::libmodle_internal
          Modle::stats)

target_link_system_libraries(
  test_main
  PRIVATE
  absl::flat_hash_set
  absl::hash
  absl::span
  absl::strings
  Boost::filesystem
  Boost::headers
  Catch2::Catch2WithMain
  fmt::fmt
  HDF5::HDF5
  spdlog::spdlog
  thread-pool)

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
  OUTPUT_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/Testing/
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
  OUTPUT_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/Testing/
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
  OUTPUT_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/Testing/
  EXTRA_ARGS
  --success)
