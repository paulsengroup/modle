find_package(Boost REQUIRED)
find_package(Catch2 REQUIRED)

include(CTest)
include(Catch)

add_executable(catch_main)

target_sources(
  catch_main
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/catch_main.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/contact_matrix/contacts_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/correlation/common.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/correlation/correlation_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/correlation/correlation_utils_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/interval_tree/interval_tree_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/libmodle_io/bed_parser_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/libmodle_io/bed_tree_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/libmodle_io/compressed_io_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/libmodle_io/cooler_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/libmodle_io/hdf5_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_cpu/common.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_cpu/simulation_complex_unit_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_cpu/simulation_simple_unit_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_internal/extrusion_barriers_test.cpp)

target_link_libraries(
  catch_main
  PRIVATE project_warnings
          project_options
          Boost::headers
          Catch2::Catch2
          Modle::cmatrix
          Modle::correlation
          Modle::interval_tree
          Modle::libmodle
          # Modle::libmodle_internal
          Modle::libmodle_io)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Testing/)

# automatically discover tests that are defined in catch based test files you can modify the unittests. TEST_PREFIX to
# whatever you want, or use different for different binaries
catch_discover_tests(
  catch_main
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
  catch_main
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
  catch_main
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
