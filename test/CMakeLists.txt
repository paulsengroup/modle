find_package(Catch2 REQUIRED)

include(CTest)
include(Catch)

add_executable(catch_main catch_main.cpp)

target_link_libraries(
  catch_main
  PUBLIC Catch2::Catch2
  PRIVATE project_warnings
          project_options
          Modle::cmatrix_test
          Modle::correlation_test
          Modle::interval_tree_test
          Modle::libmodle_cpu_test
          Modle::libmodle_internal_test
          Modle::libmodle_io_test)

if(ENABLE_PCH)
  target_precompile_headers(catch_main INTERFACE <catch2/catch.hpp>)
endif()

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
