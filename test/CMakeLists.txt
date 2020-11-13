add_executable(modle_gtest)
target_link_libraries(modle_gtest)

target_sources(
  modle_gtest
  PRIVATE units/contacts_test.cpp
  PRIVATE units/bed_parser_test.cpp
  PRIVATE units/correlation_test.cpp
  PRIVATE units/utils_test.cpp)

target_include_directories(modle_gtest SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/googletest/include/)

target_link_libraries(
  modle_gtest
  PRIVATE project_options
  PUBLIC libmodle CONAN_PKG::gtest gtest_main)
