find_package(absl REQUIRED)
find_package(Boost REQUIRED)
find_package(Catch2 REQUIRED)
find_package(fmt REQUIRED)

add_library(modle_correlation_test INTERFACE)
add_library(Modle::correlation_test ALIAS modle_correlation_test)

target_include_directories(modle_correlation_test INTERFACE include/)

target_sources(
  modle_correlation_test
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/common.hpp ${CMAKE_CURRENT_SOURCE_DIR}/correlation_test.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/correlation_utils_test.hpp)

target_link_libraries(
  modle_correlation_test
  INTERFACE project_warnings
            project_options
            absl::span
            Catch2::Catch2
            fmt::fmt
            Modle::common
            Modle::correlation)
