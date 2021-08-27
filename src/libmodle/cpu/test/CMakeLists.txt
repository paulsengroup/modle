find_package(absl REQUIRED)
find_package(Boost REQUIRED)
find_package(Catch2 REQUIRED)

add_library(libmodle_cpu_test INTERFACE)
add_library(Modle::libmodle_cpu_test ALIAS libmodle_cpu_test)

target_include_directories(libmodle_cpu_test INTERFACE include/)

target_sources(libmodle_cpu_test INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/simulation_complex_unit_tests.hpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/simulation_simple_unit_tests.hpp)

target_link_libraries(
  libmodle_cpu_test
  INTERFACE project_warnings
            project_options
            absl::span
            absl::strings
            # Boost
            Catch2::Catch2
            Modle::libmodle_io
            Modle::libmodle)
