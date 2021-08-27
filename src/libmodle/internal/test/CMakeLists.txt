find_package(Catch2 REQUIRED)

add_library(libmodle_internal_test INTERFACE)
add_library(Modle::libmodle_internal_test ALIAS libmodle_internal_test)

target_include_directories(libmodle_internal_test INTERFACE include/)

target_sources(libmodle_internal_test INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_barriers_test.hpp)

target_link_libraries(
  libmodle_internal_test
  INTERFACE project_warnings
            project_options
            Catch2::Catch2
            Modle::common
            Modle::libmodle_internal)
