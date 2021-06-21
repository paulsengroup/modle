add_library(libmodle_internal STATIC)
add_library(Modle::modle_internal ALIAS libmodle_internal)

target_sources(
  libmodle_internal
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/genome.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/genome_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_barriers.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_factors.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_factors_impl.hpp)
target_link_libraries(
  libmodle_internal
  PRIVATE project_options
          project_warnings
          Modle::common
          Modle::cmatrix
  PUBLIC Modle::io CONAN_PKG::abseil CONAN_PKG::xxhash)

target_include_directories(libmodle_internal PUBLIC include/)
set_target_properties(libmodle_internal PROPERTIES OUTPUT_NAME "modle_internal")
