add_library(libmodle_internal)
add_library(Modle::modle_internal ALIAS libmodle_internal)

target_sources(
  libmodle_internal
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/dna.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/dna_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_barriers.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_factors.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_factors_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/setup.cpp)
target_link_libraries(
  libmodle_internal
  PRIVATE project_options project_warnings
  PUBLIC Modle::common
         Modle::io
         CONAN_PKG::abseil
         CONAN_PKG::xxhash)

target_include_directories(libmodle_internal PUBLIC include/)
