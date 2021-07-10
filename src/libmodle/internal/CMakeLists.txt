add_library(libmodle_internal INTERFACE)
add_library(Modle::libmodle_internal ALIAS libmodle_internal)

target_sources(
  libmodle_internal
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/genome.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/genome_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_barriers.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_factors.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_factors_impl.hpp)
target_link_libraries(
  libmodle_internal
  INTERFACE project_options
            project_warnings
            Modle::common
            Modle::cmatrix
  INTERFACE Modle::libmodle_io CONAN_PKG::abseil CONAN_PKG::xxhash)

target_include_directories(libmodle_internal INTERFACE include/)

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()
