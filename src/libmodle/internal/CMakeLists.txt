find_package(absl REQUIRED)
find_package(Boost REQUIRED)
find_package(fmt REQUIRED)
find_package(spdlog REQUIRED)
find_package(xxHash REQUIRED)

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
            absl::btree
            absl::optional
            absl::span
            Boost::filesystem
            fmt::fmt
            xxHash::xxhash
            spdlog::spdlog
            Modle::common
            Modle::cmatrix
            Modle::libmodle_io)

target_include_directories(libmodle_internal INTERFACE include/)

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()
