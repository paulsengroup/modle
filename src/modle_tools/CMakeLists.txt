add_executable(modle_tools main.cpp)

target_sources(
  modle_tools
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cli.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/include/modle_tools/config.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/cli.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/eval.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/filter_barriers.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/find_barrier_clusters.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/noisify.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/stats.cpp)

target_link_libraries(
  modle_tools
  PRIVATE project_warnings
          project_options
          Modle::common
          Modle::correlation
          Modle::libmodle_io
          CONAN_PKG::abseil
          CONAN_PKG::boost
          CONAN_PKG::cli11
          CONAN_PKG::readerwriterqueue)

target_include_directories(modle_tools PRIVATE include)

if(ENABLE_PCH)
  target_precompile_headers(
    modle_tools
    PRIVATE
    <algorithm>
    <cassert>
    <cstddef>
    <cstdint>
    <filesystem>
    <stdexcept>
    <string>
    <string_view>
    <utility>
    <vector>
    <absl/container/flat_hash_set.h>
    <absl/types/span.h>
    <fmt/format.h>)
endif()

install(
  TARGETS modle_tools
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)
