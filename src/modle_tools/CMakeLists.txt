add_executable(modle_tools main.cpp)

target_sources(
  modle_tools
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cli.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/eval_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/noisify_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/tools_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/stats_impl.hpp)

target_link_libraries(
  modle_tools
  PRIVATE project_warnings
          project_options
          Modle::modle
          Modle::correlation
          Modle::io
          CGRanges::CGRanges
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
    <fmt/format.h>
    [["modle_tools/config.hpp"]]
    [["modle_tools/eval.hpp"]]
    [["modle_tools/noisify.hpp"]]
    [["modle_tools/stats.hpp"]]
    [["modle_tools/tools.hpp"]])
endif()

install(
  TARGETS modle_tools
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)
