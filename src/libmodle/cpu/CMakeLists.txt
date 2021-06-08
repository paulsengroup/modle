add_library(libmodle_cpu STATIC)

target_sources(
  libmodle_cpu
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/simulation.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_correct_moves.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_detect_collisions.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_scheduler.cpp)

target_include_directories(libmodle_cpu PUBLIC include/)

target_link_libraries(
  libmodle_cpu
  PRIVATE project_warnings
          project_options
          Modle::cmatrix
          Modle::common
          CONAN_PKG::cpp-sort
  PUBLIC Modle::modle_internal
         CONAN_PKG::abseil
         CONAN_PKG::boost
         CONAN_PKG::concurrentqueue
         # Deal with std::filesystem quirkiness on older compilers
         $<$<AND:$<BOOL:${USING_LIBCXX}>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lc++fs>
         $<$<AND:$<BOOL:${USING_STDLIBCXX}>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lstdc++fs>)

set_target_properties(libmodle_cpu PROPERTIES OUTPUT_NAME modle_cpu)

if(ENABLE_PCH)
  target_precompile_headers(
    libmodle_cpu
    PUBLIC
    <algorithm>
    <cstddef>
    <cstdint>
    <filesystem>
    <limits>
    <memory>
    <string>
    <string_view>
    <utility>
    <vector>
    <absl/types/span.h>
    <fmt/format.h>)
endif()

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()
