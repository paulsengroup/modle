add_library(libmodle_cpu INTERFACE)

target_sources(
  libmodle_cpu
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/simulation.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_correct_moves.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_detect_collisions.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_scheduler.cpp)

target_include_directories(libmodle_cpu INTERFACE include/)

target_link_libraries(
  libmodle_cpu
  INTERFACE project_warnings
            project_options
            Modle::cmatrix
            Modle::common
            CONAN_PKG::cpp-sort
  INTERFACE Modle::modle_internal
            CONAN_PKG::abseil
            CONAN_PKG::boost
            CONAN_PKG::concurrentqueue
            # Deal with boost::filesystem quirkiness on older compilers
            $<$<AND:$<BOOL:${USING_LIBCXX}>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lc++fs>
            $<$<AND:$<BOOL:${USING_STDLIBCXX}>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lstdc++fs>)

if(ENABLE_PCH)
  target_precompile_headers(
    libmodle_cpu
    INTERFACE
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
