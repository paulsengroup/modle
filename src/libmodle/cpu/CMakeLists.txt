add_library(libmodle_cpu INTERFACE)

target_sources(
  libmodle_cpu
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/dna_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_barriers_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_factors_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_correct_moves_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_detect_collisions_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_scheduler_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_setup_impl.hpp)

target_include_directories(libmodle_cpu INTERFACE include/)

target_link_libraries(
  libmodle_cpu
  INTERFACE project_warnings
            project_options
            Modle::cmatrix
            Modle::config
            Modle::io
            Modle::modle_io_internal
            $<$<BOOL:${ENABLE_XOSHIRO}>:XoshiroCpp::XoshiroCpp>
            CONAN_PKG::abseil
            CONAN_PKG::boost
            CONAN_PKG::concurrentqueue
            CONAN_PKG::cpp-sort
            CONAN_PKG::xxhash
            # Deal with std::filesystem quirkiness on older compilers
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
    <fmt/format.h>
    [["modle/common.hpp"]]
    [["modle/config.hpp"]])
endif()

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()
