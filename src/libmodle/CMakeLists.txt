add_library(libmodle INTERFACE)
add_library(Modle::Modle ALIAS libmodle)

target_sources(
  libmodle
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/dna_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_barriers_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_factors_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_impl.hpp)

target_include_directories(libmodle INTERFACE include/)

target_link_libraries(
  libmodle
  INTERFACE
    project_warnings
    project_options
    Modle::cmatrix
    Modle::io
    $<$<BOOL:${ENABLE_XOSHIRO}>:XoshiroCpp::XoshiroCpp>
    CONAN_PKG::abseil
    CONAN_PKG::boost
    CONAN_PKG::concurrentqueue
    CONAN_PKG::cpp-sort
    # Deal with std::filesystem quirkiness on older compilers
#    $<$<BOOL:${STDLIBCXX_REQUIRES_LSTDCXXFS_LINK_FLAG}>:-lstdc++fs>
#    $<$<BOOL:${LIBCXX_REQUIRES_LCXXFS_LINK_FLAG}>:-lc++fs>)
  $<$<AND:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lc++fs>
  $<$<AND:$<CXX_COMPILER_ID:GNU>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lstdc++fs>)

if(ENABLE_PCH)
  target_precompile_headers(
    libmodle
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
