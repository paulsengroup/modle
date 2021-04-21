add_library(libio INTERFACE)
add_library(Modle::io ALIAS libio)
target_sources(
  libio
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/bed_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/bigwig_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/chrom_sizes_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/cooler_impl.hpp)
add_dependencies(libio libBigWig_lib)

target_include_directories(libio INTERFACE include/)

target_link_libraries(
  libio
  INTERFACE project_warnings
            project_options
            Modle::cmatrix
            libBigWig::libBigWig
            CONAN_PKG::abseil
            CONAN_PKG::boost
            CONAN_PKG::fmt
            CONAN_PKG::hdf5
            CONAN_PKG::readerwriterqueue
            CONAN_PKG::zlib
            # Deal with std::filesystem quirkiness on older compilers
            $<$<AND:$<BOOL:${USING_LIBCXX}>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lc++fs>
            $<$<AND:$<BOOL:${USING_STDLIBCXX}>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lstdc++fs>)

target_compile_definitions(libio INTERFACE -DNOCURL)

if(ENABLE_PCH)
  target_precompile_headers(
    libio
    INTERFACE
    <algorithm>
    <cassert>
    <cstddef>
    <cstdint>
    <filesystem>
    <fstream>
    <limits>
    <memory>
    <stdexcept>
    <string>
    <string_view>
    <vector>
    <utility>
    <absl/container/flat_hash_map.h>
    <fmt/format.h>
    <H5Cpp.h>
    [["modle/bed.hpp"]]
    [["modle/bigwig.hpp"]]
    [["modle/chrom_sizes.hpp"]]
    [["modle/cooler.hpp"]]
    [["modle/hdf5.hpp"]])
endif()

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()
