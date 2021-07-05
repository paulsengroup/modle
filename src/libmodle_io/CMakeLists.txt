add_library(libmodle_io)
add_library(Modle::libmodle_io ALIAS libmodle_io)
target_sources(
  libmodle_io
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/bed.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/bigwig.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/bigwig_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/chrom_sizes.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/cooler.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/cooler_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/hdf5.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/hdf5_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/compressed_io.cpp)

add_dependencies(libmodle_io libBigWig)
target_include_directories(libmodle_io PUBLIC include/)

target_link_libraries(
  libmodle_io
  PRIVATE project_warnings
          project_options
          Modle::common
          CONAN_PKG::xxhash
  PUBLIC libBigWig::libbigwig
         Modle::cmatrix
         Modle::interval_tree
         CONAN_PKG::abseil
         CONAN_PKG::boost
         CONAN_PKG::fmt
         CONAN_PKG::hdf5
         CONAN_PKG::libarchive
         CONAN_PKG::readerwriterqueue
         CONAN_PKG::zlib
         # Deal with boost::filesystem quirkiness on older compilers
         $<$<AND:$<BOOL:${USING_LIBCXX}>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lc++fs>
         $<$<AND:$<BOOL:${USING_STDLIBCXX}>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lstdc++fs>)

target_compile_definitions(libmodle_io PRIVATE -DNOCURL)
set_target_properties(libmodle_io PROPERTIES OUTPUT_NAME modle_io)

if(ENABLE_PCH)
  target_precompile_headers(
    libmodle_io
    PUBLIC
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
    <H5Cpp.h>)
endif()

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()

install(
  TARGETS libmodle_io
  ARCHIVE
  INCLUDES
  LIBRARY)
