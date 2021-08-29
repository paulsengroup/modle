find_package(absl REQUIRED)
find_package(Boost REQUIRED)
find_package(Catch2 REQUIRED)
find_package(fmt REQUIRED)
find_package(HDF5 REQUIRED)
find_package(LibArchive REQUIRED)
find_package(readerwriterqueue REQUIRED)
find_package(spdlog REQUIRED)
find_package(ZLIB REQUIRED)
find_package(xxHash REQUIRED)

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
          Boost::headers
          xxHash::xxhash
          Modle::common
  PUBLIC absl::btree
         absl::flat_hash_map
         absl::flat_hash_set
         absl::hash
         absl::span
         absl::strings
         absl::time
         absl::variant
         Catch2::Catch2
         Boost::filesystem
         Boost::iostreams
         fmt::fmt
         HDF5::HDF5
         LibArchive::LibArchive
         readerwriterqueue::readerwriterqueue
         spdlog::spdlog
         ZLIB::ZLIB
         Modle::cmatrix
         Modle::interval_tree
         libBigWig::libbigwig)

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

install(
  TARGETS libmodle_io
  ARCHIVE
  INCLUDES
  LIBRARY)
