# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl QUIET REQUIRED)
find_package(Boost QUIET REQUIRED)
find_package(fmt QUIET REQUIRED)
find_package(HDF5 QUIET REQUIRED)
find_package(LibArchive QUIET REQUIRED)
find_package(readerwriterqueue QUIET REQUIRED)
find_package(spdlog QUIET REQUIRED)
find_package(ZLIB QUIET REQUIRED)
find_package(xxHash QUIET REQUIRED)

# BED
add_library(libmodle_io_bed)
add_library(Modle::io_bed ALIAS libmodle_io_bed)
target_sources(libmodle_io_bed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/bed.cpp ${CMAKE_CURRENT_SOURCE_DIR}/bed_impl.hpp)

target_include_directories(libmodle_io_bed PUBLIC include/bed)

target_link_libraries(
  libmodle_io_bed
  PRIVATE project_warnings project_options absl::strings
  PUBLIC absl::btree
         absl::flat_hash_map
         absl::hash
         absl::span
         Boost::filesystem
         fmt::fmt
         xxHash::xxHash
         Modle::common
         Modle::interval_tree
         Modle::io_compressed)

set_target_properties(libmodle_io_bed PROPERTIES OUTPUT_NAME modle_io_bed)
install(
  TARGETS libmodle_io_bed
  ARCHIVE
  INCLUDES
  LIBRARY)

# BigWig
add_library(libmodle_io_bigwig)
add_library(Modle::io_bigwig ALIAS libmodle_io_bigwig)
target_sources(libmodle_io_bigwig PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/bigwig.cpp
                                          ${CMAKE_CURRENT_SOURCE_DIR}/bigwig_impl.hpp)

target_include_directories(libmodle_io_bigwig PUBLIC include/bigwig)

target_link_libraries(
  libmodle_io_bigwig
  PRIVATE project_warnings project_options
  PUBLIC absl::span
         Boost::filesystem
         fmt::fmt
         Modle::common
         libBigWig::libbigwig)

add_dependencies(libmodle_io_bigwig libBigWig)
target_compile_definitions(libmodle_io_bigwig PRIVATE -DNOCURL)

set_target_properties(libmodle_io_bigwig PROPERTIES OUTPUT_NAME modle_io_bigwig)
install(
  TARGETS libmodle_io_bigwig
  ARCHIVE
  INCLUDES
  LIBRARY)

# ChromSizes
add_library(libmodle_io_chrom_sizes)
add_library(Modle::io_chrom_sizes ALIAS libmodle_io_chrom_sizes)
target_sources(libmodle_io_chrom_sizes PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/chrom_sizes.cpp)

target_include_directories(libmodle_io_chrom_sizes PUBLIC include/chrom_sizes)

target_link_libraries(
  libmodle_io_chrom_sizes
  PRIVATE project_warnings
          project_options
          absl::flat_hash_set
          absl::strings
          fmt::fmt
          Modle::io_bed
  PUBLIC Boost::filesystem Modle::common Modle::io_compressed)

set_target_properties(libmodle_io_chrom_sizes PROPERTIES OUTPUT_NAME modle_io_chrom_sizes)
install(
  TARGETS libmodle_io_chrom_sizes
  ARCHIVE
  INCLUDES
  LIBRARY)

# Compressed IO
add_library(libmodle_io_compressed)
add_library(Modle::io_compressed ALIAS libmodle_io_compressed)
target_sources(libmodle_io_compressed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/compressed_io.cpp)

target_include_directories(libmodle_io_compressed PUBLIC include/compressed_io)

target_link_libraries(
  libmodle_io_compressed
  PRIVATE project_warnings
          project_options
          absl::strings
          fmt::fmt
  PUBLIC Boost::filesystem
         Boost::iostreams
         LibArchive::LibArchive
         Modle::common)

set_target_properties(libmodle_io_compressed PROPERTIES OUTPUT_NAME modle_io_compressed)
install(
  TARGETS libmodle_io_compressed
  ARCHIVE
  INCLUDES
  LIBRARY)

# Cooler
add_library(libmodle_io_cooler INTERFACE)
add_library(Modle::io_cooler ALIAS libmodle_io_cooler)
target_sources(
  libmodle_io_cooler
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cooler_impl.hpp ${CMAKE_CURRENT_SOURCE_DIR}/cooler_read_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/cooler_write_impl.hpp)

target_include_directories(libmodle_io_cooler INTERFACE include/cooler)

target_link_libraries(
  libmodle_io_cooler
  INTERFACE project_warnings
            project_options
            absl::span
            absl::strings
            absl::time
            absl::variant
            Boost::filesystem
            fmt::fmt
            readerwriterqueue::readerwriterqueue
            spdlog::spdlog
            Modle::common
            Modle::io_hdf5)

install(
  TARGETS libmodle_io_cooler
  ARCHIVE
  INCLUDES
  LIBRARY)

# HDF5
add_library(libmodle_io_hdf5)
add_library(Modle::io_hdf5 ALIAS libmodle_io_hdf5)
target_sources(libmodle_io_hdf5 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/hdf5_impl.hpp ${CMAKE_CURRENT_SOURCE_DIR}/hdf5.cpp)

target_include_directories(libmodle_io_hdf5 PUBLIC include/hdf5)

target_link_libraries(
  libmodle_io_hdf5
  PRIVATE project_warnings project_options absl::strings
  PUBLIC absl::variant
         fmt::fmt
         Modle::common
         HDF5::HDF5)

set_target_properties(libmodle_io_hdf5 PROPERTIES OUTPUT_NAME modle_io_hdf5)
install(
  TARGETS libmodle_io_hdf5
  ARCHIVE
  INCLUDES
  LIBRARY)
