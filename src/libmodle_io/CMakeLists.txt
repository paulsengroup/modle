# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(LibArchive CONFIG REQUIRED)
find_package(BigWig CONFIG REQUIRED)
find_package(readerwriterqueue CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(ZLIB CONFIG REQUIRED)
find_package(xxHash CONFIG REQUIRED)

# BED
add_library(libmodle_io_bed)
add_library(Modle::io_bed ALIAS libmodle_io_bed)
target_sources(libmodle_io_bed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/bed.cpp ${CMAKE_CURRENT_SOURCE_DIR}/bed_impl.hpp)

target_include_directories(libmodle_io_bed PUBLIC include/bed)

target_link_libraries(
  libmodle_io_bed
  PRIVATE project_warnings project_options
  PUBLIC Modle::common Modle::interval_tree Modle::io_compressed)

target_link_system_libraries(
  libmodle_io_bed
  PRIVATE
  absl::strings
  PUBLIC
  absl::btree
  absl::flat_hash_map
  absl::hash
  absl::span
  xxHash::xxHash)

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
  PUBLIC Modle::common fmt::fmt)

target_link_system_libraries(
  libmodle_io_bigwig
  PUBLIC
  absl::span
  BigWig::BigWig)

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
  PRIVATE project_warnings project_options Modle::io_bed
  PUBLIC Modle::common Modle::io_compressed fmt::fmt)

target_link_system_libraries(
  libmodle_io_chrom_sizes
  PRIVATE
  absl::flat_hash_set
  absl::strings)

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
  PRIVATE project_warnings project_options
  PUBLIC Modle::common fmt::fmt)

target_link_system_libraries(
  libmodle_io_compressed
  PRIVATE
  project_warnings
  project_options
  absl::strings
  PUBLIC
  Boost::iostreams
  LibArchive::LibArchive)

set_target_properties(libmodle_io_compressed PROPERTIES OUTPUT_NAME modle_io_compressed)
install(
  TARGETS libmodle_io_compressed
  ARCHIVE
  INCLUDES
  LIBRARY)

# Contact matrix IO
add_library(libmodle_io_contact_matrix INTERFACE)
add_library(Modle::io_contact_matrix ALIAS libmodle_io_contact_matrix)
target_sources(libmodle_io_contact_matrix INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/contact_matrix_dense_io_impl.hpp)

target_include_directories(libmodle_io_contact_matrix INTERFACE include/contact_matrix_io)

target_link_libraries(
  libmodle_io_contact_matrix
  INTERFACE project_warnings
            project_options
            Modle::cmatrix
            Modle::common
            fmt::fmt)

target_link_system_libraries(
  libmodle_io_contact_matrix
  INTERFACE
  project_warnings
  project_options
  coolerpp)

set_target_properties(libmodle_io_contact_matrix PROPERTIES OUTPUT_NAME modle_io_contact_matrix)
install(
  TARGETS libmodle_io_contact_matrix
  ARCHIVE
  INCLUDES
  LIBRARY)
