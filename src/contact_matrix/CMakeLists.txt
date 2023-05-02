# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED COMPONENTS iostreams serialization)
find_package(bshoshany-thread-pool CONFIG REQUIRED)
find_package(cpp-sort CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(libcuckoo CONFIG REQUIRED)
find_package(xxHash CONFIG REQUIRED)

add_library(cmatrix INTERFACE)
add_library(Modle::cmatrix ALIAS cmatrix)

target_sources(
  cmatrix
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/contact_matrix_dense_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/contact_matrix_dense_safe_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/contact_matrix_dense_unsafe_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/contact_matrix_internal_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/contact_matrix_sparse_impl.hpp)

target_include_directories(cmatrix INTERFACE include/)

target_link_libraries(
  cmatrix
  INTERFACE project_warnings
            project_options
            Modle::common
            Modle::io_compressed
            Modle::stats
            fmt::fmt)

target_link_system_libraries(
  cmatrix
  INTERFACE
  absl::span
  absl::strings
  Boost::headers
  Boost::iostreams
  Boost::serialization
  bshoshany-thread-pool::bshoshany-thread-pool
  cpp-sort::cpp-sort
  libcuckoo::libcuckoo
  xxHash::xxhash)
