# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl QUIET REQUIRED)
find_package(Boost QUIET REQUIRED)
find_package(cpp-sort QUIET REQUIRED)
find_package(fmt QUIET REQUIRED)
find_package(xxHash QUIET REQUIRED)

add_library(cmatrix INTERFACE)
add_library(Modle::cmatrix ALIAS cmatrix)

target_sources(
  cmatrix INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/contacts_impl.hpp ${CMAKE_CURRENT_SOURCE_DIR}/contacts_safe_impl.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/contacts_unsafe_impl.hpp)

target_include_directories(cmatrix INTERFACE include/)

target_link_libraries(
  cmatrix
  INTERFACE project_warnings
            project_options
            absl::span
            absl::strings
            Boost::filesystem
            Boost::headers
            cpp-sort::cpp-sort
            fmt::fmt
            xxHash::xxHash
            Modle::common
            Modle::io_compressed
            Modle::stats
            thread-pool)

if(ENABLE_PCH)
  target_precompile_headers(
    cmatrix
    INTERFACE
    <cstdint>
    <utility>
    <vector>
    <boost/dynamic_bitset/dynamic_bitset.hpp>
    [["modle/contacts.hpp"]])
endif()
