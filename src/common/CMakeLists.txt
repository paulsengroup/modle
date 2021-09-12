# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl REQUIRED)
find_package(Boost REQUIRED)
find_package(fmt REQUIRED)
find_package(xxHash REQUIRED)

add_library(modle_common INTERFACE)

add_library(Modle::common ALIAS modle_common)

target_sources(
  modle_common
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/common.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/config.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/genextreme_value_distribution.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/math.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/random.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/random_sampling.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/suppress_compiler_warnings.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/utils.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/math_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/utils_impl.hpp)

target_link_libraries(
  modle_common
  INTERFACE project_options
            project_warnings
            absl::strings
            absl::str_format
            Boost::headers
            Boost::filesystem
            Boost::random
            fmt::fmt
            $<$<NOT:$<BOOL:${MODLE_USE_MERSENNE_TWISTER}>>:XoshiroCpp::XoshiroCpp>
            xxHash::xxhash
            mscharconv)

target_include_directories(modle_common INTERFACE include/)
