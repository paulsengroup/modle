# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)

add_library(modle_config)

add_library(Modle::config ALIAS modle_config)

target_sources(
  modle_config
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/config/git.hpp
         ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/config/version.hpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/git.cpp ${CMAKE_CURRENT_SOURCE_DIR}/version.cpp)

target_link_libraries(modle_config PRIVATE project_options project_warnings Modle::common)
target_link_system_libraries(
  modle_config
  PUBLIC
  absl::strings
  fmt::fmt)

target_include_directories(modle_config PUBLIC include/)
