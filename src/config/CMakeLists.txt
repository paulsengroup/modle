# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)

add_library(modle_config)

add_library(Modle::config ALIAS modle_config)

add_dependencies(modle_config check_git)

# Sometimes the hook that is supposed to generate these files (usually git.cpp) does not work. Other times the file is
# created at config time, but is not updated at build time. Touching both files seems to address both issues.
file(TOUCH ${CMAKE_CURRENT_BINARY_DIR}/git.cpp ${CMAKE_CURRENT_BINARY_DIR}/version.cpp)

target_sources(
  modle_config
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/config/git.hpp
         ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/config/version.hpp
  PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/git.cpp ${CMAKE_CURRENT_BINARY_DIR}/version.cpp)

target_link_libraries(
  modle_config
  PRIVATE project_options
          project_warnings
          Modle::common
          fmt::fmt)
target_link_system_libraries(modle_config PUBLIC absl::strings)

target_include_directories(modle_config PUBLIC include/)
