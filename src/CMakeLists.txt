# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

option(MODLE_BUILD_UTILS "Enable compilation of MoDLE's utilities (such as modle_tools)" ON)

if(MODLE_BUILD_UTILS)
  message(STATUS "Building MoDLE utils")
  add_subdirectory(modle_tools)
endif()

add_subdirectory(common EXCLUDE_FROM_ALL)
add_subdirectory(contact_matrix EXCLUDE_FROM_ALL)
add_subdirectory(correlation EXCLUDE_FROM_ALL)
add_subdirectory(interval_tree EXCLUDE_FROM_ALL)
add_subdirectory(libmodle_io EXCLUDE_FROM_ALL)
add_subdirectory(libmodle EXCLUDE_FROM_ALL)
add_subdirectory(modle)
