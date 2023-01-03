# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)

# cmake-format: off
FetchContent_Declare(
        coolerpp
        GIT_REPOSITORY https://github.com/robomics/coolerpp.git
        GIT_TAG 7349bf97d92e16bb173b1591c1e399fdca1f3ab2
)
# cmake-format: on

set(COOLERPP_ENABLE_TESTING OFF)
set(COOLERPP_BUILD_EXAMPLES OFF)

if(BUILD_SHARED)
  FetchContent_MakeAvailable(coolerpp)
else()
  # Do not install coolerpp when using static linking
  FetchContent_GetProperties(coolerpp)
  if(NOT coolerpp_POPULATED)
    FetchContent_Populate(coolerpp)
    add_subdirectory(
      ${coolerpp_SOURCE_DIR}
      ${coolerpp_BINARY_DIR}
      EXCLUDE_FROM_ALL
      SYSTEM)
  endif()
endif()
