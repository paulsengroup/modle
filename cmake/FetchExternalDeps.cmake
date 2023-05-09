# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

# cmake-format: off
FetchContent_Declare(
        coolerpp
        GIT_REPOSITORY https://github.com/robomics/coolerpp.git
        GIT_TAG 2e49aba0ef56184bc71879c551f71b8b1e45061b
)
# cmake-format: on

set(COOLERPP_ENABLE_TESTING OFF)
set(COOLERPP_BUILD_EXAMPLES OFF)

FetchContent_MakeAvailable(coolerpp)
