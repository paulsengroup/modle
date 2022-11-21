# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)

# cmake-format: off
FetchContent_Declare(
        coolerpp
        GIT_REPOSITORY https://github.com/robomics/coolerpp.git
        GIT_TAG d4e5d6222a21feedea7c101469238b451144047c
)
# cmake-format: on

SET(ENABLE_TESTING OFF)
SET(BUILD_EXAMPLES OFF)
FetchContent_MakeAvailable(coolerpp)
