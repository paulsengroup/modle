# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)

# cmake-format: off
FetchContent_Declare(
        coolerpp
        GIT_REPOSITORY https://github.com/robomics/coolerpp.git
        GIT_TAG 9e03846216b0dcae96c190ccc24b60ed23c16fde
)
# cmake-format: on

SET(ENABLE_TESTING OFF)
SET(BUILD_EXAMPLES OFF)
FetchContent_MakeAvailable(coolerpp)
