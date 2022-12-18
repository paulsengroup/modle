# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)

# cmake-format: off
FetchContent_Declare(
        coolerpp
        GIT_REPOSITORY https://github.com/robomics/coolerpp.git
        GIT_TAG 54ee8f5063eab4a4f1a54f113cb01d6b2c01a499
)
# cmake-format: on

SET(COOLERPP_ENABLE_TESTING OFF)
SET(COOLERPP_BUILD_EXAMPLES OFF)
FetchContent_MakeAvailable(coolerpp)
