# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)

# cmake-format: off
FetchContent_Declare(
        coolerpp
        GIT_REPOSITORY https://github.com/robomics/coolerpp.git
        GIT_TAG d12325ae25121524d514e952c6f35a6350ee7207
)
# cmake-format: on

set(COOLERPP_ENABLE_TESTING OFF)
set(COOLERPP_BUILD_EXAMPLES OFF)

FetchContent_MakeAvailable(coolerpp)
