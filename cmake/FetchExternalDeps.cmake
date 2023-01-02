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
FetchContent_MakeAvailable(coolerpp)
