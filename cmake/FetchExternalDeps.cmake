# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)

# cmake-format: off
FetchContent_Declare(
        hictk
        URL ${CMAKE_CURRENT_SOURCE_DIR}/external/hictk-v2.1.4.tar.xz
        URL_HASH SHA256=21abacc4322408e4611341e47dd3198aa5c7e9177d0d39e6addb615adf238a97
        EXCLUDE_FROM_ALL
        OVERRIDE_FIND_PACKAGE
        SYSTEM
)
# cmake-format: on

set(HICTK_BUILD_BENCHMARKS OFF)
set(HICTK_BUILD_EXAMPLES OFF)
set(HICTK_BUILD_TOOLS OFF)
set(HICTK_ENABLE_FUZZY_TESTING OFF)
set(HICTK_ENABLE_TESTING OFF)
set(HICTK_INSTALL OFF)
set(HICTK_WITH_ARROW OFF)
set(HICTK_WITH_EIGEN OFF)
set(HICTK_ENABLE_GIT_VERSION_TRACKING OFF)

FetchContent_MakeAvailable(hictk)
