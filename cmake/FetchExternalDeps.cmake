# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

# cmake-format: off
FetchContent_Declare(
        hictk
        URL ${CMAKE_CURRENT_SOURCE_DIR}/external/hictk-24f75ae.tar.xz
        URL_HASH SHA256=452519a7b4e6d64f45d84ac6eacff9fc7ef119fa4aeb10e5d57ca6251710969c
)
# cmake-format: on

set(HICTK_ENABLE_TESTING OFF)
set(HICTK_BUILD_EXAMPLES OFF)
set(HICTK_BUILD_BENCHMARKS OFF)
set(HICTK_BUILD_TOOLS OFF)

FetchContent_MakeAvailable(hictk)
