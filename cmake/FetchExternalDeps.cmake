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
        URL ${CMAKE_CURRENT_SOURCE_DIR}/external/hictk-7a1a5d7.tar.xz
        URL_HASH SHA256=ebbf369cc4ab01af64a4c320851b577d2fd6df636e87aeb36fb6c86524066e6e
)
# cmake-format: on

set(HICTK_ENABLE_TESTING OFF)
set(HICTK_BUILD_EXAMPLES OFF)
set(HICTK_BUILD_BENCHMARKS OFF)
set(HICTK_BUILD_TOOLS OFF)

FetchContent_MakeAvailable(hictk)
