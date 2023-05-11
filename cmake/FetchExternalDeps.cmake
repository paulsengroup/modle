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
        URL ${CMAKE_CURRENT_SOURCE_DIR}/external/coolerpp-0b19a5b.tar.xz
        URL_HASH SHA512=b8df4b89afba941851573033ddf72a5725faadc3476ff4ac9a1d9abb7070d6a503be4bdf91ace0fa0ec9b87239d801d45f0bc6170b284cc21b99f45dda55e507
)
# cmake-format: on

set(COOLERPP_ENABLE_TESTING OFF)
set(COOLERPP_BUILD_EXAMPLES OFF)

FetchContent_MakeAvailable(coolerpp)
