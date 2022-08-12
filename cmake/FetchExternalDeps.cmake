# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)

# cmake-format: off
FetchContent_Declare(
        bitflags
        URL ${CMAKE_CURRENT_SOURCE_DIR}/external/bitflags-1.5.0.tar.xz
        URL_HASH SHA512=918b73fd40ce6180c237caae221c0d8bea74b203d75a77ee2c399cbf1e063894f1df4a838d13fd87ca0870551ff83449462f8705e61e1c21c7f3c1e47ba07b71
)
FetchContent_Declare(
        libBigWig
        GIT_REPOSITORY https://github.com/dpryan79/libBigWig.git
        GIT_TAG b77db88cb8d3f3bfa27df5d0c4f3623fd4445046
)
# cmake-format: on

set(WITH_CURL
    OFF
    CACHE INTERNAL "")

FetchContent_MakeAvailable(bitflags libBigWig)
