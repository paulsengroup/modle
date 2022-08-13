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
        URL https://github.com/dpryan79/libBigWig/archive/refs/tags/0.4.7.tar.gz
        URL_HASH SHA512=52f1b7c8e21e16238b3bb07baef6aa3611797b1b5ff44d912c874f8f4527c516a0676877fad21c103b8b25a733e84bef48530f28dc224a79d43f7764eae7ed40
)
# cmake-format: on

set(WITH_CURL
    OFF
    CACHE INTERNAL "")

FetchContent_MakeAvailable(bitflags libBigWig)
