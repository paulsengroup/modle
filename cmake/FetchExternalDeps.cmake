# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)

# cmake-format: off
FetchContent_Declare(
  libBigWig
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/libBigWig-0.4.6.tar.xz
  URL_HASH SHA512=11cfb35da7fa99fe8f73d219654d2cfb838a74b7c6ba5610b826265251bdfcfb4feb2c6e2fc2377edd73154dd971f9604ea0e270bd5830d1f28509f84ad49f7e
)
FetchContent_Declare(
  mscharconv
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/mscharconv.tar.xz
  URL_HASH SHA512=4378f5be5336c726c3c9d104941a45379484519cb34d61001c408c496b47b53547c9af3a82c9cd0eb332df0c97d244e2f6a20a34aa437cb304980e567e364c2c
)
FetchContent_Declare(
  thread-pool
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/thread-pool-2.0.0.tar.xz
  URL_HASH SHA512=71fa39216842c4759a6eb1a68b37b30a92c017c3753ad21869010aa4a898ea573aeaec85deb641fbb837b9fe1a4641813d97db72f89430abeb24663be8bda4dd
)
FetchContent_Declare(
  Xoshiro-cpp
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/Xoshiro-cpp-1.1.tar.xz
  URL_HASH SHA512=fb584cae675ebdb181801237a1462b0931478cb3123987b06dee8cbb4b6d823fcfa148f38aef184dd3192c985f6fe1984339bb2d5e1399db40501ab81a92ecfb
)
# cmake-format: on

FetchContent_GetProperties(libBigWig)
FetchContent_GetProperties(mscharconv)
FetchContent_GetProperties(thread-pool)
FetchContent_GetProperties(Xoshiro)

set(WITH_CURL
    OFF
    CACHE INTERNAL "")
if(NOT ${libbigwig}_POPULATED)
  FetchContent_Populate(libBigWig)
endif()
add_subdirectory(${libbigwig_SOURCE_DIR} ${libbigwig_BINARY_DIR} EXCLUDE_FROM_ALL)

if(NOT ${mscharconv}_POPULATED)
  FetchContent_Populate(mscharconv)
endif()
add_subdirectory(${mscharconv_SOURCE_DIR} ${mscharconv_BINARY_DIR})

if(NOT ${thread-pool}_POPULATED)
  FetchContent_Populate(thread-pool)
endif()
add_subdirectory(${thread-pool_SOURCE_DIR} ${thread-pool_BINARY_DIR})

if(NOT MODLE_USE_MERSENNE_TWISTER)
  if(NOT ${xoshiro-cpp}_POPULATED)
    FetchContent_Populate(Xoshiro-cpp)
  endif()
  add_subdirectory(${xoshiro-cpp_SOURCE_DIR} ${xoshiro-cpp_BINARY_DIR})
else()
  target_compile_definitions(project_options INTERFACE MODLE_USE_MERSENNE_TWISTER=1)
endif()
