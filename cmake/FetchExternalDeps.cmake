include(FetchContent)

# cmake-format: off
FetchContent_Declare(
  libBigWig
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/libBigWig-0.4.6.tar.xz
  URL_HASH SHA512=a76cc69e4d3b5f8a42a801d304b6559837723d0a4ff9a13aaa85717e5d1406d5adfc0b2fa44e7e7ba31b5da5909004d71e1be96d209152d30fbfbdd6330f72e2
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
FetchContent_GetProperties(thread-pool)
FetchContent_GetProperties(Xoshiro)

set(WITH_CURL
    OFF
    CACHE INTERNAL "")
if(NOT ${libbigwig}_POPULATED)
  FetchContent_Populate(libBigWig)
endif()
add_subdirectory(${libbigwig_SOURCE_DIR} ${libbigwig_BINARY_DIR} EXCLUDE_FROM_ALL)

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
