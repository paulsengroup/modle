include(FetchContent)

# cmake-format: off
FetchContent_Declare(
  libBigWig
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/libBigWig-0.4.6.tar.xz
  URL_HASH SHA512=2708f5c14966749b3fe0d80196acd229c48b47833c5d3af2b2d456393ee14cbbf454faa9970b98f4b4793564e45194ab47ea4aeb450206be647997ad50de0262
)
FetchContent_Declare(
  thread-pool
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/thread-pool-2.0.0.tar.xz
  URL_HASH SHA512=6a24e20870361c24a2f104356c7b772c9628355c8f7b58fa961e67b7e67e75369dd39ab2ad49b41b93131662fe075f3d750216f7d93cff15c19a347d7cc865dd
)
FetchContent_Declare(
  Xoshiro-cpp
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/Xoshiro-cpp-1.1.tar.xz
  URL_HASH SHA512=b8f7786f56733a284cf678142d0df232dba80d46e427d80bc360cc853c8617d80da3c49bf4654fea6635119fb4b86532a777ca66226262363537bd8625530354
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
