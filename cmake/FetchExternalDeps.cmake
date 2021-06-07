include(FetchContent)
FetchContent_Declare(
  libBigWig
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/libBigWig-0.4.6.tar.xz
  URL_HASH
    SHA512=5347edca4aeaf5c1ddc9c6e1c4b0383bf8bb40f1ea6d60253e67eeabc43e7f70383afab8faafab16622cf02a0ca761452df75ab322f0c29216b7b898a75f2279
)
FetchContent_Declare(
  Xoshiro-cpp
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/Xoshiro-cpp-1.1.tar.xz
  URL_HASH
    SHA512=95f2af121f6b062f13c1ed150508bcacc36e399839c923672f13971f54cfe136067d630539a204c828f54c77a04eba325f5f0f5a65b6e3146ef9cd75a6bafec9
)

FetchContent_GetProperties(libBigWig)
FetchContent_GetProperties(Xoshiro)

if(NOT ${libbigwig}_POPULATED)
  FetchContent_Populate(libBigWig)
endif()
add_subdirectory(${libbigwig_SOURCE_DIR} ${libbigwig_BINARY_DIR})

if(ENABLE_XOSHIRO)
  target_compile_definitions(project_options INTERFACE USE_XOSHIRO=1)
  if(NOT ${xoshiro-cpp}_POPULATED)
    FetchContent_Populate(Xoshiro-cpp)
  endif()
  add_subdirectory(${xoshiro-cpp_SOURCE_DIR} ${xoshiro-cpp_BINARY_DIR})
endif()
