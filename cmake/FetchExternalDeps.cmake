include(FetchContent)
FetchContent_Declare(
  cgranges
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/cgranges-10-06-2021.tar.xz
  URL_HASH
    SHA512=287006ec02a9e8b4e67f5bcf31b4c039123e33f4ff4e14405e44b053c24556b1117c02cd2a1c2c4c205f81a1989852ae2b7219ac8be0af92275f9dc9747303d2
)
FetchContent_Declare(
  libBigWig
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/libBigWig-0.4.6.tar.xz
  URL_HASH
    SHA512=ec99f5a1e5cecdd85444e7960fa27c12bf0bfef4d997bd268b7272b4e0310996e206e444f65cf35dc856d9debb2d6ec9f321de87941af50850c8a139d5d4f85a
)
FetchContent_Declare(
  Xoshiro-cpp
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/Xoshiro-cpp-1.1.tar.xz
  URL_HASH
    SHA512=b8f7786f56733a284cf678142d0df232dba80d46e427d80bc360cc853c8617d80da3c49bf4654fea6635119fb4b86532a777ca66226262363537bd8625530354
)

FetchContent_GetProperties(cgranges)
FetchContent_GetProperties(libBigWig)
FetchContent_GetProperties(Xoshiro)

if(MODLE_BUILD_UTILS)
  if(NOT ${cgranges}_POPULATED)
    FetchContent_Populate(cgranges)
  endif()
  add_subdirectory(${cgranges_SOURCE_DIR} ${cgranges_BINARY_DIR})
endif()

if(NOT ${libbigwig}_POPULATED)
  FetchContent_Populate(libBigWig)
endif()
add_subdirectory(${libbigwig_SOURCE_DIR} ${libbigwig_BINARY_DIR})

if(NOT USE_MERSENNE_TWISTER)
  if(NOT ${xoshiro-cpp}_POPULATED)
    FetchContent_Populate(Xoshiro-cpp)
  endif()
  add_subdirectory(${xoshiro-cpp_SOURCE_DIR} ${xoshiro-cpp_BINARY_DIR})
else()
  target_compile_definitions(project_options INTERFACE USE_MERSENNE_TWISTER=1)
endif()
