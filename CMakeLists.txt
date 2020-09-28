cmake_minimum_required(VERSION 3.12)
project(modle)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(FetchContent)

FetchContent_Declare(abseil-cpp
        GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
        GIT_TAG        20200923
        GIT_SHALLOW    ON
        SOURCE_DIR     ${CMAKE_CURRENT_SOURCE_DIR}/external/abseil-cpp)

FetchContent_MakeAvailable(abseil-cpp)

# find_package(TBB 2020.1 REQUIRED)

add_executable(modle)

target_sources(modle
        PRIVATE src/main.cpp
        PRIVATE src/config.cpp
        PRIVATE src/dna.cpp
        PRIVATE src/extr_barrier.cpp
        PRIVATE src/genome.cpp
        PRIVATE src/lefs.cpp
        PRIVATE src/parsers.cpp
        )

target_include_directories(modle
        SYSTEM
        PUBLIC    ${CMAKE_CURRENT_SOURCE_DIR}/src/include/
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/abseil-cpp/include/)

target_link_libraries(modle
        absl::base
        absl::flat_hash_map
        absl::strings
        absl::str_format
        absl::time)