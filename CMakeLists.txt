cmake_minimum_required(VERSION 3.12)
project(modle)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(FetchContent)
include(CheckIPOSupported)

FetchContent_Declare(
        abseil-cpp
        GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
        GIT_TAG 20200923
        GIT_SHALLOW ON
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/abseil-cpp
)

FetchContent_MakeAvailable(abseil-cpp)

set(MODLE_DEBUG_OPTIONS "-Wall;-Wextra;-ggdb;-Og;-g3")
set(MODLE_RELEASE_OPTIONS "-Wall;-Wextra;-O3")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

check_ipo_supported(RESULT IPO_AVAILABLE)
if (IPO_AVAILABLE)
    message(STATUS "IPO/LTO available")
else ()
    message(STATUS "IPO/LTO not available")
endif ()

IF (CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "Debug mode detected: enabling tests")
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/include/)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)
endif ()

add_library(
        libmodle
        src/config.cpp
        src/dna.cpp
        src/extr_barrier.cpp
        src/genome.cpp
        src/lefs.cpp
        src/parsers.cpp
        src/contacts.cpp
)

target_include_directories(
        libmodle
        SYSTEM
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/include/
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/abseil-cpp/include/
)

target_link_libraries(
        libmodle
        absl::base
        absl::flat_hash_map
        absl::strings
        absl::str_format
        absl::time
)


target_compile_options(libmodle PUBLIC "$<$<CONFIG:DEBUG>:${MODLE_DEBUG_OPTIONS}>")
target_compile_options(libmodle PUBLIC "$<$<CONFIG:RELEASE>:${MODLE_RELEASE_OPTIONS}>")

add_executable(modle src/main.cpp)

target_link_libraries(
        modle
        libmodle
        absl::base
        absl::strings
        absl::str_format
        absl::time
)

target_include_directories(
        modle
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/include/
        SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/abseil-cpp/include/
)

target_compile_options(modle PUBLIC "$<$<CONFIG:DEBUG>:${MODLE_DEBUG_OPTIONS}>")
target_compile_options(modle PUBLIC "$<$<CONFIG:RELEASE>:${MODLE_RELEASE_OPTIONS}>")

