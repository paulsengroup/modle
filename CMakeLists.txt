cmake_minimum_required(VERSION 3.12)
project(modle
        HOMEPAGE_URL "https://github.com/robomics/modle"
        LANGUAGES CXX
        VERSION 0.0.1
        DESCRIPTION "Software for in-silico modeling of DNA loop extrusion.")

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to \"Release\" as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif ()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if (NOT MODLE_BUILD_DOCS)
    set(MODLE_BUILD_DOCS OFF)
endif ()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
include(CheckIPOSupported)

find_package(ZLIB REQUIRED)

FetchContent_Declare(
        abseil-cpp
        GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
        GIT_TAG 20200923
        GIT_SHALLOW ON
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/abseil-cpp
)

FetchContent_Declare(
        CLI11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG v1.9.1
        GIT_SHALLOW ON
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/CLI11
)

FetchContent_MakeAvailable(abseil-cpp CLI11)

if (MODLE_BUILD_DOCS)
    message(STATUS "Will build ModLE's documentation")
    add_subdirectory(docs)
endif ()

set(MODLE_DEBUG_OPTIONS "-Wall;-Wextra;-ggdb;-Og;-g3")
set(MODLE_RELEASE_OPTIONS "-Wall;-Wextra;-O3")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

check_ipo_supported(RESULT IPO_AVAILABLE)
if (IPO_AVAILABLE)
    message(STATUS "IPO/LTO available")
else ()
    message(STATUS "IPO/LTO not available")
endif ()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Debug mode detected: enabling tests")
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/include/)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)
endif ()

add_library(
        libmodle
        src/dna.cpp
        src/extr_barrier.cpp
        src/genome.cpp
        src/lefs.cpp
        src/parsers.cpp
)

target_include_directories(
        libmodle
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/include/
        PRIVATE ${ZLIB_INCLUDE_DIRS}/
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/abseil-cpp/include/
)

target_link_libraries(
        libmodle
        absl::base
        absl::flat_hash_map
        absl::strings
        absl::str_format
        absl::time
        ZLIB::ZLIB
)

set_target_properties(libmodle
        PROPERTIES OUTPUT_NAME modle)


target_compile_options(libmodle PUBLIC "$<$<CONFIG:DEBUG>:${MODLE_DEBUG_OPTIONS}>")
target_compile_options(libmodle PUBLIC "$<$<CONFIG:RELEASE>:${MODLE_RELEASE_OPTIONS}>")

add_executable(
        modle
        src/main.cpp
)

target_link_libraries(
        modle
        libmodle
        absl::base
        absl::strings
        absl::str_format
        absl::time
        CLI11::CLI11
)

target_include_directories(
        modle
        PUBLIC libmodle
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/CLI11/include/
)

target_compile_options(modle PUBLIC "$<$<CONFIG:DEBUG>:${MODLE_DEBUG_OPTIONS}>")
target_compile_options(modle PUBLIC "$<$<CONFIG:RELEASE>:${MODLE_RELEASE_OPTIONS}>")

install(TARGETS libmodle
        LIBRARY
        COMPONENT Libraries
        NAMELINK_COMPONENT Development
        PUBLIC_HEADER
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/modle
        COMPONENT Development)

install(TARGETS modle
        DESTINATION ${CMAKE_INSTALL_BINDIR})
