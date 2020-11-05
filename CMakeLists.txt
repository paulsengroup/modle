cmake_minimum_required(VERSION 3.10)
project(modle
        HOMEPAGE_URL "https://github.com/robomics/modle"
        LANGUAGES CXX
        VERSION 0.0.1
        DESCRIPTION "Software for in-silico modeling of DNA loop extrusion.")

cmake_policy(VERSION 3.10)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

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

find_package(BZip2 REQUIRED)
find_package(TBB REQUIRED)
find_package(Boost REQUIRED COMPONENTS iostreams filesystem)
find_package(Threads REQUIRED)

FetchContent_Declare(
        abseil-cpp
        GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
        GIT_TAG "20200923.2"
        GIT_SHALLOW ON
        SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/abseil-cpp"
)

FetchContent_Declare(
        CLI11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG v1.9.1
        GIT_SHALLOW ON
        SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/CLI11"
)

FetchContent_MakeAvailable(abseil-cpp CLI11)

if (MODLE_BUILD_DOCS)
    message(STATUS "Will build ModLE's documentation")
    add_subdirectory(docs)
endif ()

set(MODLE_DEBUG_OPTIONS "-Wall;-Wextra;-ggdb;-Og;-g3")
set(MODLE_RELEASE_OPTIONS "-Wall;-Wextra;-O3")
set(MODLE_RELEASE_W_TSAN_OPTIONS "-Wall;-Wextra;-O2;-fsanitize=thread;-g")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

check_ipo_supported(RESULT IPO_AVAILABLE)
if (IPO_AVAILABLE)
    message(STATUS "IPO/LTO available")
else ()
    message(STATUS "IPO/LTO not available")
endif ()

include_directories(
        "${CMAKE_CURRENT_SOURCE_DIR}/src/include"
        "${BZIP2_INCLUDE_DIRS}"
        "${CMAKE_CURRENT_SOURCE_DIR}/external/abseil-cpp/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/external/CLI11/include"
        "${Boost_INCLUDE_DIRS}")


if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Debug mode detected: enabling tests")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)
endif ()

add_library(
        libmodle
        src/dna.cpp
        src/extr_barrier.cpp
        src/genome.cpp
        src/lefs.cpp
        src/parsers.cpp)

target_link_libraries(
        libmodle
        absl::base
        absl::flat_hash_map
        absl::strings
        absl::str_format
        absl::time
        BZip2::BZip2
        Boost::iostreams
        Boost::filesystem
        TBB::tbb
        Threads::Threads
)

set_target_properties(libmodle
        PROPERTIES OUTPUT_NAME modle)


target_compile_options(libmodle PUBLIC "$<$<CONFIG:DEBUG>:${MODLE_DEBUG_OPTIONS}>")
target_compile_options(libmodle PUBLIC "$<$<CONFIG:RELEASE>:${MODLE_RELEASE_OPTIONS}>")
target_compile_options(libmodle PUBLIC "$<$<CONFIG:TSAN>:${MODLE_RELEASE_W_TSAN_OPTIONS}>")
target_link_options(libmodle PUBLIC "$<$<CONFIG:TSAN>:-fsanitize=thread>")

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
)

add_custom_command(
        TARGET modle POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_SOURCE_DIR}/utils/plot_contact_matrix.py
        ${CMAKE_CURRENT_BINARY_DIR}/modle_plot_contact_matrix.py
)

target_compile_options(modle PUBLIC "$<$<CONFIG:DEBUG>:${MODLE_DEBUG_OPTIONS}>")
target_compile_options(modle PUBLIC "$<$<CONFIG:RELEASE>:${MODLE_RELEASE_OPTIONS}>")
target_compile_options(modle PUBLIC "$<$<CONFIG:TSAN>:${MODLE_RELEASE_W_TSAN_OPTIONS}>")
target_link_options(modle PUBLIC "$<$<CONFIG:TSAN>:-fsanitize=thread>")

add_subdirectory(utils/compare_contact_matrix)

install(TARGETS libmodle
        LIBRARY
        COMPONENT Libraries
        NAMELINK_COMPONENT Development
        PUBLIC_HEADER
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/modle
        COMPONENT Development)

install(TARGETS modle
        DESTINATION ${CMAKE_INSTALL_BINDIR})

install(PROGRAMS utils/plot_contact_matrix.py
        TYPE BIN
        RENAME "modle_plot_contact_matrix.py")