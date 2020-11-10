cmake_minimum_required(VERSION 3.10)
project(modle_tools
        HOMEPAGE_URL "https://github.com/robomics/modle"
        LANGUAGES CXX
        VERSION 0.0.1
        DESCRIPTION "Helper tool to process the contact matrix produced by modle.")

cmake_policy(VERSION 3.10)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to \"Release\" as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif ()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CheckIPOSupported)
find_package(TBB REQUIRED)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

check_ipo_supported(RESULT IPO_AVAILABLE)
if (IPO_AVAILABLE)
    message(STATUS "IPO/LTO available")
else ()
    message(STATUS "IPO/LTO not available")
endif ()

include_directories(modle)

add_executable(modle_tools main.cpp)

target_link_libraries(modle_tools
        libmodle
        CLI11::CLI11
        absl::base
        absl::time
        absl::strings)

set(DEBUG_OPTIONS "-Wall;-Wextra;-ggdb;-Og;-g3")
set(RELEASE_OPTIONS "-Wall;-Wextra;-O3")
target_compile_options(modle_tools PUBLIC "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")
target_compile_options(modle_tools PUBLIC "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")

install(TARGETS modle_tools
        DESTINATION ${CMAKE_INSTALL_BINDIR})