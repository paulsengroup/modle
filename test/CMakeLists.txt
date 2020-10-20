cmake_minimum_required(VERSION 3.11.4)
project(modle_test)

cmake_policy(VERSION 3.11.4)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(FetchContent)

FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG 477998eefac3d55831e7f3974671ba22b8739b93
        GIT_SHALLOW OFF # GIT_SHALLOW does not allow using commit hashes
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/googletest
)

set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(BUILD_GTEST ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

add_executable(modle_test)

target_sources(
        modle_test
        PRIVATE units/contacts_test.cpp
)

target_include_directories(
        modle_test
        SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/googletest/include/
)

target_link_libraries(
        modle_test
        libmodle
        gtest_main
)

set(MODLE_DEBUG_OPTIONS "-Wall;-Wextra;-ggdb;-Og;-g3")
set(MODLE_RELEASE_OPTIONS "-Wall;-Wextra;-O3")

target_compile_options(modle_test PUBLIC "$<$<CONFIG:DEBUG>:${MODLE_DEBUG_OPTIONS}>")
target_compile_options(modle_test PUBLIC "$<$<CONFIG:RELEASE>:${MODLE_RELEASE_OPTIONS}>")
