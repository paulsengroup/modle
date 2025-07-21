# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.25)
cmake_policy(VERSION 3.25...4.0)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

# Not ideal to use this global variable, but necessary to make sure that tooling and projects use the same version
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_C_STANDARD 11)

# strongly encouraged to enable this globally to avoid conflicts between -Wpedantic being enabled and -std=c++20 and
# -std=gnu++20 for example when compiling with PCH enabled
set(CMAKE_CXX_EXTENSIONS OFF)

set(ENABLE_DEVELOPER_MODE
    OFF
    CACHE BOOL "Enable 'developer mode'")

if(BUILD_SHARED_LIBS)
  message(
    FATAL_ERROR
      "Building project with dynamic linking is currently not supported! Please remove option -DBUILD_SHARED_LIBS=ON and re-run cmake."
  )
endif()

if(NOT CMAKE_FIND_PACKAGE_PREFER_CONFIG)
  set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ON)
endif()

set(MODLE_PROJECT_VERSION_MAJOR 1)
set(MODLE_PROJECT_VERSION_MINOR 1)
set(MODLE_PROJECT_VERSION_PATCH 0)
set(MODLE_PROJECT_VERSION_SUFFIX "")

project(
  MoDLE
  LANGUAGES C CXX
  VERSION "${MODLE_PROJECT_VERSION_MAJOR}.${MODLE_PROJECT_VERSION_MINOR}.${MODLE_PROJECT_VERSION_PATCH}"
  HOMEPAGE_URL https://github.com/paulsengroup/modle
  DESCRIPTION "MoDLE: High-performance stochastic modeling of DNA loop extrusion interactions.")

include(FetchContent)
# gersemi: off
FetchContent_Declare(
        _modle_project_options
        URL "${CMAKE_CURRENT_SOURCE_DIR}/external/project_options-v0.36.6.tar.xz"
        URL_HASH SHA256=b4df2a4107847248d4711ee48477cdf93c1a1c7d51a931937759eb9ffba67af8
)
# gersemi: on

FetchContent_MakeAvailable(_modle_project_options)
include("${_modle_project_options_SOURCE_DIR}/Index.cmake")

get_property(BUILDING_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(BUILDING_MULTI_CONFIG)
  if(NOT CMAKE_BUILD_TYPE)
    # Make sure that all supported configuration types have their associated conan packages available. You can reduce
    # this list to only the configuration types you use, but only if one is not forced-set on the command line for VS
    message(TRACE "Setting up multi-config build types")
    set(
      CMAKE_CONFIGURATION_TYPES
      Debug
      Release
      RelWithDebInfo
      CACHE STRING
      "Enabled build types"
      FORCE
    )
  else()
    message(TRACE "User chose a specific build type, so we are using that")
    set(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE} CACHE STRING "Enabled build types" FORCE)
  endif()
endif()

include("${_modle_project_options_SOURCE_DIR}/src/DynamicProjectOptions.cmake")
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/CompilerWarnings.cmake")

# dynamic_project_options sets recommended defaults and provides user and developer modes and full GUI support for
# choosing options at configure time

# for more flexibility, look into project_options() macro

# Any default can be overridden set(<feature_name>_DEFAULT <value>) - set default for both user and developer modes
# set(<feature_name>_DEVELOPER_DEFAULT <value>) - set default for developer mode set(<feature_name>_USER_DEFAULT
# <value>) - set default for user mode

# Initialize project_options variable related to this project This overwrites `project_options` and sets
# `project_warnings` uncomment the options to enable them:

set(ENABLE_CACHE_DEFAULT ON)
set(ENABLE_COMPILE_COMMANDS_SYMLINK_DEFAULT OFF)
set(ENABLE_CONAN_DEFAULT OFF)
set(ENABLE_CLANG_TIDY_DEFAULT OFF)
set(ENABLE_CPPCHECK_DEFAULT OFF)
set(ENABLE_DOXYGEN_DEFAULT OFF)
set(ENABLE_INTERPROCEDURAL_OPTIMIZATION_DEFAULT ON)
set(ENABLE_NATIVE_OPTIMIZATION_DEFAULT OFF)
set(ENABLE_PCH_DEFAULT OFF)

set(ENABLE_SANITIZER_ADDRESS_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_ADDRESS_DEVELOPER_DEFAULT ON)
set(ENABLE_SANITIZER_LEAK_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_LEAK_DEVELOPER_DEFAULT ON)
set(ENABLE_SANITIZER_POINTER_COMPARE_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_POINTER_COMPARE_DEVELOPER_DEFAULT OFF)
set(ENABLE_SANITIZER_POINTER_SUBTRACT_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_POINTER_SUBTRACT_DEVELOPER_DEFAULT ON)
set(ENABLE_SANITIZER_UNDEFINED_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_UNDEFINED_DEVELOPER_DEFAULT ON)

dynamic_project_options(
  PREFIX
  modle
  # CLANG_WARNINGS # Override the defaults for the CLANG warnings GCC_WARNINGS   # Override the defaults for the GCC
  CPPCHECK_OPTIONS
  --enable=performance,portability,style,warning
  --inline-suppr
  # We cannot act on a bug/missing feature of cppcheck
  --suppress=internalAstError
  # if a file does not have an internalAstError, we get an unmatchedSuppression error
  --suppress=unmatchedSuppression
  --suppress=passedByValue
  --inconclusive
  MSVC_WARNINGS
  "${MSVC_WARNINGS}"
  CLANG_WARNINGS
  "${CLANG_WARNINGS}"
  GCC_WARNINGS
  "${GCC_WARNINGS}"
  CUDA_WARNINGS
  "${CUDA_WARNINGS}"
)

# TODO fixme
# This is used so that source files can include "modle/version/version.hpp" and "modle/version/git.hpp" which is
# automatically generated by CMake
include_directories("${CMAKE_CURRENT_BINARY_DIR}/src/version/include")

target_compile_features(modle_project_options INTERFACE "cxx_std_${CMAKE_CXX_STANDARD}")

target_compile_definitions(
  modle_project_options
  INTERFACE
    # Tweak fmt
    FMT_HEADER_ONLY
    FMT_USE_FULL_CACHE_DRAGONBOX
    # Tweak spdlog
    "SPDLOG_ACTIVE_LEVEL=$<IF:$<CONFIG:Debug>,SPDLOG_LEVEL_TRACE,SPDLOG_LEVEL_INFO>"
    SPDLOG_CLOCK_COARSE
    SPDLOG_FMT_EXTERNAL
    SPDLOG_NO_THREAD_ID
    SPDLOG_NO_ATOMIC_LEVELS
    # Tweak xxHash
    XXH_INLINE_ALL
)

include(CMakeDependentOption)
option(MODLE_WITH_BOOST_RANDOM "Use Boost's random library instead of that from the STL (recommended)" ON)
option(MODLE_ENABLE_TESTING "Build MoDLE's unit tests" ON)

cmake_dependent_option(
  MODLE_DOWNLOAD_TEST_DATASET
  "Download datasets required by unit and integration tests"
  ON
  "MODLE_ENABLE_TESTING"
  OFF
)

if(MODLE_WITH_BOOST_RANDOM)
  target_compile_definitions(modle_project_options INTERFACE MODLE_WITH_BOOST_RANDOM=1)
endif()

include(cmake/FetchExternalDeps.cmake)

add_subdirectory(src)

if(MODLE_ENABLE_TESTING)
  enable_testing()
  target_compile_definitions(modle_project_options INTERFACE MODLE_ENABLE_TESTING)
  add_subdirectory(test)
endif()
