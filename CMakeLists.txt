# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.16)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

# Not ideal to use this global variable, but necessary to make sure that tooling and projects use the same version
set(CMAKE_CXX_STANDARD 17)

# strongly encouraged to enable this globally to avoid conflicts between -Wpedantic being enabled and -std=c++20 and
# -std=gnu++20 for example when compiling with PCH enabled
set(CMAKE_CXX_EXTENSIONS OFF)

set(ENABLE_DEVELOPER_MODE
    OFF
    CACHE BOOL "Enable 'developer mode'")

include(FetchContent)
# cmake-format: off
FetchContent_Declare(
        _project_options
        URL https://github.com/cpp-best-practices/project_options/archive/refs/tags/v0.17.0.tar.gz
        URL_HASH SHA512=8478eb239a651ff773ee7b5424d33b523a8c8dcbb3c530e906c244d1c34fddd515a58b87d1435afa2bb483e165b937b9bf304fd6c9520554f085d8ea0417a07f
)
# cmake-format: on

FetchContent_MakeAvailable(_project_options)
include(${_project_options_SOURCE_DIR}/Index.cmake)

project(
  MoDLE
  LANGUAGES CXX
  VERSION 0.0.1
  HOMEPAGE_URL https://github.com/robomics/modle
  DESCRIPTION "MoDLE: stochastic modeling of DNA loop extrusion.")

get_property(BUILDING_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(BUILDING_MULTI_CONFIG)
  if(NOT CMAKE_BUILD_TYPE)
    # Make sure that all supported configuration types have their associated conan packages available. You can reduce
    # this list to only the configuration types you use, but only if one is not forced-set on the command line for VS
    message(TRACE "Setting up multi-config build types")
    set(CMAKE_CONFIGURATION_TYPES
        Debug Release RelWithDebInfo
        CACHE STRING "Enabled build types" FORCE)
  else()
    message(TRACE "User chose a specific build type, so we are using that")
    set(CMAKE_CONFIGURATION_TYPES
        ${CMAKE_BUILD_TYPE}
        CACHE STRING "Enabled build types" FORCE)
  endif()
endif()

include(${_project_options_SOURCE_DIR}/src/DynamicProjectOptions.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CompilerWarnings.cmake)

# dynamic_project_options sets recommended defaults and provides user and developer modes and full GUI support for
# choosing options at configure time

# for more flexibility, look into project_options() macro

# Any default can be overridden set(<feature_name>_DEFAULT <value>) - set default for both user and developer modes
# set(<feature_name>_DEVELOPER_DEFAULT <value>) - set default for developer mode set(<feature_name>_USER_DEFAULT
# <value>) - set default for user mode

# Initialize project_options variable related to this project This overwrites `project_options` and sets
# `project_warnings` uncomment the options to enable them:
dynamic_project_options(
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
  ${MSVC_WARNINGS}
  CLANG_WARNINGS
  ${CLANG_WARNINGS}
  GCC_WARNINGS
  ${GCC_WARNINGS}
  CUDA_WARNINGS
  ${CUDA_WARNINGS})

set(ENABLE_CACHE_DEVELOPER ON)
set(ENABLE_CACHE_USER OFF)
set(ENABLE_CLANG_TIDY_DEFAULT OFF)
set(ENABLE_CONAN_DEFAULT ON)
set(ENABLE_CPPCHECK_DEFAULT OFF)
set(ENABLE_DOXYGEN_USER OFF)
set(ENABLE_DOXYGEN_DEVELOPER ON)
set(ENABLE_IPO_DEFAULT OFF)
set(ENABLE_PCH_DEFAULT OFF)

option(MODLE_USE_MERSENNE_TWISTER "Use Mersenne Twister instead of Xoshiro as default PRNG engine" OFF)
option(WITH_BOOST_RANDOM "Use Boost's random library instead of that from the STL (recommended)" ON)
option(OPTIMIZE_FOR_PROFILING
       "Compile project in RelWithDebInfo and with less aggressive optimizations to aid profiling" OFF)
option(
  ENABLE_ASSERTIONS
  "Enable assertions and various other runtime checks (this is done regardless of the type passed to CMAKE_BUILD_TYPE)"
  OFF)

target_compile_features(project_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

include(cmake/TestStdlibFeatures.cmake)

# Define macros based on library feature availability
target_compile_definitions(
  project_options
  INTERFACE $<$<BOOL:${CHARCONV_INT_AVAILABLE}>:MODLE_CHARCONV_INT_AVAILABLE>
            $<$<BOOL:${CHARCONV_FP_AVAILABLE}>:MODLE_CHARCONV_FP_AVAILABLE>
            $<$<BOOL:${CHARCONV_VARIANT_AVAILABLE}>:MODLE_CHARCONV_VARIANT_AVAILABLE>)

# cmake-format: off
# Define preproc. variables for proj. versioning
target_compile_definitions(
  project_options INTERFACE MODLE_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
                            MODLE_VERSION_MINOR=${PROJECT_VERSION_MINOR}
                            MODLE_VERSION_PATCH=${PROJECT_VERSION_PATCH})
# cmake-format: on

# Reduce level of optimization to improve the code profiling experience
if(OPTIMIZE_FOR_PROFILING)
  if(NOT
     "${uppercase_CMAKE_BUILD_TYPE}"
     STREQUAL
     "RelWithDebInfo")
    message(
      WARNING
        "Switching build type from ${CMAKE_BUILD_TYPE} to RelWithDebInfo because OPTIMIZE_FOR_PROFILING was set to ON by the user"
    )
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
  endif()

  target_compile_options(project_options INTERFACE -O1 -fno-omit-frame-pointer)
endif()

# Source:
# https://github.com/llvm/llvm-project/blob/656ebd519e3fd52050e1c8abeacafcf94d1fa260/llvm/cmake/modules/HandleLLVMOptions.cmake#L54
if(ENABLE_ASSERTIONS)
  if(NOT MSVC)
    add_definitions(-D_DEBUG)
  endif()
  # On non-Debug builds cmake automatically defines NDEBUG, so we explicitly undefine it:
  if(NOT
     uppercase_CMAKE_BUILD_TYPE
     STREQUAL
     "DEBUG")
    # NOTE: use `add_compile_options` rather than `add_definitions` since `add_definitions` does not support generator
    # expressions.
    target_compile_options(project_options INTERFACE $<$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>:-UNDEBUG>)

    # Also remove /D NDEBUG to avoid MSVC warnings about conflicting defines.
    foreach(
      flags_var_to_scrub
      CMAKE_CXX_FLAGS_RELEASE
      CMAKE_CXX_FLAGS_RELWITHDEBINFO
      CMAKE_CXX_FLAGS_MINSIZEREL
      CMAKE_C_FLAGS_RELEASE
      CMAKE_C_FLAGS_RELWITHDEBINFO
      CMAKE_C_FLAGS_MINSIZEREL)
      string(
        REGEX
        REPLACE "(^| )[/-]D *NDEBUG($| )"
                " "
                "${flags_var_to_scrub}"
                "${${flags_var_to_scrub}}")
    endforeach()
  endif()
endif()

if(WITH_BOOST_RANDOM)
  target_compile_definitions(project_options INTERFACE MODLE_WITH_BOOST_RANDOM=1)
endif()

message(STATUS "CPU arch detected by CMake: ${CMAKE_SYSTEM_PROCESSOR}")
if(ENABLE_PSO AND ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64" OR "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL AMD64))
  if(MSVC)
    target_compile_options(project_options INTERFACE /arch:native)
  else()
    target_compile_options(project_options INTERFACE -march=native)
  endif()
  if(ENABLE_PCH)
    message(WARNING "Disabling precompiled headers as they tend to cause problems when combined with -DENABLE_PSO=ON")
    set(ENABLE_PCH OFF)
  endif()
endif()

# Prepare external dependencies
include(cmake/FetchExternalDeps.cmake)

# Tweak spdlog
target_compile_definitions(
  project_options
  INTERFACE SPDLOG_CLOCK_COARSE
  INTERFACE SPDLOG_DISABLE_DEFAULT_LOGGER
  INTERFACE SPDLOG_NO_ATOMIC_LEVELS
  INTERFACE SPDLOG_PREVENT_CHILD_FD)

# Tweak fmt
target_compile_definitions(project_options INTERFACE FMT_ENFORCE_COMPILE_STRING)

# Tweak boost
find_package(Boost CONFIG REQUIRED)

if(Boost_VERSION_STRING VERSION_GREATER_EQUAL 1.77.0)
  target_compile_definitions(project_options INTERFACE BOOST_FILESYSTEM_VERSION=4)
elseif(Boost_VERSION_STRING VERSION_GREATER 1.44.0)
  target_compile_definitions(project_options INTERFACE BOOST_FILESYSTEM_VERSION=3)
endif()

if(ENABLE_TESTING)
  enable_testing()
  message("Building unit tests.")
  target_compile_definitions(project_options INTERFACE ENABLE_TESTING)
  add_subdirectory(test)
endif()

add_subdirectory(src)

include(cmake/Packaging.cmake)
