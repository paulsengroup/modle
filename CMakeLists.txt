# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.11)

project(
  modle
  LANGUAGES CXX
  VERSION 0.4.1
  DESCRIPTION "MoDLE: stochastic modeling of DNA loop extrusion.")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

include(cmake/StandardProjectSettings.cmake)
include(cmake/PreventInSourceBuilds.cmake)
include(cmake/TestStdlibFeatures.cmake)
include(CheckCXXCompilerFlag)

# Link this 'library' to set the c++ standard / compile-time options requested
add_library(project_options INTERFACE)

if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
  option(ENABLE_BUILD_WITH_TIME_TRACE "Enable -ftime-trace to generate time tracing .json files on clang" OFF)
  if(ENABLE_BUILD_WITH_TIME_TRACE)
    add_compile_definitions(project_options INTERFACE -ftime-trace)
  endif()
endif()

check_cxx_compiler_flag("-std=c++20" MODLE_CXX20_COMPILER_AVAILABLE)
if(MODLE_CXX20_COMPILER_AVAILABLE)
  set(MODLE_CXX20_COMPILER_AVAILABLE ON)
else()
  set(MODLE_CXX20_COMPILER_AVAILABLE OFF)
endif()

if(MODLE_ENABLE_CXX20 AND MODLE_CXX20_COMPILER_AVAILABLE)
  set(CMAKE_CXX_STANDARD 20)
  target_compile_features(project_options INTERFACE cxx_std_20)
else()
  set(CMAKE_CXX_STANDARD 17)
  target_compile_features(project_options INTERFACE cxx_std_17)
endif()

if(ENABLE_TESTING)
  target_compile_definitions(project_options INTERFACE ENABLE_TESTING)
endif()

# Define macros to signal library feature availability
target_compile_definitions(
  project_options
  INTERFACE $<$<BOOL:${CHARCONV_INT_AVAILABLE}>:MODLE_CHARCONV_INT_AVAILABLE>
            $<$<BOOL:${CHARCONV_FP_AVAILABLE}>:MODLE_CHARCONV_FP_AVAILABLE>
            $<$<BOOL:${CHARCONV_VARIANT_AVAILABLE}>:MODLE_CHARCONV_VARIANT_AVAILABLE>)

# Tweak spdlog
target_compile_definitions(
  project_options
  INTERFACE SPDLOG_CLOCK_COARSE
  INTERFACE SPDLOG_DISABLE_DEFAULT_LOGGER
  INTERFACE SPDLOG_NO_ATOMIC_LEVELS
  INTERFACE SPDLOG_PREVENT_CHILD_FD)

# Tweak fmt
target_compile_definitions(project_options INTERFACE FMT_ENFORCE_COMPILE_STRING)

# Define preproc. variables for proj. versioning
target_compile_definitions(
  project_options INTERFACE MODLE_VERSION_MAJOR=${PROJECT_VERSION_MAJOR} MODLE_VERSION_MINOR=${PROJECT_VERSION_MINOR}
                            MODLE_VERSION_PATCH=${PROJECT_VERSION_PATCH})

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

  target_compile_options(project_options INTERFACE -O1)
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

# Link this 'library' to use the warnings specified in CompilerWarnings.cmake
add_library(project_warnings INTERFACE)

# enable cache system
include(cmake/Cache.cmake)

# Add linker configuration
include(cmake/Linker.cmake)
configure_linker(project_options)

# standard compiler warnings
include(cmake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

# sanitizer options if supported by compiler
include(cmake/Sanitizers.cmake)
enable_sanitizers(project_options)

# enable doxygen
include(cmake/Doxygen.cmake)
enable_doxygen()

# allow for static analysis options
include(cmake/StaticAnalyzers.cmake)

option(ENABLE_CONAN "Use Conan for dependency management" ON)
if(ENABLE_CONAN)
  include(cmake/Conan.cmake)
  run_conan()
endif()

# Prepare external dependencies
include(cmake/FetchExternalDeps.cmake)

if(ENABLE_TESTING)
  enable_testing()
  message(STATUS "Building tests")
  add_subdirectory(test)
endif()

#[[
if(ENABLE_FUZZING)
  message("Building Fuzz Tests, using fuzzing sanitizer https://www.llvm.org/docs/LibFuzzer.html")
  add_subdirectory(fuzz_test)
endif()
#]]

add_subdirectory(src)

include(cmake/Packaging.cmake)
