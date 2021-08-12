cmake_minimum_required(VERSION 3.11)

project(
  modle
  LANGUAGES CXX
  VERSION 0.0.1
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
  target_compile_definitions(project_options INTERFACE ENABLE_TESTING=1)
endif()

target_compile_definitions(project_options INTERFACE $<$<CONFIG:Release>:BOOST_STACKTRACE_USE_NOOP>)

# Define macros to signal library feature availability
if(NOT STD_FILESYSTEM_AVAILABLE)
  target_compile_definitions(project_options INTERFACE MODLE_WITH_BOOST_FILESYSTEM=1)
endif()

if(CHARCONV_INT_AVAILABLE)
  target_compile_definitions(project_options INTERFACE MODLE_CHARCONV_INT_AVAILABLE=1)
endif()

if(CHARCONV_FP_AVAILABLE)
  target_compile_definitions(project_options INTERFACE MODLE_CHARCONV_FP_AVAILABLE=1)
endif()

if(VARIANT_AVAILABLE)
  target_compile_definitions(project_options INTERFACE MODLE_VARIANT_AVAILABLE=1)
endif()

# Tweak spdlog
target_compile_definitions(project_options INTERFACE SPDLOG_CLOCK_COARSE=1)
target_compile_definitions(project_options INTERFACE SPDLOG_DISABLE_DEFAULT_LOGGER=1)
target_compile_definitions(project_options INTERFACE SPDLOG_NO_ATOMIC_LEVELS=1)
target_compile_definitions(project_options INTERFACE SPDLOG_PREVENT_CHILD_FD=1)

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

# Link this 'library' to use the warnings specified in CompilerWarnings.cmake
add_library(project_warnings INTERFACE)

# enable cache system
include(cmake/Cache.cmake)

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

include(cmake/Conan.cmake)
run_conan()

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
