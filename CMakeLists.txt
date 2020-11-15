cmake_minimum_required(VERSION 3.11)

# Set the project name to your project name, my project isn't very descriptive
project(
  modle
  LANGUAGES CXX
  VERSION 0.0.1
  DESCRIPTION "Software for in-silico modeling of DNA loop extrusion.")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

include(cmake/StandardProjectSettings.cmake)
include(cmake/PreventInSourceBuilds.cmake)

# Link this 'library' to set the c++ standard / compile-time options requested
add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_17)

if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
  option(ENABLE_BUILD_WITH_TIME_TRACE "Enable -ftime-trace to generate time tracing .json files on clang" OFF)
  if(ENABLE_BUILD_WITH_TIME_TRACE)
    add_compile_definitions(project_options INTERFACE -ftime-trace)
  endif()
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

option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" OFF)
option(ENABLE_TESTING "Enable Test Builds" ON)

# Very basic PCH example
option(ENABLE_PCH "Enable Precompiled Headers" OFF)
if(ENABLE_PCH)
  # This sets a global PCH parameter, each project will build its own PCH, which is a good idea if any #define's change
  #
  # consider breaking this out per project as necessary
  target_precompile_headers(
    project_options
    INTERFACE
    <vector>
    <string>
    <utility>)
endif()

# Make sure compiler support execution policies (C++17)
include(CheckIncludeFileCXX)
set(CMAKE_REQUIRED_INCLUDES "${CMAKE_INCLUDE_PATH}")
check_include_file_cxx("execution" HAVE_EXECUTION)
if(NOT HAVE_EXECUTION)
  message(
    FATAL_ERROR
      "Your compiler C++ library does not seem to support execution policies (C++17, currently only implemented in libstdc++ >=9)."
  )
endif()

# Set up some extra Conan dependencies based on our needs before loading Conan
set(CONAN_EXTRA_REQUIRES "")
set(CONAN_EXTRA_OPTIONS "")

set(CONAN_EXTRA_REQUIRES
    ${CONAN_EXTRA_REQUIRES}
    tbb/2020.2
    zlib/1.2.11
    bzip2/1.0.8
    boost/1.74.0
    abseil/20200923.2
    fmt/7.1.2
    spdlog/1.8.1
    cli11/1.9.1)
set(CONAN_EXTRA_OPTIONS
    ${CONAN_EXTRA_OPTIONS}
    bzip2:build_executable=False
    zlib:minizip=False
    boost:system_no_deprecated=True
    boost:asio_no_deprecated=True
    boost:filesystem_no_deprecated=True
    boost:zlib=True
    boost:bzip2=True)
set(CONAN_EXTRA_SETTINGS ${CONAN_EXTRA_SETTINGS} compiler.cppstd=17)

if(ENABLE_TESTING)
  set(CONAN_EXTRA_REQUIRES ${CONAN_EXTRA_REQUIRES} catch2/2.11.1)
endif()

include(cmake/Conan.cmake)
run_conan()

if(ENABLE_TESTING)
  enable_testing()
  message("Building Tests")
  add_subdirectory(test)
endif()

#[[
if(ENABLE_FUZZING)
  message("Building Fuzz Tests, using fuzzing sanitizer https://www.llvm.org/docs/LibFuzzer.html")
  add_subdirectory(fuzz_test)
endif()
#]]

add_subdirectory(src)
