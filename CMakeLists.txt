cmake_minimum_required(VERSION 3.11)

project(
  modle
  LANGUAGES CXX
  VERSION 0.0.1
  DESCRIPTION "Software for in-silico modeling of DNA loop extrusion.")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

include(cmake/StandardProjectSettings.cmake)
include(cmake/PreventInSourceBuilds.cmake)
include(CheckCXXCompilerFlag)

# Link this 'library' to set the c++ standard / compile-time options requested
add_library(project_options INTERFACE)

if(ENABLE_CXX20)
  check_cxx_compiler_flag("-std=c++20" CXX20_COMPILER_AVAILABLE)
endif()

if(ENABLE_CXX20 AND CXX20_COMPILER_AVAILABLE)
  target_compile_features(project_options INTERFACE cxx_std_20)
else()
  target_compile_features(project_options INTERFACE cxx_std_17)
endif()

if(ENABLE_TESTING)
  target_compile_definitions(project_options INTERFACE ENABLE_TESTING=1)
endif()

target_compile_definitions(project_options INTERFACE $<$<CONFIG:Release>:BOOST_STACKTRACE_USE_NOOP>)

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

# figure out how to properly link to <filesystem>
include(cmake/Filesystem.cmake)

# detect which standard library is being used
include(cmake/DetectStdlib.cmake)

option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" OFF)
option(ENABLE_TESTING "Enable Test Builds" ON)
option(ENABLE_PCH "Enable Precompiled Headers" OFF)

if(ENABLE_PCH)
  if(${CMAKE_VERSION} VERSION_LESS "3.16.0")
    message(WARNING "CMake 3.16 or newer is required to enable precompiled headers. Ignoring -DENABLE_PCH.")
    set(ENABLE_PCH OFF)
  endif()
endif()

#[[
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
]]

# Set up some extra Conan dependencies based on our needs before loading Conan
set(CONAN_EXTRA_REQUIRES "")
set(CONAN_EXTRA_OPTIONS "")

set(CONAN_EXTRA_REQUIRES
    ${CONAN_EXTRA_REQUIRES}
    abseil/20210324.0
    boost/1.75.0
    cli11/1.9.1
    concurrentqueue/1.0.2
    cpp-sort/1.10.0
    fmt/7.1.3
    hdf5/1.12.0
    readerwriterqueue/1.0.3
    spdlog/1.8.5
    zlib/1.2.11)

set(CONAN_EXTRA_OPTIONS
    ${CONAN_EXTRA_OPTIONS}
    boost:system_no_deprecated=True
    boost:asio_no_deprecated=True
    boost:filesystem_no_deprecated=True
    boost:zlib=True
    boost:bzip2=False)

if(CXX20_COMPILER_AVAILABLE)
  set(CONAN_EXTRA_SETTINGS ${CONAN_EXTRA_SETTINGS} compiler.cppstd=20 compiler.libcxx=libstdc++11)
else()
  set(CONAN_EXTRA_SETTINGS ${CONAN_EXTRA_SETTINGS} compiler.cppstd=17 compiler.libcxx=libstdc++11)
endif()

if(ENABLE_TESTING)
  set(CONAN_EXTRA_REQUIRES ${CONAN_EXTRA_REQUIRES} catch2/2.13.4)
endif()

include(cmake/Conan.cmake)
run_conan()

# import external dependencies
add_subdirectory(external/libBigWig-0.4.6)
if(ENABLE_XOSHIRO)
  target_compile_definitions(project_options INTERFACE USE_XOSHIRO=1)
  add_subdirectory(external/Xoshiro-cpp-1.1)
endif()

if(ENABLE_TESTING)
  enable_testing()
  message(STATUS "Building tests")
  add_subdirectory(test)
endif()

# Reduce level of optimization to improve the code profiling experience
if(OPTIMIZE_FOR_PROFILING)
  if(NOT
     "${CMAKE_BUILD_TYPE}"
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

#[[
if(ENABLE_FUZZING)
  message("Building Fuzz Tests, using fuzzing sanitizer https://www.llvm.org/docs/LibFuzzer.html")
  add_subdirectory(fuzz_test)
endif()
#]]

add_subdirectory(src)
