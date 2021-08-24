# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE
      RelWithDebInfo
      CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui, ccmake
  set_property(
    CACHE CMAKE_BUILD_TYPE
    PROPERTY STRINGS
             "Debug"
             "Release"
             "MinSizeRel"
             "RelWithDebInfo")
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" uppercase_CMAKE_BUILD_TYPE)

# Generate compile_commands.json to make it easier to work with clang based tools
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" OFF)
option(ENABLE_TESTING "Enable unit tests" ON)
option(ENABLE_PCH "Enable Precompiled Headers" OFF)
option(ENABLE_IPO "Enable Interprocedural Optimization, aka Link Time Optimization (LTO)" OFF)
option(ENABLE_PSO "Enable platform specific optimization" OFF)
option(MODLE_USE_MERSENNE_TWISTER "Use Mersenne Twister instead of Xoshiro as default PRNG engine" OFF)
option(WITH_BOOST_RANDOM "Use Boost's random library instead of that from the STL (recommended)" ON)
option(ENABLE_CXX20 "Compile project under C++20 if this is supported by the compiler" OFF)
option(OPTIMIZE_FOR_PROFILING
       "Compile project in RelWithDebInfo and with less aggressive optimizations to aid profiling" OFF)
option(
  ENABLE_ASSERTIONS
  "Enable assertions and various other runtime checks (this is done regardless of the type passed to CMAKE_BUILD_TYPE)"
  OFF)

if(ENABLE_PCH)
  if(${CMAKE_VERSION} VERSION_LESS "3.16.0")
    message(WARNING "CMake 3.16 or newer is required to enable precompiled headers. Ignoring -DENABLE_PCH.")
    set(ENABLE_PCH OFF)
  endif()
endif()

if(ENABLE_IPO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT result OUTPUT output)
  if(result)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
  else()
    message(SEND_ERROR "IPO is not supported: ${output}")
  endif()
endif()
