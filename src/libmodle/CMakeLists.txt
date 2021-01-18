add_library(libmodle INTERFACE)
add_library(Modle::Modle ALIAS libmodle)

if (ENABLE_XOSHIRO)
  include(ExternalProject)
  ExternalProject_Add(
        Xoshiro-cpp
        GIT_REPOSITORY https://github.com/Reputeless/Xoshiro-cpp
        GIT_TAG "v1.1"
        CONFIGURE_COMMAND cp "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists_XoshiroCpp.txt" "<SOURCE_DIR>/CMakeLists.txt"
        COMMAND ${CMAKE_COMMAND} <SOURCE_DIR> -Wno-dev
        INSTALL_COMMAND ""
)
  add_dependencies(libmodle Xoshiro-cpp)
  ExternalProject_Get_Property(Xoshiro-cpp SOURCE_DIR)
  set(XOSHIRO_INCLUDE_DIR ${SOURCE_DIR})
  unset(SOURCE_DIR)
endif()


target_sources(
  libmodle
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/dna_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/extr_barrier_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/genome_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/lefs_impl.hpp)

target_include_directories(libmodle
        INTERFACE
        include/
        ${XOSHIRO_INCLUDE_DIR})

target_link_libraries(
  libmodle
  INTERFACE
    project_warnings
    project_options
    Modle::cmatrix
    Modle::io
    CONAN_PKG::abseil
    CONAN_PKG::boost
    # Deal with std::filesystem quirkiness on older compilers
    $<$<AND:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<VERSION_LESS_EQUAL:$<CXX_COMPILER_VERSION>,9.0.0>>:-lc++fs>
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS_EQUAL:$<CXX_COMPILER_VERSION>,9.1.0>>:-lstdc++fs>)

if(ENABLE_PCH AND ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0")
  target_precompile_headers(
    libmodle
    INTERFACE
    [["modle/common.hpp"]]
    [["modle/config.hpp"]]
    [["modle/dna.hpp"]]
    [["modle/extr_barrier.hpp"]]
    [["modle/genome.hpp"]]
    [["modle/lefs.hpp"]])
endif()
