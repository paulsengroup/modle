# find_package(Threads REQUIRED)

add_library(libmodle INTERFACE)
add_library(Modle::Modle ALIAS libmodle)

target_sources(
  libmodle
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/dna_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/extr_barrier_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/genome_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/lefs_impl.hpp)

target_include_directories(libmodle INTERFACE include/)

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
