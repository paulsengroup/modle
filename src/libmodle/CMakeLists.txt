add_library(libmodle INTERFACE)
add_library(Modle::Modle ALIAS libmodle)

target_sources(libmodle INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/genome_impl.hpp)

target_include_directories(libmodle INTERFACE include/ ${XOSHIRO_INCLUDE_DIR})

target_link_libraries(
  libmodle
  INTERFACE
    project_warnings
    project_options
    Modle::cmatrix
    Modle::io
    $<$<BOOL:${ENABLE_XOSHIRO}>:XoshiroCpp::XoshiroCpp>
    CONAN_PKG::abseil
    CONAN_PKG::boost
    CONAN_PKG::concurrentqueue
    CONAN_PKG::cpp-sort
    # Deal with std::filesystem quirkiness on older compilers
    $<$<AND:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<VERSION_LESS_EQUAL:$<CXX_COMPILER_VERSION>,9.0.0>>:-lc++fs>
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS_EQUAL:$<CXX_COMPILER_VERSION>,9.1.0>>:-lstdc++fs>)

if(ENABLE_PCH AND ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0")
  target_precompile_headers(
    libmodle
    INTERFACE
    <absl/container/btree_set.h>
    <absl/container/flat_hash_map.h>
    <absl/container/flat_hash_set.h>
    <absl/time/clock.h>
    <absl/time/time.h>
    <absl/types/span.h>
    <absl/strings/match.h>
    <absl/strings/str_cat.h>
    <absl/strings/str_join.h>
    <absl/strings/str_split.h>
    <absl/strings/strip.h>
    <boost/asio/thread_pool.hpp>
    <fmt/format.h>)
endif()

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()
