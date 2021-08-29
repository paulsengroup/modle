find_package(absl REQUIRED)
find_package(Boost REQUIRED)

add_library(modle_correlation INTERFACE)
add_library(Modle::correlation ALIAS modle_correlation)
target_include_directories(modle_correlation INTERFACE include/)
target_link_libraries(modle_correlation INTERFACE project_warnings project_options absl::span)

if(ENABLE_PCH)
  target_precompile_headers(
    modle_correlation
    INTERFACE
    <algorithm>
    <cstddef>
    <cstdint>
    <type_traits>
    <utility>
    <vector>
    <absl/types/span.h>
    <boost/dynamic_bitset/dynamic_bitset.hpp>
    [["modle/correlation.hpp"]])
endif()
