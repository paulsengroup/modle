add_library(modle_utils INTERFACE)
add_library(Modle::utils ALIAS modle_utils)
target_include_directories(modle_utils INTERFACE include/)
target_link_libraries(modle_utils INTERFACE CONAN_PKG::abseil CONAN_PKG::boost CONAN_PKG::fmt)

if(ENABLE_PCH)
  target_precompile_headers(
    modle_utils
    INTERFACE
    <cstddef>
    <cstdint>
    <string>
    <string_view>
    <type_traits>
    <utility>
    <vector>
    [["modle/suppress_compiler_warnings.hpp"]]
    [["modle/utils.hpp"]])
endif()
