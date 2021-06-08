add_library(cmatrix INTERFACE)
add_library(Modle::cmatrix ALIAS cmatrix)

target_sources(cmatrix INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/contacts_impl.hpp)

target_include_directories(cmatrix INTERFACE include/)

target_link_libraries(
  cmatrix
  INTERFACE project_warnings
            project_options
            Modle::common
            CONAN_PKG::boost
            CONAN_PKG::cpp-sort)

if(ENABLE_PCH)
  target_precompile_headers(
    cmatrix
    INTERFACE
    <cstdint>
    <utility>
    <vector>
    <boost/dynamic_bitset/dynamic_bitset.hpp>
    [["modle/contacts.hpp"]])
endif()

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()
