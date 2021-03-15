add_library(cmatrix INTERFACE)
add_library(Modle::cmatrix ALIAS cmatrix)
target_include_directories(cmatrix INTERFACE include/)

target_link_libraries(
  cmatrix
  INTERFACE
    project_warnings
    project_options
    Modle::utils
    $<$<BOOL:${ENABLE_XOSHIRO}>:XoshiroCpp::XoshiroCpp>
    CONAN_PKG::boost
    # Deal with std::filesystem quirkiness on older compilers
    $<$<AND:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<VERSION_LESS_EQUAL:$<CXX_COMPILER_VERSION>,9.0.0>>:-lc++fs>
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS_EQUAL:$<CXX_COMPILER_VERSION>,9.1.0>>:-lstdc++fs>)

if(ENABLE_PCH AND ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0")
  target_precompile_headers(cmatrix INTERFACE [["modle/contacts.hpp"]])
endif()

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()
