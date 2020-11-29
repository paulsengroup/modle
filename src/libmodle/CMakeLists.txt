find_package(Threads REQUIRED)

add_library(
  libmodle
  dna.cpp
  extr_barrier.cpp
  genome.cpp
  lefs.cpp)
add_library(Modle::Modle ALIAS libmodle)

target_include_directories(libmodle PUBLIC include/)

target_link_libraries(
  libmodle
  PRIVATE
    project_warnings
    project_options
    Modle::io
    CONAN_PKG::boost
    # Deal with std::filesystem quirkiness on older compilers
    $<$<AND:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<VERSION_LESS_EQUAL:$<CXX_COMPILER_VERSION>,9.0.0>>:-lc++fs>
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS_EQUAL:$<CXX_COMPILER_VERSION>,9.1.0>>:-lstdc++fs>
  PUBLIC Threads::Threads Modle::cmatrix CONAN_PKG::abseil)

set_target_properties(libmodle PROPERTIES OUTPUT_NAME modle)

#[[
install(
  TARGETS libmodle
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)
]]
