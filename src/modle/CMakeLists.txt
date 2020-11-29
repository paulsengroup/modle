find_package(Threads REQUIRED)

add_executable(modle main.cpp)

target_link_libraries(
  modle
  PRIVATE
    project_warnings
    project_options
    Modle::Modle
    Modle::io
    CONAN_PKG::abseil
    CONAN_PKG::boost
    CONAN_PKG::cli11
    # Deal with std::filesystem quirkiness on older compilers
    $<$<AND:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<VERSION_LESS_EQUAL:$<CXX_COMPILER_VERSION>,9.0.0>>:-lc++fs>
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS_EQUAL:$<CXX_COMPILER_VERSION>,9.1.0>>:-lstdc++fs>)

install(
  TARGETS modle
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)
