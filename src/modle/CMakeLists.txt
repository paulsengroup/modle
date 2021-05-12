add_executable(modle main.cpp)

target_link_libraries(
  modle
  PRIVATE project_warnings
          project_options
          Modle::modle
          Modle::io
          CONAN_PKG::abseil
          CONAN_PKG::boost
          CONAN_PKG::cli11
          # Deal with std::filesystem quirkiness on older compilers
          $<$<AND:$<BOOL:${USING_LIBCXX}>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lc++fs>
          $<$<AND:$<BOOL:${USING_STDLIBCXX}>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lstdc++fs>)

install(
  TARGETS modle
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)
