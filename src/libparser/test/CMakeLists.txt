add_library(modle_parser_test STATIC bed_parser_test.cpp)
add_library(Modle::parser_test ALIAS modle_parser_test)
target_link_libraries(
  modle_parser_test
  PUBLIC CONAN_PKG::catch2
  PRIVATE project_warnings
          project_options
          Modle::utils
          Modle::parser
          CONAN_PKG::abseil)
