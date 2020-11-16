add_library(modle_cmatrix_test STATIC contacts_test.cpp)
add_library(Modle::cmatrix_test ALIAS modle_cmatrix_test)
target_link_libraries(
  modle_cmatrix_test
  PUBLIC CONAN_PKG::catch2
  PRIVATE project_warnings
          project_options
          Modle::cmatrix
          Modle::utils
          CONAN_PKG::abseil
          CONAN_PKG::boost)
