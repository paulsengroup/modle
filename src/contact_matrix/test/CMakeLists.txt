add_library(modle_cmatrix_test INTERFACE)
add_library(Modle::cmatrix_test ALIAS modle_cmatrix_test)

target_include_directories(modle_cmatrix_test INTERFACE include/)

target_sources(modle_cmatrix_test INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/contacts_test.hpp)

target_link_libraries(
  modle_cmatrix_test
  INTERFACE project_warnings
            project_options
            Modle::cmatrix
            Modle::utils
            CONAN_PKG::abseil
            CONAN_PKG::boost
            CONAN_PKG::catch2)
