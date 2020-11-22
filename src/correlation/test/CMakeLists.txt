find_package(Threads REQUIRED)

add_library(modle_correlation_test STATIC correlation_test.cpp correlation_utils_test.cpp common.hpp)
add_library(Modle::correlation_test ALIAS modle_correlation_test)
target_link_libraries(
  modle_correlation_test
  PUBLIC CONAN_PKG::catch2
  PRIVATE project_warnings
          project_options
          Modle::correlation
          Modle::utils
          CONAN_PKG::abseil
          CONAN_PKG::boost
          Threads::Threads)
