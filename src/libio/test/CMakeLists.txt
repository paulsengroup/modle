add_library(
  modle_io_test STATIC
  bed_parser_test.cpp
  cooler_test.cpp
  hdf5_test.cpp
  juicer_contacts_test.cpp
  cooler_test.cpp)
add_library(Modle::io_test ALIAS modle_io_test)
target_link_libraries(
  modle_io_test
  PUBLIC CONAN_PKG::catch2
  PRIVATE project_warnings
          project_options
          Modle::utils
          Modle::io
          CONAN_PKG::abseil)
