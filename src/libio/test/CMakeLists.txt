add_library(modle_io_test test.cpp)

add_library(Modle::io_test ALIAS modle_io_test)

target_sources(modle_io_test
  INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/bed_parser_test.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cooler_test.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/hdf5_test.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/juicer_contacts_test.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cooler_test.hpp)

target_link_libraries(modle_io_test
  PRIVATE project_warnings
          project_options
  PUBLIC
        Modle::io
        Modle::utils
        CONAN_PKG::abseil
        CONAN_PKG::catch2)
