find_package(absl REQUIRED)
find_package(Catch2 REQUIRED)
find_package(HDF5 REQUIRED)

add_library(modle_io_test INTERFACE)

add_library(Modle::libmodle_io_test ALIAS modle_io_test)

target_include_directories(modle_io_test INTERFACE include/)

target_sources(
  modle_io_test
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/bed_parser_test.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/cooler_test.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/hdf5_test.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/cooler_test.hpp)

target_link_libraries(
  modle_io_test
  INTERFACE project_warnings
            project_options
            absl::span
            absl::strings
            Catch2::Catch2
            HDF5::HDF5
            Modle::common
            Modle::libmodle_io)
