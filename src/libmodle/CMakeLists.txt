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
  PRIVATE project_warnings
          project_options
          Modle::parser
          CONAN_PKG::boost
          CONAN_PKG::tbb
  PUBLIC Threads::Threads Modle::cmatrix CONAN_PKG::abseil)

set_target_properties(libmodle PROPERTIES OUTPUT_NAME modle)

install(
  TARGETS libmodle
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)
