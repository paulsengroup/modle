find_package(Threads REQUIRED)

add_library(
  libmodle
  dna.cpp
  extr_barrier.cpp
  genome.cpp
  lefs.cpp
  parsers.cpp)

target_include_directories(libmodle PUBLIC include/)

target_link_libraries(
  libmodle
  PRIVATE project_warnings
          project_options
          CONAN_PKG::boost
          CONAN_PKG::tbb
  PUBLIC Threads::Threads CONAN_PKG::abseil)

set_target_properties(libmodle PROPERTIES OUTPUT_NAME modle)

install(
  TARGETS libmodle
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)
