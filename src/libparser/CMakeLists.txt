add_library(libparser parsers.cpp)
add_library(Modle::parser ALIAS libparser)

target_include_directories(libparser PUBLIC include/)

target_link_libraries(
  libparser
  PRIVATE project_warnings project_options
  PUBLIC Modle::cmatrix CONAN_PKG::abseil CONAN_PKG::fmt)

set_target_properties(libparser
        PROPERTIES OUTPUT_NAME parser)


if(ENABLE_TESTING)
  add_subdirectory(test)
endif()
