add_library(libmodle_test INTERFACE)
add_library(Modle::libmodle_test ALIAS libmodle_test)

target_include_directories(libmodle_test INTERFACE include/)

target_sources(libmodle_test INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/simulation_test.hpp)

target_link_libraries(
  libmodle_test
  INTERFACE project_warnings
            project_options
            Modle::io
            Modle::Modle
            CONAN_PKG::abseil
            CONAN_PKG::boost
            CONAN_PKG::catch2)
