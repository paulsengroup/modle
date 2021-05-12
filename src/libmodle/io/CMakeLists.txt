add_library(libmodle_io INTERFACE)
add_library(Modle::modle_io_internal ALIAS libmodle_io)

target_sources(libmodle_io INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/libmodle_io.hpp ${CMAKE_CURRENT_SOURCE_DIR}/libmodle_io_impl.hpp)
target_link_libraries(
  libmodle_io
  INTERFACE project_options
            project_warnings
            Modle::io
            CONAN_PKG::abseil)

target_include_directories(libmodle_io INTERFACE include/)
