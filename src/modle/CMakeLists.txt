find_package(Threads REQUIRED)

add_executable(modle main.cpp)

target_link_libraries(
  modle
  PRIVATE project_warnings
          project_options
          Modle::Modle
          Modle::io
          CONAN_PKG::abseil
          CONAN_PKG::boost
          CONAN_PKG::cli11)

install(
  TARGETS modle
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)
