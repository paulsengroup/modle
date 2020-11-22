find_package(Threads REQUIRED)

add_executable(
  modle_tools
  main.cpp
  convert.cpp
  eval.cpp
  tools.cpp)

target_link_libraries(
  modle_tools
  PRIVATE project_warnings
          project_options
          Modle::Modle
          Modle::correlation
          Modle::io
          CONAN_PKG::abseil
          CONAN_PKG::boost
          CONAN_PKG::cli11
          CONAN_PKG::range-v3)

target_include_directories(modle_tools PRIVATE include)

install(
  TARGETS modle_tools
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)
