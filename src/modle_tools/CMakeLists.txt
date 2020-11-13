find_package(Threads REQUIRED)

add_executable(modle_tools main.cpp convert.cpp eval.cpp)

target_link_libraries(
  modle_tools
  PRIVATE libmodle
          project_warnings
          project_options
          CONAN_PKG::abseil
          CONAN_PKG::cli11
          CONAN_PKG::boost)

target_include_directories(modle_tools PRIVATE include)
