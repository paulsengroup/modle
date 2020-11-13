find_package(Threads REQUIRED)

add_executable(modle main.cpp)

target_link_libraries(
  modle
  PRIVATE project_warnings
          project_options
          libmodle
          CONAN_PKG::abseil
          CONAN_PKG::boost
          CONAN_PKG::cli11)
