add_library(modle_utils INTERFACE)
add_library(Modle::utils ALIAS modle_utils)
target_include_directories(modle_utils INTERFACE include/)
target_link_libraries(modle_utils INTERFACE CONAN_PKG::abseil CONAN_PKG::fmt)
