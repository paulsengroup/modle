add_library(modle_correlation INTERFACE)
add_library(Modle::correlation ALIAS modle_correlation)
target_include_directories(modle_correlation INTERFACE include/ CONAN_PKG::range-v3)
target_link_libraries(modle_correlation INTERFACE project_warnings project_options CONAN_PKG::range-v3)

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()
