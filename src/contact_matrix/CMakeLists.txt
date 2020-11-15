add_library(cmatrix INTERFACE)
add_library(Modle::cmatrix ALIAS cmatrix)
target_include_directories(cmatrix INTERFACE include/)
target_link_libraries(cmatrix INTERFACE Modle::utils project_warnings project_options)

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()
