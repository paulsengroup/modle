add_library(cmatrix INTERFACE)
add_library(Modle::cmatrix ALIAS cmatrix)
target_include_directories(cmatrix INTERFACE include/)
target_link_libraries(
  cmatrix
  INTERFACE project_warnings
            project_options
            Modle::utils
            CONAN_PKG::boost)

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()

#[[
install(
  TARGETS cmatrix
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)
]]
