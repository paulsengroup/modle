add_library(modle_correlation INTERFACE)
add_library(Modle::correlation ALIAS modle_correlation)
target_include_directories(modle_correlation INTERFACE include/)
target_link_libraries(modle_correlation INTERFACE project_warnings project_options CONAN_PKG::boost)

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()

#[[
install(
  TARGETS modle_correlation
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)
]]
