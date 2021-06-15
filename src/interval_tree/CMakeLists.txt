add_library(interval_tree INTERFACE)
add_library(Modle::interval_tree ALIAS interval_tree)

target_sources(interval_tree INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/interval_tree_impl.hpp)

target_include_directories(interval_tree INTERFACE include/)

target_link_libraries(
  interval_tree
  INTERFACE project_warnings
            project_options
            CONAN_PKG::cpp-sort)

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()

