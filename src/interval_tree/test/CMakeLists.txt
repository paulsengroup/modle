find_package(cpp-sort REQUIRED)
find_package(Catch2 REQUIRED)

add_library(modle_interval_tree_test INTERFACE)
add_library(Modle::interval_tree_test ALIAS modle_interval_tree_test)

target_include_directories(modle_interval_tree_test INTERFACE include/)

target_sources(modle_interval_tree_test INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/interval_tree_test.hpp)

target_link_libraries(
  modle_interval_tree_test
  INTERFACE project_warnings
            project_options
            Catch2::Catch2
            cpp-sort::cpp-sort
            Modle::interval_tree)
