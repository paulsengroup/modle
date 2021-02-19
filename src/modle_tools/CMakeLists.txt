# find_package(Threads REQUIRED)

add_executable(modle_tools main.cpp)

target_sources(
  modle_tools
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cli.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/eval_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/tools_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/stats_impl.hpp)

target_link_libraries(
  modle_tools
  PRIVATE project_warnings
          project_options
          Modle::Modle
          Modle::correlation
          Modle::io
          CONAN_PKG::abseil
          CONAN_PKG::boost
          CONAN_PKG::cli11)

target_include_directories(modle_tools PRIVATE include)

if(ENABLE_PCH AND ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16.0")
  target_precompile_headers(
    modle_tools
    PRIVATE
    [["modle_tools/config.hpp"]]
    [["modle_tools/eval.hpp"]]
    [["modle_tools/tools.hpp"]]
    [["modle_tools/stats.hpp"]]
    [["modle_tools/utils.hpp"]])
endif()

install(
  TARGETS modle_tools
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)
