add_library(modle_common INTERFACE)
add_library(modle_common_cuda INTERFACE)

add_library(Modle::common ALIAS modle_common)

target_sources(modle_common INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/config.hpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common.hpp)

target_link_libraries(
  modle_common
  INTERFACE project_options
            project_warnings
            Modle::utils
            $<$<BOOL:${ENABLE_XOSHIRO}>:XoshiroCpp::XoshiroCpp>)

target_include_directories(modle_common INTERFACE include/)
