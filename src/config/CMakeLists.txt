add_library(modle_config INTERFACE)
add_library(modle_config_cuda INTERFACE)

add_library(Modle::config ALIAS modle_config)
add_library(Modle::config_cuda_only ALIAS modle_config_cuda)

target_sources(modle_config_cuda INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/config_cuda.hpp)
target_sources(modle_config INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/config.hpp)

target_link_libraries(modle_config_cuda INTERFACE project_warnings project_options)
target_link_libraries(modle_config INTERFACE project_warnings project_options modle_config_cuda)

target_include_directories(modle_config_cuda INTERFACE include/)
target_include_directories(modle_config INTERFACE include/)
