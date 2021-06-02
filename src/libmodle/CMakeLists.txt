add_subdirectory(internal)
add_subdirectory(cpu)

add_library(libmodle INTERFACE)
add_library(Modle::modle ALIAS libmodle)

target_link_libraries(libmodle INTERFACE libmodle_cpu)
