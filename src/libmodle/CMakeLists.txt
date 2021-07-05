add_subdirectory(internal)
add_subdirectory(cpu)

add_library(libmodle)
add_library(Modle::libmodle ALIAS libmodle)

target_link_libraries(libmodle PUBLIC libmodle_cpu)
set_target_properties(libmodle PROPERTIES OUTPUT_NAME modle)

install(
  TARGETS libmodle
  ARCHIVE
  INCLUDES
  LIBRARY)
