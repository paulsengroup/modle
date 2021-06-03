add_subdirectory(internal)
add_subdirectory(cpu)

add_library(libmodle)
add_library(Modle::modle ALIAS libmodle)

target_link_libraries(libmodle PUBLIC libmodle_cpu)
set_target_properties(libmodle PROPERTIES OUTPUT_NAME modle)

if(ENABLE_SHARED)
  install(
    TARGETS libmodle
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES
    DESTINATION include)
endif()
