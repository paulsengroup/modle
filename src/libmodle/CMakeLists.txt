add_subdirectory(io)
add_subdirectory(cpu)
if(ENABLE_CUDA)
  add_subdirectory(cu)
endif()

add_library(libmodle INTERFACE)
add_library(Modle::modle ALIAS libmodle)

target_link_libraries(libmodle INTERFACE libmodle_cpu $<$<BOOL:${ENABLE_CUDA}>:libmodle_cuda>)
