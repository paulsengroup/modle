option(MODLE_BUILD_UTILS "Enable compilation of Modle's utilities (such as modle_tools)" ON)

if(MODLE_BUILD_UTILS)
  message("Building ModLE utils")
  add_subdirectory(modle_tools)
endif()

add_subdirectory(libmodle)
add_subdirectory(modle)
