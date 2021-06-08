option(MODLE_BUILD_UTILS "Enable compilation of ModLE's utilities (such as modle_tools)" ON)

if(MODLE_BUILD_UTILS)
  message(STATUS "Building ModLE utils")
  add_subdirectory(modle_tools)
endif()

add_subdirectory(common)
add_subdirectory(contact_matrix)
add_subdirectory(correlation)
add_subdirectory(libio)
add_subdirectory(libmodle)
add_subdirectory(modle)
