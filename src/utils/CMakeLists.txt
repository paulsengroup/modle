add_library(modle_utils INTERFACE)
add_library(Modle::utils ALIAS modle_utils)
target_include_directories(modle_utils INTERFACE include/)
