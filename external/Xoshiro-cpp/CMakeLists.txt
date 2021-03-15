project(xoshiro-cpp CXX)

add_library(xoshiro-cpp INTERFACE)
add_library(XoshiroCpp::XoshiroCpp ALIAS xoshiro-cpp)

target_sources(xoshiro-cpp PRIVATE include/Xoshiro-cpp/XoshiroCpp.hpp)

target_include_directories(xoshiro-cpp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(xoshiro-cpp INTERFACE project_warnings project_options)
