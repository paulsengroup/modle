find_package(absl REQUIRED)
find_package(Boost REQUIRED)
find_package(concurrentqueue REQUIRED)
find_package(cpp-sort REQUIRED)
find_package(spdlog REQUIRED)

add_library(libmodle_cpu)
add_library(Modle::libmodle_cpu ALIAS libmodle_cpu)
set_target_properties(libmodle_cpu PROPERTIES OUTPUT_NAME modle_cpu)

target_sources(
  libmodle_cpu
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/simulation.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_correct_moves.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_detect_collisions.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/scheduler_common.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/scheduler_simulate.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/scheduler_perturbate.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/scheduler_replay.cpp)

target_include_directories(libmodle_cpu PUBLIC include/)

target_link_libraries(
  libmodle_cpu
  PRIVATE project_warnings
          project_options
          absl::span
          cpp-sort::cpp-sort
          Modle::cmatrix
          Modle::common
  PUBLIC absl::btree
         absl::fixed_array
         absl::flat_hash_map
         absl::flat_hash_set
         absl::strings
         absl::time
         Boost::filesystem
         concurrentqueue::concurrentqueue
         spdlog::spdlog
         Modle::libmodle_io
         Modle::libmodle_internal
         thread-pool)

if(ENABLE_PCH)
  target_precompile_headers(
    libmodle_cpu
    PRIVATE
    <algorithm>
    <cstddef>
    <cstdint>
    <filesystem>
    <limits>
    <memory>
    <string>
    <string_view>
    <utility>
    <vector>
    <absl/types/span.h>
    <fmt/format.h>)
endif()

install(
  TARGETS libmodle_cpu
  ARCHIVE
  INCLUDES
  LIBRARY)
