# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl REQUIRED)
find_package(Boost REQUIRED)
find_package(CLI11 REQUIRED)
find_package(fmt REQUIRED)
find_package(spdlog REQUIRED)
find_package(xxHash REQUIRED)

add_executable(modle)

target_sources(modle PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cli.cpp ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

target_link_libraries(
  modle
  PRIVATE project_warnings
          project_options
          absl::failure_signal_handler
          absl::strings
          absl::symbolize
          absl::time
          Boost::filesystem
          Boost::headers
          CLI11::CLI11
          fmt::fmt
          spdlog::spdlog
          Modle::libmodle
          Modle::common)
# Modle::libmodle_io)

install(
  TARGETS modle
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)
