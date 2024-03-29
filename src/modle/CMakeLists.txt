# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED COMPONENTS headers)
find_package(CLI11 CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(tomlplusplus CONFIG REQUIRED)

add_executable(modle)

target_sources(modle PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cli.cpp ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

target_link_libraries(
  modle
  PRIVATE project_warnings
          project_options
          Modle::common
          Modle::config
          Modle::libmodle)

target_link_system_libraries(
  modle
  PRIVATE
  absl::failure_signal_handler
  absl::flat_hash_map
  absl::span
  absl::strings
  absl::symbolize
  absl::time
  Boost::headers
  CLI11::CLI11
  spdlog::spdlog
  tomlplusplus::tomlplusplus
  hictk)

include(GNUInstallDirs)
install(
  TARGETS modle
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${PROJECT_SOURCE_DIR}/LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/modle/)
