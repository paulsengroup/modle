# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(absl QUIET REQUIRED)
find_package(Boost QUIET REQUIRED)
find_package(CLI11 QUIET REQUIRED)
find_package(fmt QUIET REQUIRED)
find_package(spdlog QUIET REQUIRED)
find_package(tomlplusplus QUIET REQUIRED)
find_package(xxHash QUIET REQUIRED)

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
          tomlplusplus::tomlplusplus
          Modle::libmodle
          Modle::common)

include(GNUInstallDirs)
install(
  TARGETS modle
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${PROJECT_SOURCE_DIR}/LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/modle/)
