# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_subdirectory(internal)
add_subdirectory(cpu)

add_library(Modle::libmodle ALIAS libmodle_cpu)
