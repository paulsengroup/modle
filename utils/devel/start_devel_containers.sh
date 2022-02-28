#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u

clang_versions=(6.0 7 8 9 10 11 12)
gcc_versions=(7 8 9 10)

port=2222
for ver in "${clang_versions[@]}"; do
  label="modle-devel-clang-${ver}:latest"
  name="modle-devel-clang-${ver}"

  echo "Starting $label on port $port..."
  sudo docker run -d --cap-add sys_ptrace  \
                     -p 127.0.0.1:$port:22 \
                     --name "$name"        \
                     "$label"
  port=$(( port + 1 ))
done

for ver in "${gcc_versions[@]}"; do
  label="modle-devel-gcc-${ver}:latest"
  name="modle-devel-gcc-${ver}"

  echo "Starting $label on port $port..."
  sudo docker run -d --cap-add sys_ptrace  \
                     -p 127.0.0.1:$port:22 \
                     --name "$name"        \
                     "$label"
  port=$(( port + 1 ))
done
