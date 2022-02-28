#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -x

# Usage: prime_clion_remote_dev.sh user@hostname remote_path local_path
# Example: ./prime_clion_remote_dev.sh myuser@myserver modle/cmake-builds/cmake-build-release-gcc8 modle/cmake-builds/cmake-build-release-gcc8

host="$1"
remote_dir="$2"
local_dir="$3"


# Copy directory layout only
mkdir -p "$local_dir/CMakeFiles/"
rsync -av -f"+ */" -f"- *" "$host:$remote_dir/CMakeFiles/*.*.*" "$local_dir/CMakeFiles/"

rsync -aPv "$host:$remote_dir/CMakeCache.txt" "$local_dir"
rsync -aPv "$host:$remote_dir/MakeFiles/TargetDirectories.txt" "$local_dir/CMakeFiles/"
rsync -aPv "$host:$remote_dir/CMakeFiles/Makefile.cmake" "$local_dir/CMakeFiles/"
rsync -aPv "$host:$remote_dir/CMakeFiles/*.*.*/CMakeCXXCompiler.cmake" "$local_dir/CMakeFiles/"*.*.*/
