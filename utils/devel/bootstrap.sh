#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# Requires python3 and a compiler (gcc, clang or apple clang) installed
# TODO: allow users to specify a different /tmp dir.
# TODO: setup conan to write everything (including data/) inside the tmp dir
set -e


wd=/tmp/modle_tmp/

if [ $# -gt 2 ]; then
  echo "Usage: $0 /path/to/tmp/dir"
fi

if [ $# -eq 2 ]; then
  wd="$1"
  echo "Will use \"$wd\" as tmp directory"
fi

modle_src="$(readlink -f "$(dirname "$0")/..")"

rm -rf "$wd"
mkdir -p "$wd"
python3 -m venv "$wd/modle_venv"
PATH="$wd/modle_venv/bin:$PATH"
pip3 install --upgrade pip
pip3 --no-cache-dir install conan
conan profile new default --detect
conan profile update settings.compiler.libcxx=libstdc++11 default
echo -e "[requires]\ncmake/3.19.1\n\n[generators]\ncmake\n" > "$wd/conanfile.txt"
mkdir "$wd/build_dir" && cd "$wd/build_dir"
conan install "$wd/conanfile.txt"
cmake_root="$(grep 'CMAKE_ROOT=' conanbuildinfo.txt)"
PATH="$cmake_root/bin:$PATH"

mkdir cmake-build-dir && cd cmake-build-dir
cmake -DCMAKE_BUILD_TYPE=Release \
      -DENABLE_IPO=ON \
      -DENABLE_PCH=OFF \
      -DENABLE_TESTING=OFF \
      -DWARNINGS_AS_ERRORS=OFF \
      -G "CodeBlocks - Unix Makefiles" \
      "$modle_src"

make -j "$(nproc)"
# make install

