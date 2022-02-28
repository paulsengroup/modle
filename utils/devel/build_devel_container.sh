#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -x
set -u

trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM EXIT

pkey="$1"
clang_versions=(6.0 7 8 9 10 11 12)
gcc_versions=(7 8 9 10)

mkdir -p containers

build_image () {
  set -e
  set -x
  name="$1"
  ver="$2"
  pkey="$3"
  label="modle-devel-${name}-${ver}:latest"
  docker build -t "$label"                          \
               -q                                   \
               -f devel.Dockerfile                  \
               --build-arg "COMPILER_NAME=$name"    \
               --build-arg "COMPILER_VERSION=$ver"  \
               --build-arg "AUTHORIZED_KEYS=$pkey"  \
               .

  imgname="containers/$(basename "$label" ':latest').tar.zst"
  rm -f "$imgname"
  docker save "$label" | zstd -T0 --adapt -o "$imgname"
}

export -f build_image

for ver in "${clang_versions[@]}"; do
  DOCKER_BUILDKIT=0 \
  sudo -E bash -c "$(declare -f build_image); build_image 'clang' \"$ver\" \"$pkey\"" &
done

for ver in "${gcc_versions[@]}"; do
  DOCKER_BUILDKIT=0 \
  sudo -E bash -c "$(declare -f build_image); build_image 'gcc' \"$ver\" \"$pkey\"" &
done

wait
