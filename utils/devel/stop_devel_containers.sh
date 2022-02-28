#!/usr/bin/env bash

set -e
set -u

clang_versions=(6.0 7 8 9 10 11 12)
gcc_versions=(7 8 9 10)

for ver in "${clang_versions[@]}"; do
  label="modle-devel-clang-${ver}:latest"
  name="modle-devel-clang-${ver}"

  echo "Stopping $label..."
  sudo docker stop "$name"
  sudo docker rm "$name"
done

for ver in "${gcc_versions[@]}"; do
  label="modle-devel-gcc-${ver}:latest"
  name="modle-devel-gcc-${ver}"

  echo "Stopping $label..."
  sudo docker stop "$name"
  sudo docker rm "$name"
done
