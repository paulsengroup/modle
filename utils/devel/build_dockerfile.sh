#!/usr/bin/env bash

# Copyright (c) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail

img_name='modle'
ver="$(git rev-parse --short HEAD)"
if ! git diff-index --quiet HEAD --; then
  ver+="-dirty"
fi

echo "Building \"$img_name:$ver\""

sudo docker pull docker.io/library/ubuntu:22.04

sudo docker build \
  --build-arg "BUILD_BASE_IMAGE=ghcr.io/paulsengroup/ci-docker-images/modle/ubuntu-22.04-cxx-clang-14:latest" \
  --build-arg "TEST_BASE_IMAGE=ghcr.io/paulsengroup/ci-docker-images/modle/ubuntu-22.04-cxx-clang-14:latest" \
  --build-arg "FINAL_BASE_IMAGE=docker.io/library/ubuntu" \
  --build-arg "FINAL_BASE_IMAGE_TAG=22.04" \
  --build-arg "FINAL_BASE_IMAGE_DIGEST=$(sudo docker inspect --format='{{index .RepoDigests 0}}' docker.io/library/ubuntu:22.04 | grep -o '[[:alnum:]:]\+$')" \
  --build-arg "C_COMPILER=clang-14" \
  --build-arg "CXX_COMPILER=clang++-14" \
  --build-arg "GIT_HASH=$(git rev-parse HEAD)" \
  --build-arg "GIT_SHORT_HASH=$(git rev-parse --short HEAD)" \
  --build-arg "CREATION_DATE=$(date --iso-8601)" \
  -t "$img_name:latest" \
  -t "$img_name:$(date --iso-8601 | tr -d '\-' )" \
  -t "$img_name:$ver" \
  "$(git rev-parse --show-toplevel)"

 # sudo singularity build -F "${img_name}_v${ver}.sif" \
 #                           "docker-daemon://${img_name}:${ver}"
