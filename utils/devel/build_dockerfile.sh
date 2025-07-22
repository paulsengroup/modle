#!/usr/bin/env bash

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail

IMAGE_NAME='modle'

if [ "$(uname)" == "Darwin" ]; then
  BUILD_USER="$USER"
else
  BUILD_USER='root'
fi

GIT_HASH="$(git rev-parse HEAD)"
GIT_SHORT_HASH="$(git rev-parse --short HEAD)"
GIT_TAG="$(git for-each-ref 'refs/tags/v*.*.*' --count 1 --sort=-v:refname --format "%(refname:short)" --points-at HEAD)"
CREATION_DATE="$(date -I)"

if [[ $(git status --porcelain -uno) ]]; then
  GIT_IS_DIRTY=1
else
  GIT_IS_DIRTY=0
fi

IMAGE_TAG="sha-$GIT_SHORT_HASH"
if [ $GIT_IS_DIRTY -ne 0 ]; then
  IMAGE_TAG+='-dirty'
fi

if [ -z "$GIT_TAG" ]; then
  GIT_TAG="sha-$GIT_SHORT_HASH"
fi

2>&1 echo "Building \"$IMAGE_NAME:$IMAGE_TAG\"..."

sudo -u "$BUILD_USER" docker pull docker.io/library/ubuntu:24.04
FINAL_BASE_IMAGE_DIGEST="$(sudo -u "$BUILD_USER" docker inspect --format='{{index .RepoDigests 0}}' docker.io/library/ubuntu:24.04 | grep -o '[[:alnum:]:]\+$')"

BUILD_BASE_IMAGE='ghcr.io/paulsengroup/ci-docker-images/ubuntu-24.04-cxx-clang-20:latest'

sudo -u "$BUILD_USER" docker pull "$BUILD_BASE_IMAGE"

# sudo -u "$BUILD_USER" docker buildx build --platform linux/amd64,linux/arm64 \
sudo -u "$BUILD_USER" docker buildx build \
  --build-arg "BUILD_BASE_IMAGE=$BUILD_BASE_IMAGE" \
  --build-arg "FINAL_BASE_IMAGE=docker.io/library/ubuntu" \
  --build-arg "FINAL_BASE_IMAGE_TAG=24.04" \
  --build-arg "FINAL_BASE_IMAGE_DIGEST=$FINAL_BASE_IMAGE_DIGEST" \
  --build-arg "C_COMPILER=clang" \
  --build-arg "CXX_COMPILER=clang++" \
  --build-arg "GIT_HASH=$GIT_HASH" \
  --build-arg "GIT_SHORT_HASH=$GIT_SHORT_HASH" \
  --build-arg "GIT_TAG=$GIT_TAG" \
  --build-arg "GIT_IS_DIRTY=$GIT_IS_DIRTY" \
  --build-arg "CREATION_DATE=$CREATION_DATE" \
  --build-arg 'LIBCXX=libc++' \
  -t "$IMAGE_NAME:latest" \
  -t "$IMAGE_NAME:$(echo "$CREATION_DATE" | tr -d '\-' )" \
  -t "$IMAGE_NAME:$IMAGE_TAG" \
  "$(git rev-parse --show-toplevel)" \
  --progress=plain

 # sudo -u "$BUILD_USER" apptainer build -F "${img_name}_v${ver}.sif" \
 #                           "docker-daemon://${img_name}:${ver}"
