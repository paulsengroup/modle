#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail

repo_root="$(git rev-parse --show-toplevel)"
dockerfile="$repo_root/Dockerfile"

C_COMPILER='clang-14'
CXX_COMPILER='clang++-14'
FINAL_BASE_IMAGE_TAG='22.04'

BUILD_BASE_IMAGE="ghcr.io/paulsengroup/ci-docker-images/modle/ubuntu-$FINAL_BASE_IMAGE_TAG-cxx-$C_COMPILER:latest"
TEST_BASE_IMAGE="ghcr.io/paulsengroup/ci-docker-images/modle/ubuntu-$FINAL_BASE_IMAGE_TAG-cxx-$C_COMPILER:latest"
FINAL_BASE_IMAGE='docker.io/library/ubuntu'

FINAL_BASE_IMAGE_DIGEST="$(sudo docker inspect --format='{{index .RepoDigests 0}}' "$FINAL_BASE_IMAGE:$FINAL_BASE_IMAGE_TAG" | grep -o '[[:alnum:]:]\+$')"

GIT_HASH="$(git rev-parse HEAD)"
GIT_SHORT_HASH="$(git rev-parse --short HEAD)"
CREATION_DATE="$(date --iso-8601)"
#VERSION="0.0.1" # TODO changeme


sudo docker build -f "$dockerfile" \
    --build-arg "BUILD_BASE_IMAGE=$BUILD_BASE_IMAGE"               \
    --build-arg "TEST_BASE_IMAGE=$TEST_BASE_IMAGE"                 \
    --build-arg "FINAL_BASE_IMAGE=$FINAL_BASE_IMAGE"               \
    --build-arg "FINAL_BASE_IMAGE_TAG=$FINAL_BASE_IMAGE_TAG"       \
    --build-arg "FINAL_BASE_IMAGE_DIGEST=$FINAL_BASE_IMAGE_DIGEST" \
    --build-arg "C_COMPILER=$C_COMPILER"                           \
    --build-arg "CXX_COMPILER=$CXX_COMPILER"                       \
    --build-arg "GIT_HASH=$GIT_HASH"                               \
    --build-arg "GIT_SHORT_HASH=$GIT_SHORT_HASH"                   \
    --build-arg "CREATION_DATE=$CREATION_DATE"                     \
    -t 'modle:latest'                                              \
    -t "modle:${VERSION:-sha-$GIT_SHORT_HASH}"                     \
    -t "modle:${CREATION_DATE//-/}"                                \
    "$repo_root"

