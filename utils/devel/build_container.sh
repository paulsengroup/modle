#!/usr/bin/env bash

# Copyright (c) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e

if ! systemctl is-active --quiet docker.service; then
  echo "Starting docker.service..."
  sudo systemctl start --quiet docker.service
  stop_docker=true
else
  stop_docker=false
fi

ver="$(git rev-parse --short HEAD)"
if ! git diff-index --quiet HEAD --; then
  ver+="-dirty"
fi

img_name='modle'

sudo docker build \
     --compress \
     -t "robomics/${img_name}:latest" \
     -t "robomics/${img_name}:$ver" \
     -f "$(pwd)/Dockerfile" \
     --build-arg "ver=$ver" \
     .

 sudo singularity build -F "${img_name}_v${ver}.sif" \
                           "docker-daemon://robomics/${img_name}:${ver}"

# Restore docker.service to the original state
if $stop_docker ; then
  echo "Stopping docker.service..."
  sudo systemctl stop --quiet docker.service
fi
