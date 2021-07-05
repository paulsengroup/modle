#!/usr/bin/env bash

set -e

memory='8G'
cpus="$(nproc)"
name=modle
skip_tests=false

for arg in "$@"; do

  key=$(echo "$arg" | cut -f1 -d=)
  value=$(echo "$arg" | cut -f2 -d=)

  case "$key" in
  --name)        name=${value}   ;;
  --version)     ver=${value}    ;;
  --build-type)  type=${value}   ;;
  --cpus)        cpus=${value}   ;;
  --memory)      memory=${value} ;;
  --no-test)     skip_tests=true ;;
  *) ;;
  esac
done

if [[ $# -eq 1 && ("$1" == *help || "$1" == *h ) ]] ; then
  echo "Usage: ./build_container [ --name=modle --version=a.b.c --build-type=Release|Debug --cpus $(nproc) --memory 8G ]"
  exit 1
fi

if ! systemctl is-active --quiet docker.service; then
  echo "Starting docker.service..."
  sudo systemctl start --quiet docker.service
  stop_docker=true
else
  stop_docker=false
fi

if [ -z "$ver" ]; then
  ver="$(git rev-parse --short HEAD)"
  if ! git diff-index --quiet HEAD --; then
    ver+="-dirty"
    fi

fi

if [ -z "$type" ]; then
  type="Release"
fi

case "${type,,}" in
  "")                type=Release; echo "No build type provided. Defaulting to Release"   ;;
  debug)             type=Debug             ;;
  release)           type=Release           ;;
  relwithdebinfo)    type=RelWithDebInfo    ;;
  *)
    echo "Unknown build type '$type'. Allowed build types are Release, Debug and RelWithDebInfo"
    return 1
    ;;
esac

img_name="${name,,}"
if [ "$type" != "Release" ]; then
  img_name="$img_name-${type,,}"
fi

sudo docker build --memory="${memory}" \
     --compress \
     -t "robomics/${img_name}:latest" \
     -t "robomics/${img_name}:$ver" \
     -f "$(pwd)/Dockerfile" \
     --cpu-period="100000" \
     --cpu-quota="$((100000 * cpus))" \
     --build-arg "cpus=$cpus" \
     --build-arg "ver=$ver" \
     --build-arg "build_type=$type" \
     --build-arg "skip-tests=$skip_tests" \
     .

 sudo singularity build -F "${img_name}_v${ver}.sif" \
                           "docker-daemon://robomics/${img_name}:${ver}"

# Restore docker.service to the original state
if $stop_docker ; then
  echo "Stopping docker.service..."
  sudo systemctl stop --quiet docker.service
fi
