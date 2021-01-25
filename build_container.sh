#!/usr/bin/env bash

set -e

memory='8G'
cpus="$(nproc)"
march='x86-64'

for arg in "$@"; do

  key=$(echo "$arg" | cut -f1 -d=)
  value=$(echo "$arg" | cut -f2 -d=)

  case "$key" in
  --name)        name=${value}   ;;
  --version)     ver=${value}    ;;
  --build-type)  type=${value}   ;;
  --cpus)        cpus=${value}   ;;
  --memory)      memory=${value} ;;
  --march)       march=${value}  ;;
  *) ;;
  esac
done

if [ -z "$name" ] || [ -z "$type" ]; then
  echo "Usage: ./build_container --name=modle --version=a.b.c --build-type=Release|Debug [ --cpus $(nproc) --memory 8G --march='x86-64' ]"
  exit 1
fi

sudo systemctl start docker.service

if [ -z "$ver" ]; then
  ver="$(git rev-parse --short HEAD)"
fi

if [ "${type,,}" = "debug" ]; then
  sudo docker build --memory="${memory}" \
    --compress \
    -t "robomics/${name,,}_dbg:latest" \
    -t "robomics/${name,,}_dbg:$ver" \
    -f "$(pwd)/Dockerfile.debug" \
    --cpu-period="100000" \
    --cpu-quota="$((100000 * $cpus))" \
    --build-arg "march=$march" \
    --build-arg "cpus=$cpus" \
    --build-arg "ver=$ver" \
    .

   sudo singularity build -F "${name}_v${ver}_dbg-${march}.sif" \
                          "docker-daemon://robomics/${name,,}_dbg:${ver}"
else
  sudo docker build --memory="${memory}" \
    --compress \
    -t "robomics/${name,,}:latest" \
    -t "robomics/${name,,}:$ver" \
    -f "$(pwd)/Dockerfile.release" \
    --cpu-period="100000" \
    --cpu-quota="$((100000 * $cpus))" \
    --build-arg "march=$march" \
    --build-arg "cpus=$cpus" \
    .
   sudo singularity build -F "${name}_v${ver}-${march}.sif" \
                          "docker-daemon://robomics/${name,,}:${ver}"
fi

# sudo systemctl stop docker.service
