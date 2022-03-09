<!--
Copyright (C) 2022 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# MoDLE

[![Unit tests Ubuntu](https://github.com/paulsengroup/modle/actions/workflows/unit-tests-ubuntu.yml/badge.svg?branch=main)](https://github.com/paulsengroup/modle/actions/workflows/unit-tests-ubuntu.yml)
[![Unit tests Ubuntu](https://github.com/paulsengroup/modle/actions/workflows/unit-tests-ubuntu.yml/badge.svg?branch=devel)](https://github.com/paulsengroup/modle/actions/workflows/unit-tests-ubuntu.yml)
[![Build Docker image](https://github.com/paulsengroup/modle/actions/workflows/build-docker-image.yml/badge.svg)](https://github.com/paulsengroup/modle/actions/workflows/build-docker-image.yml)

## Using MoDLE

The recommended way to run MoDLE is using the official containers hosted on [dockerhub](https://hub.docker.com/repository/docker/robomics/modle) and [ghcr.io](https://github.com/robomics/modle/pkgs/container/modle)

```bash
sudo docker run robomics/modle:main --help
```

## Building MoDLE

### Pre-requisites

Building MoDLE requires a compiler toolchain supporting C++17, such as:

- GCC 7 and newer
- Clang 7 and newer
- Apple-Clang 10.0 and newer

MoDLE is being developed on Linux and should run on most UNIX-like OS (including macOS).

It is in theory possible to compile and run MoDLE on Windows, but we don't officially support this.

In addition to a C++17 compiler, building MoDLE requires the following tools:

- CMake >= 3.11
- Conan (recommended v1.30 or newer)

### Installing Conan

Conan can be installed with one of the following methods:
- `pip3 install "conan>=1.30"`
- `brew install conan`

### Compiling MoDLE

Run the following command from inside MoDLE's source tree.
The instructions assume the build machine has 8 CPU cores.
Feel free to adjust `-j 8` to match the number of CPU cores available on your machine.

```bash
mkdir build
cd build

cmake ..

cmake --build . -j 8
```

MoDLE's binaries will be located inside the `build/modle/` folder.

### Testing MoDLE

```bash
ctest -j 8                 \
      --test-dir .         \
      --schedule-random    \
      --output-on-failure  \
      --no-tests=error
```

### Installing MoDLE

This will install MoDLE files unde the prefix specified through `-DCMAKE_INSTALL_PREFIX` (`/usr/local` by default)

```
cmake --install .
```

### Troubleshooting

#### Incorrect or incomplete Conan profile
In this case the build process will fail during project configuration (i.e. when running CMake) with an error message similar to the following:

```
ERROR: libBigWig/0.4.6: 'settings.compiler' value not defined
CMake Error at build/conan.cmake:631 (message):
  Conan install failed='1'
```

This issue is usually fixed by forcing Conan to re-detect compiler information:

```bash
# Backup old profile
mv ~/.conan/profiles/default ~/.conan/profiles/default.bak
conan profile new ~/.conan/profiles/default --detect
```

If you see a warning about `GCC OLD ABI COMPATIBILITY` run:
```bash
conan profile update settings.compiler.libcxx=libstdc++11 default
```

On a Linux x86_64 machine with GCC 11 installed, the default profile should look similar to:
```
[settings]
os=Linux
os_build=Linux
arch=x86_64
arch_build=x86_64
compiler=gcc
compiler.version=11
compiler.libcxx=libstdc++11
build_type=Release
[options]
[build_requires]
[env]
```

On a Mac with Intel CPU the profile will be like:
```
[settings]
os=Macos
os_build=Macos
arch=x86_64
arch_build=x86_64
compiler=apple-clang
compiler.version=11.0
compiler.libcxx=libc++
build_type=Release
[options]
[build_requires]
[env]
```
