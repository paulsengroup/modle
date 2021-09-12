# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM conanio/clang11 AS builder

ARG src_dir='/home/conan/modle'
ARG build_dir='/home/conan/modle/build'
ARG staging_dir='/home/conan/modle/staging'
ARG install_dir='/usr/local'

ARG LIBBIGWIG_VER=0.4.6
ARG THREAD_POOL_VER=2.0.0
ARG XOSHIRO_CPP_VER=1.1
ARG SCIPY_VER=1.7.1

ARG CONAN_V2=1
ARG CONAN_REVISIONS_ENABLED=1
ARG CONAN_NON_INTERACTIVE=1
ARG CONAN_CMAKE_GENERATOR=Ninja

# Update system repo and install required tools
RUN sudo apt-get update                            \
    && sudo apt-get install -y ninja-build         \
    && pip install scipy=="${SCIPY_VER}"

RUN mkdir -p "$src_dir" "$build_dir"

COPY conanfile.py "$src_dir"

RUN sudo chown -R conan "$src_dir"

RUN cd "$build_dir"                              \
    && conan install "$src_dir/conanfile.py"     \
                  --build outdated               \
                  -s compiler.cppstd=17          \
                  -s build_type=RelWithDebInfo   \
                  -s compiler.libcxx=libstdc++11 \
                  -o enable_testing=ON

COPY LICENSE                "$src_dir/LICENSE"
COPY "external/libBigWig-$LIBBIGWIG_VER.tar.xz"                \
     "$src_dir/external/libBigWig-$LIBBIGWIG_VER.tar.xz"
COPY "external/mscharconv.tar.xz"                              \
     "$src_dir/external/mscharconv.tar.xz"
COPY "external/thread-pool-$THREAD_POOL_VER.tar.xz"            \
     "$src_dir/external/thread-pool-$THREAD_POOL_VER.tar.xz"
COPY "external/Xoshiro-cpp-$XOSHIRO_CPP_VER.tar.xz"            \
     "$src_dir/external/Xoshiro-cpp-$XOSHIRO_CPP_VER.tar.xz"
COPY cmake                  "$src_dir/cmake"
COPY CMakeLists.txt         "$src_dir/CMakeLists.txt"
COPY test                   "$src_dir/test"
COPY src                    "$src_dir/src"

RUN sudo chown -R conan "$src_dir"
RUN cd "$build_dir"                                \
    && cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo     \
             -DENABLE_IPO=ON                       \
             -DWARNINGS_AS_ERRORS=ON               \
             -DENABLE_TESTING=ON                   \
             -DCMAKE_INSTALL_PREFIX="$staging_dir" \
             -G Ninja                              \
             "$src_dir"

RUN cd "$build_dir"                  \
    && cmake --build . -j "$(nproc)" \
    && ctest -j "$(nproc)"           \
             --test-dir .            \
             --schedule-random       \
             --output-on-failure     \
             --no-tests=error        \
    && cmake --install .

FROM ubuntu:bionic AS base

ARG staging_dir='/home/conan/modle/staging'
ARG install_dir='/usr/local'
ARG ver

COPY --from=builder "$staging_dir" "$install_dir"

LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version="$ver"
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/modle"]

RUN modle --help
