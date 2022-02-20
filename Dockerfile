# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ubuntu:20.04 AS builder

ARG src_dir='/root/modle'
ARG build_dir='/root/modle/build'
ARG staging_dir='/root/modle/staging'
ARG install_dir='/usr/local'

ARG BITFLAGS_VER=1.5.0
ARG LIBBIGWIG_VER=0.4.6
ARG THREAD_POOL_VER=2.0.0
ARG XOSHIRO_CPP_VER=1.1

ARG CONAN_VERSION='1.45.*'
ENV CONAN_V2=1
ENV CONAN_REVISIONS_ENABLED=1
ENV CONAN_NON_INTERACTIVE=1
ENV CONAN_CMAKE_GENERATOR=Ninja

ENV CC=clang-12
ENV CXX=clang++-12

RUN apt-get update                            \
&& apt-get install -y --no-install-recommends \
                   cmake                      \
                   clang-12                   \
                   make                       \
                   ninja-build                \
                   python3-pip

RUN pip3 install "conan==$CONAN_VERSION"

RUN mkdir -p "$src_dir" "$build_dir"

COPY conanfile.py "$src_dir"


RUN conan profile new "$HOME/.conan/profiles/default" --detect        \
&& conan profile update settings.compiler.libcxx=libstdc++11 default  \
&& conan profile update settings.compiler.cppstd=17 default

RUN cd "$build_dir"                          \
&& conan install "$src_dir/conanfile.py"     \
              --build outdated               \
              -s build_type=Release

COPY LICENSE                "$src_dir/LICENSE"
COPY "external/bitflags-$BITFLAGS_VER.tar.xz"                  \
     "$src_dir/external/bitflags-$BITFLAGS_VER.tar.xz"
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

RUN cd "$build_dir"                            \
&& cmake -DCMAKE_BUILD_TYPE=Release            \
         -DENABLE_DEVELOPER_MODE=OFF           \
         -DENABLE_TESTING=ON                   \
         -DCMAKE_INSTALL_PREFIX="$staging_dir" \
         -G Ninja                              \
         "$src_dir"

RUN cd "$build_dir"               \
&& cmake --build . -j "$(nproc)"  \
&& cmake --install .

FROM ubuntu:20.04 AS testing

ARG SCIPY_VER="1.5.1"
ARG WCORR_VER="1.9.5"
ARG src_dir="/root/modle"

RUN ln -snf /usr/share/zoneinfo/CET /etc/localtime \
&& echo CET | tee /etc/timezone > /dev/null

RUN apt-get update -q \
&& apt-get install -y -q --no-install-recommends cmake                      \
                                                 curl                       \
                                                 dirmngr                    \
                                                 python3-pip                \
                                                 software-properties-common

RUN curl -L 'https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc' | \
    tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc > /dev/null \
&& add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
&& add-apt-repository -y ppa:c2d4u.team/c2d4u4.0+

RUN apt-get install -y -q --no-install-recommends r-base        \
                                                  r-base-dev    \
                                                  r-cran-minqa  \
                                                  r-cran-mnormt \
                                                  r-cran-rcpparmadillo

RUN echo "options(Ncpus = $(nproc))" | tee "$HOME/.Rprofile" > /dev/null \
&& Rscript --no-save -e 'install.packages("wCorr", dependencies=c("Depends", "Imports", "LinkingTo"), repos="https://cloud.r-project.org")' \
&& Rscript --no-save -e 'quit(status=!library("wCorr", character.only=T, logical.return=T), save="no")'

RUN pip3 install "scipy==${SCIPY_VER}"

COPY --from=builder "$src_dir" "$src_dir"
RUN cd "$src_dir/build"           \
&& ctest -j "$(nproc)"            \
         --test-dir .             \
         --schedule-random        \
         --output-on-failure      \
         --no-tests=error         \
         --timeout 60             \
         --repeat after-timeout:3 \
&& rm -rf "$src_dir/test/Testing"

FROM ubuntu:20.04 AS base

ARG staging_dir='/root/modle/staging'
ARG install_dir='/usr/local'
ARG ver

COPY --from=testing "$staging_dir" "$install_dir"

LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version="$ver"
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/modle"]

RUN modle --help
RUN modle_tools --help
