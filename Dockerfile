# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

##### IMPORTANT #####
# This Dockerfile requires several build arguments to be defined through --build-arg
# See utils/devel/build_dockerfile.sh for an example of how to build this Dockerfile
#####################


ARG BUILD_BASE_IMAGE=ghcr.io/paulsengroup/ci-docker-images/alpine-3.22-cxx
ARG FINAL_BASE_IMAGE=ubuntu
ARG FINAL_BASE_IMAGE_DIGEST=unknown

FROM ${BUILD_BASE_IMAGE} AS builder

ARG src_dir='/root/modle'
ARG build_dir='/root/modle/build'
ARG staging_dir='/root/modle/staging'
ARG install_dir='/usr/local'

ARG CONAN_CMAKE_GENERATOR=Ninja
ARG CMAKE_POLICY_VERSION_MINIMUM=3.5

# Install b2 using Conan
RUN printf '[requires]\nb2/5.3.3\n[options]\nb2*:toolset=%s' \
           "$(basename "$(which "$CC")" | cut -f 1 -d -)" > /tmp/conanfile.txt

RUN conan install /tmp/conanfile.txt                 \
                 --build='*'                         \
                 -pr:b="$CONAN_DEFAULT_PROFILE_PATH" \
                 -pr:h="$CONAN_DEFAULT_PROFILE_PATH"

# Build MoDLE's deps using Conan
RUN mkdir -p "$src_dir" "$build_dir"

COPY conanfile.py "$src_dir/conanfile.py"

# Make sure we're linking to libc++
RUN sed -i "s/^compiler\.libcxx.*$/compiler.libcxx=libc++/" "$CONAN_DEFAULT_PROFILE_PATH"

RUN conan install "$src_dir/conanfile.py"            \
                 --build='missing'                   \
                 -pr:b="$CONAN_DEFAULT_PROFILE_PATH" \
                 -pr:h="$CONAN_DEFAULT_PROFILE_PATH" \
                 -s build_type=Release               \
                 -s compiler.cppstd=20               \
                 --output-folder="$build_dir"

# Copy source files
COPY LICENSE "$src_dir/"
COPY external "$src_dir/external/"
COPY cmake "$src_dir/cmake/"
COPY test/CMakeLists.txt "$src_dir/test/"
COPY CMakeLists.txt "$src_dir/"
COPY src "$src_dir/src/"
COPY test/units "$src_dir/test/units/"

ARG GIT_HASH=0000000000000000000000000000000000000000
ARG GIT_SHORT_HASH=00000000
ARG GIT_TAG=unknown
ARG GIT_IS_DIRTY=false
ARG CMAKE_EXE_LINKER_FLAGS='-static -stdlib=libc++ LINKER:-lc++,-lc++abi'

# Configure project
RUN cmake -DCMAKE_BUILD_TYPE=Release                        \
         -DCMAKE_CXX_FLAGS='-stdlib=libc++'                 \
         -DCMAKE_EXE_LINKER_FLAGS="$CMAKE_EXE_LINKER_FLAGS" \
         -DCMAKE_LINKER_TYPE=LLD \
         -DCMAKE_C_COMPILER_AR="$(which ar)" \
         -DCMAKE_C_COMPILER_RANLIB="$(which ranlib)" \
         -DCMAKE_CXX_COMPILER_AR="$(which ar)" \
         -DCMAKE_CXX_COMPILER_RANLIB="$(which ranlib)" \
         -DCMAKE_PREFIX_PATH="$build_dir"                   \
         -DENABLE_DEVELOPER_MODE=OFF                        \
         -DMODLE_ENABLE_TESTING=ON                          \
         -DMODLE_DOWNLOAD_TEST_DATASET=OFF                  \
         -DMODLE_GIT_RETRIEVED_STATE=true                   \
         -DMODLE_GIT_TAG="$GIT_TAG"                         \
         -DMODLE_GIT_IS_DIRTY="$GIT_IS_DIRTY"               \
         -DMODLE_GIT_HEAD_SHA1="$GIT_HASH"                  \
         -DMODLE_GIT_DESCRIBE="$GIT_SHORT_HASH"             \
         -DCMAKE_INSTALL_PREFIX="$staging_dir"              \
         -G Ninja                                           \
         -S "$src_dir"                                      \
         -B "$build_dir"

# Build and install project
RUN cmake --build "$build_dir" -j "$(nproc)"  \
&& cmake --install "$build_dir"

ARG BUILD_BASE_IMAGE=ghcr.io/paulsengroup/ci-docker-images/alpine-3.22-cxx
FROM ${BUILD_BASE_IMAGE} AS unit-testing

ARG src_dir="/root/modle"

COPY --from=builder "$src_dir" "$src_dir"
COPY test/data/modle_test_data.tar.zst "$src_dir/test/data/"

RUN zstd -dc "$src_dir/test/data/modle_test_data.tar.zst" |  tar -xf - -C "$src_dir/"

RUN ctest -j "$(nproc)"               \
          --test-dir "$src_dir/build" \
          --schedule-random           \
          --output-on-failure         \
          --no-tests=error            \
          --timeout 900               \
&& rm -rf "$src_dir/test/Testing"


ARG FINAL_BASE_IMAGE='ubuntu:24.04'

FROM ${FINAL_BASE_IMAGE} AS integration-testing

ARG src_dir="/root/modle"
ARG staging_dir='/root/modle/staging'
ARG install_dir='/usr/local'

RUN apt-get update \
&& apt-get install -q -y --no-install-recommends \
                   gcc                           \
                   libdigest-sha-perl            \
                   python3-dev                   \
                   python3-pip                   \
                   python3-venv                  \
                   xz-utils                      \
                   zstd

RUN python3 -m venv /tmp/venv --upgrade \
&& /tmp/venv/bin/pip install 'cooler>=0.10.3' 'pyBigWig>=0.3.24'

COPY --from=unit-testing "$staging_dir" "$staging_dir"
COPY test/data/modle_test_data.tar.zst "$src_dir/test/data/"
COPY test/scripts/modle*integration_test.sh "$src_dir/test/scripts/"
COPY test/scripts/compare_bwigs.py "$src_dir/test/scripts/"

RUN tar -xf "$src_dir/test/data/modle_test_data.tar.zst" -C "$src_dir/"

ARG PATH="/tmp/venv/bin:$PATH"

RUN "$src_dir/test/scripts/modle_integration_test.sh" "$staging_dir/bin/modle"
RUN "$src_dir/test/scripts/modle_tools_transform_integration_test.sh" "$staging_dir/bin/modle_tools"
RUN "$src_dir/test/scripts/modle_tools_eval_integration_test.sh" "$staging_dir/bin/modle_tools"
RUN "$src_dir/test/scripts/modle_tools_annotate_barriers_integration_test.sh" "$staging_dir/bin/modle_tools"

ARG FINAL_BASE_IMAGE=ubuntu

FROM ${FINAL_BASE_IMAGE} AS base

ARG staging_dir='/root/modle/staging'
ARG install_dir='/usr/local'

ARG BUILD_BASE_IMAGE=ghcr.io/paulsengroup/ci-docker-images/alpine-3.22-cxx
ARG FINAL_BASE_IMAGE=ubuntu
ARG FINAL_BASE_IMAGE_DIGEST=unknown

ARG GIT_HASH=0000000000000000000000000000000000000000
ARG GIT_SHORT_HASH=00000000
ARG GIT_TAG=unknown
ARG GIT_IS_DIRTY=false
ARG VERSION=latest
ARG CREATION_DATE=2000-01-01

# Export project binaries to the final build stage
COPY --from=integration-testing "$staging_dir" "$install_dir"

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/modle"]

RUN modle --help
RUN modle_tools --help
RUN modle --version

# https://github.com/opencontainers/image-spec/blob/main/annotations.md#pre-defined-annotation-keys
LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/modle'
LABEL org.opencontainers.image.documentation='https://github.com/paulsengroup/modle'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/modle'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title='MoDLE'
LABEL org.opencontainers.image.description='High-performance stochastic modeling of DNA loop extrusion interactions'
LABEL org.opencontainers.image.base.digest="$FINAL_BASE_IMAGE_DIGEST"
LABEL org.opencontainers.image.base.name="$FINAL_BASE_IMAGE"
LABEL paulsengroup.modle.image.build-base="$BUILD_BASE_IMAGE"

LABEL org.opencontainers.image.revision="$GIT_HASH"
LABEL org.opencontainers.image.created="$CREATION_DATE"
LABEL org.opencontainers.image.version="${VERSION:-sha-$GIT_SHORT_HASH}"
