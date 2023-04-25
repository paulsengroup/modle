# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

##### IMPORTANT #####
# This Dockerfile requires several build arguments to be defined through --build-arg
# See utils/devel/build_dockerfile.sh for an example of how to build this Dockerfile
#####################

ARG BUILD_BASE_IMAGE
ARG TEST_BASE_IMAGE
ARG FINAL_BASE_IMAGE
ARG FINAL_BASE_IMAGE_DIGEST

FROM "$BUILD_BASE_IMAGE" AS builder

ARG src_dir='/root/modle'
ARG build_dir='/root/modle/build'
ARG staging_dir='/root/modle/staging'
ARG install_dir='/usr/local'

ENV CONAN_CMAKE_GENERATOR=Ninja


ARG C_COMPILER
ARG CXX_COMPILER

RUN if [ -z "$C_COMPILER" ]; then echo "Missing C_COMPILER --build-arg" && exit 1; fi \
&&  if [ -z "$CXX_COMPILER" ]; then echo "Missing CXX_COMPILER --build-arg" && exit 1; fi

ENV CC="$C_COMPILER"
ENV CXX="$CXX_COMPILER"

# Build MoDLE's deps using Conan
RUN mkdir -p "$src_dir" "$build_dir"

COPY conanfile.txt "$src_dir"
RUN cd "$build_dir"                                  \
&& conan install "$src_dir/conanfile.txt"            \
                 --build=missing                     \
                 --build=cascade                     \
                 --update                            \
                 -pr:b="$CONAN_DEFAULT_PROFILE_PATH" \
                 -pr:h="$CONAN_DEFAULT_PROFILE_PATH" \
                 -s build_type=Release               \
                 -s compiler.libcxx=libstdc++11      \
                 -s compiler.cppstd=17               \
                 --output-folder="$build_dir"

# Copy source files
COPY LICENSE "$src_dir/"
COPY external "$src_dir/external/"
COPY cmake "$src_dir/cmake/"
COPY test/CMakeLists.txt "$src_dir/test/"
COPY CMakeLists.txt "$src_dir/"
COPY src "$src_dir/src/"
COPY test/units "$src_dir/test/units/"

ARG GIT_HASH
ARG GIT_SHORT_HASH
ARG GIT_TAG
ARG GIT_IS_DIRTY

RUN if [ -z "$GIT_HASH" ]; then echo "Missing GIT_HASH --build-arg" && exit 1; fi \
&&  if [ -z "$GIT_SHORT_HASH" ]; then echo "Missing GIT_SHORT_HASH --build-arg" && exit 1; fi \
&&  if [ -z "$GIT_IS_DIRTY" ]; then echo "Missing GIT_IS_DIRTY --build-arg" && exit 1; fi \
&&  if [ -z "$GIT_TAG" ]; then echo "Missing GIT_TAG --build-arg" && exit 1; fi

# Configure project
RUN cmake -DCMAKE_BUILD_TYPE=Release            \
          -DCMAKE_PREFIX_PATH="$build_dir"      \
          -DWARNINGS_AS_ERRORS=ON               \
          -DENABLE_DEVELOPER_MODE=OFF           \
          -DMODLE_ENABLE_TESTING=ON             \
          -DMODLE_DOWNLOAD_TEST_DATASET=OFF     \
          -DGIT_RETRIEVED_STATE=true            \
          -DGIT_TAG="$GIT_TAG"                  \
          -DGIT_IS_DIRTY="$GIT_IS_DIRTY"        \
          -DGIT_HEAD_SHA1="$GIT_HASH"           \
          -DGIT_DESCRIBE="$GIT_SHORT_HASH"      \
          -DCMAKE_INSTALL_PREFIX="$staging_dir" \
          -G Ninja                              \
          -S "$src_dir"                         \
          -B "$build_dir"

# Build and install project
RUN cmake --build "$build_dir" -j "$(nproc)"  \
&& cmake --install "$build_dir"

ARG TEST_BASE_IMAGE
FROM "$TEST_BASE_IMAGE" AS unit-testing

ARG src_dir="/root/modle"

COPY --from=builder "$src_dir" "$src_dir"
COPY test/data/modle_test_data.tar.xz "$src_dir/test/data/"

RUN tar -xf "$src_dir/test/data/modle_test_data.tar.xz" -C "$src_dir/"

RUN ctest -j "$(nproc)"               \
          --test-dir "$src_dir/build" \
          --schedule-random           \
          --output-on-failure         \
          --no-tests=error            \
          --timeout 300               \
&& rm -rf "$src_dir/test/Testing"

ARG FINAL_BASE_IMAGE
ARG FINAL_BASE_IMAGE_DIGEST
FROM "${FINAL_BASE_IMAGE}@${FINAL_BASE_IMAGE_DIGEST}" AS integration-testing

ARG src_dir="/root/modle"
ARG staging_dir='/root/modle/staging'
ARG install_dir='/usr/local'

RUN apt-get update \
&& apt-get install -q -y --no-install-recommends \
                   gcc                           \
                   libdigest-sha-perl            \
                   python3-dev                   \
                   python3-pip                   \
                   xz-utils

RUN pip3 install 'cooler>=0.9.1'

COPY --from=unit-testing "$staging_dir" "$staging_dir"
COPY test/data/modle_test_data.tar.xz "$src_dir/test/data/"
COPY test/scripts/modle*integration_test.sh "$src_dir/test/scripts/"

RUN tar -xf "$src_dir/test/data/modle_test_data.tar.xz" -C "$src_dir/"

RUN "$src_dir/test/scripts/modle_integration_test.sh" "$staging_dir/bin/modle"
RUN "$src_dir/test/scripts/modle_tools_transform_integration_test.sh" "$staging_dir/bin/modle_tools"
RUN "$src_dir/test/scripts/modle_tools_eval_integration_test.sh" "$staging_dir/bin/modle_tools"

ARG FINAL_BASE_IMAGE
ARG FINAL_BASE_IMAGE_DIGEST
FROM "${FINAL_BASE_IMAGE}@${FINAL_BASE_IMAGE_DIGEST}" AS base

ARG staging_dir='/root/modle/staging'
ARG install_dir='/usr/local'

ARG BUILD_BASE_IMAGE
ARG FINAL_BASE_IMAGE
ARG FINAL_BASE_IMAGE_DIGEST

ARG GIT_HASH
ARG GIT_SHORT_HASH
ARG VERSION
ARG CREATION_DATE

RUN if [ -z "$BUILD_BASE_IMAGE" ]; then echo "Missing BUILD_BASE_IMAGE --build-arg" && exit 1; fi \
&&  if [ -z "$FINAL_BASE_IMAGE" ]; then echo "Missing FINAL_BASE_IMAGE --build-arg" && exit 1; fi \
&&  if [ -z "$FINAL_BASE_IMAGE_DIGEST" ]; then echo "Missing FINAL_BASE_IMAGE_DIGEST --build-arg" && exit 1; fi \
&&  if [ -z "$GIT_HASH" ]; then echo "Missing GIT_HASH --build-arg" && exit 1; fi \
&&  if [ -z "$GIT_SHORT_HASH" ]; then echo "Missing GIT_SHORT_HASH --build-arg" && exit 1; fi \
&&  if [ -z "$CREATION_DATE" ]; then echo "Missing CREATION_DATE --build-arg" && exit 1; fi

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
