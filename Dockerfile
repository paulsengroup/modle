FROM fedora:33 AS modle_base

ENV SHELL=/usr/bin/bash
ENV PATH='/usr/bin:/usr/local/bin'

LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version="$ver"
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/modle"]

COPY cmake                  /tmp/modle/cmake
COPY external/Xoshiro-cpp   /tmp/modle/external/Xoshiro-cpp
COPY src                    /tmp/modle/src
COPY CMakeLists.txt         /tmp/modle/CMakeLists.txt
COPY LICENSE                /tmp/modle/LICENSE

ARG build_dir='/tmp/modle/cmake-build'
# march, cpus ver and build_type are set through --build-arg(s) at build time
ARG march
ARG cpus
ARG ver
ARG build_type

ARG CONAN_VER=1.34.1

ENV CC=/usr/bin/gcc
ENV CXX=/usr/bin/g++
ENV LD=/usr/bin/ld

# TODO: Enable tests

# Update system repo and install required tools
RUN dnf update -y \
    && dnf install -y --setopt=install_weak_deps=False --best            \
                      bash cmake gcc-c++ git make python3 zlib-devel     \
    && python3 -m venv /conan_venv --upgrade-deps                        \
    && /conan_venv/bin/pip3 --no-cache-dir install conan==${CONAN_VER}   \
    && mkdir -p /usr/local/bin                                           \
    && ln -s /conan_venv/bin/conan /usr/local/bin/conan                  \
    && mkdir "$build_dir" && cd "$build_dir"       \
    && env CC=/usr/bin/gcc                         \
           CXX=/usr/bin/g++                        \
           LD=/usr/bin/ld                          \
       cmake -DCMAKE_BUILD_TYPE=$build_type        \
             -DENABLE_IPO=ON                       \
             -DWARNINGS_AS_ERRORS=ON               \
             -DENABLE_TESTING=OFF                  \
             -DCMAKE_INSTALL_PREFIX='/usr/local'   \
             -DCMAKE_CXX_FLAGS="-march=${march}"   \
             -DCMAKE_C_FLAGS="-march=${march}"     \
             -G 'Unix Makefiles' ..                \
    && make -j "$cpus" install                     \
    &&                                             \
    if [ "$build_type" = "Debug" ]; then           \
        dnf install -y dnf-plugins-core gdb        \
     && dnf debuginfo-install -y libstdc++ zlib ;  \
    else                                           \
        unlink /usr/local/bin/conan                \
     && rm -r /conan_venv /root/.conan             \
     && dnf remove -y                              \
        cmake make gcc-c++ git zlib-devel          \
     && dnf clean all ;                            \
    fi
