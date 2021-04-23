FROM fedora:33 AS modle_base

ARG build_dir='/tmp/modle/cmake-build'
# The following args are set through --build-arg(s) at build time
ARG march
ARG cpus
ARG ver
ARG build_type
ARG skip_tests

ARG XOSHIRO_CPP_VER=1.1
ARG LIBBIGWIG_VER=0.4.6
ARG CONAN_VER=1.35.0

ENV CC=/usr/bin/gcc
ENV CXX=/usr/bin/g++
ENV LD=/usr/bin/ld
ENV SHELL=/usr/bin/bash
ENV PATH='/usr/bin:/usr/local/bin'

LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version="$ver"
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/modle"]

COPY cmake                  /tmp/modle/cmake
COPY "external/libBigWig-$LIBBIGWIG_VER"            \
     "/tmp/modle/external/libBigWig-$LIBBIGWIG_VER"
COPY "external/Xoshiro-cpp-$XOSHIRO_CPP_VER"        \
     "/tmp/modle/external/Xoshiro-cpp-$XOSHIRO_CPP_VER"
COPY src                    /tmp/modle/src
COPY test                   /tmp/modle/test
COPY CMakeLists.txt         /tmp/modle/CMakeLists.txt
COPY LICENSE                /tmp/modle/LICENSE

# Update system repo and install required tools
RUN dnf update -y \
    && dnf install -y --setopt=install_weak_deps=False --best            \
                      bash cmake gcc-c++ git make python3 python3-scipy  \
                      zlib-devel                                         \
    && python3 -m venv /conan_venv --upgrade-deps                        \
    && /conan_venv/bin/pip3 --no-cache-dir install conan==${CONAN_VER}   \
    && mkdir -p /usr/local/bin                                           \
    && ln -s /conan_venv/bin/conan /usr/local/bin/conan                  \
    && mkdir "$build_dir" && cd "$build_dir"       \
    && cmake -DCMAKE_BUILD_TYPE=$build_type        \
             -DENABLE_IPO=ON                       \
             -DWARNINGS_AS_ERRORS=ON               \
             -DENABLE_TESTING=ON                   \
             -DCMAKE_INSTALL_PREFIX='/usr/local'   \
             -DCMAKE_CXX_FLAGS="-march=${march}"   \
             -DCMAKE_C_FLAGS="-march=${march}"     \
             -G 'Unix Makefiles' ..                \
    && make modle modle_tools -j "$cpus" \
    &&                                   \
    if [ "$skip_tests" = false ]; then   \
        make test ;                      \
    fi                                   \
    && make modle modle_tools install    \
    &&                                             \
    if [ "$build_type" = "Debug" ]; then           \
        dnf install -y dnf-plugins-core gdb        \
     && dnf debuginfo-install -y libstdc++ zlib ;  \
    else                                           \
        unlink /usr/local/bin/conan                \
     && rm -r /conan_venv /root/.conan /tmp/modle  \
     && dnf remove -y                              \
        cmake make gcc-c++ git python3-scipy       \
        zlib-devel                                 \
     && dnf clean all ;                            \
    fi
