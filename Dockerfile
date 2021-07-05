FROM fedora:34 AS modle_base

ARG build_dir='/tmp/modle/cmake-build'
# The following args are set through --build-arg(s) at build time
ARG cpus
ARG ver
ARG build_type
ARG skip_tests

ARG XOSHIRO_CPP_VER=1.1
ARG LIBBIGWIG_VER=0.4.6
ARG CONAN_VER=1.38.0

# Required by new bincrafter artifactory
ARG CONAN_REVISIONS_ENABLED=1

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
COPY "external/libBigWig-$LIBBIGWIG_VER.tar.xz"                \
     "/tmp/modle/external/libBigWig-$LIBBIGWIG_VER.tar.xz"
COPY "external/Xoshiro-cpp-$XOSHIRO_CPP_VER.tar.xz"            \
     "/tmp/modle/external/Xoshiro-cpp-$XOSHIRO_CPP_VER.tar.xz"
COPY src                    /tmp/modle/src
COPY test                   /tmp/modle/test
COPY CMakeLists.txt         /tmp/modle/CMakeLists.txt
COPY conanfile.py           /tmp/modle/conanfile.py
COPY LICENSE                /tmp/modle/LICENSE

# Update system repo and install required tools
RUN dnf update -y \
    && dnf install -y --setopt=install_weak_deps=False --best            \
                      bash cmake gcc-c++ git make python3                \
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
             -G 'Unix Makefiles' ..                \
    && cmake --build "$build_dir" -j "$cpus"       \
    &&                                             \
    if [ ! "$skip_tests" = true ]; then            \
        dnf install -y python3-scipy               \
        && cd ..                                   \
        && ctest -j "$cpus"                        \
                 --test-dir "$build_dir"           \
                 --schedule-random                 \
                 --output-on-failure               \
        && dnf remove -y python3-scipy;            \
    fi                                             \
    && cmake --install "$build_dir"                \
    &&                                             \
    if [ "$build_type" = "Debug" ]; then           \
        dnf install -y dnf-plugins-core gdb        \
     && dnf debuginfo-install -y libstdc++ zlib ;  \
    else                                           \
        unlink /usr/local/bin/conan                \
     && rm -r /conan_venv /root/.conan /tmp/modle  \
     && dnf remove -y                              \
        cmake make gcc-c++ git zlib-devel          \
     && dnf clean all ;                            \
    fi
