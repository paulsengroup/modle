# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ubuntu:20.04 AS base

ENV CONAN_V2=1
ENV CONAN_REVISIONS_ENABLED=1
ENV CONAN_NON_INTERACTIVE=1
ENV CONAN_CMAKE_GENERATOR=Ninja
ARG PIP_NO_CACHE_DIR=0

ARG COMPILER_NAME
ARG COMPILER_VERSION
ARG COMPILER="$COMPILER_NAME-$COMPILER_VERSION"
ARG AUTHORIZED_KEYS

RUN ln -snf /usr/share/zoneinfo/CET /etc/localtime \
&& echo CET | tee /etc/timezone > /dev/null

RUN apt-get update                \
&& apt-get install -y "$COMPILER" \
                      clang-tidy  \
                      cmake       \
                      cppcheck    \
                      curl        \
                      gdb         \
                      ninja-build \
                      python3-pip \
                      rsync       \
                      ssh         \
                      tar         \
&& if [ $COMPILER_NAME = gcc ] ; then \
    apt-get install -y "g++-${COMPILER_VERSION}"; \
fi \
&& rm -rf /var/lib/apt/lists/*

RUN pip install "conan==1.45.*"

RUN if [ $COMPILER_NAME = gcc ] ; then \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-$COMPILER_VERSION 10  \
&&  update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-$COMPILER_VERSION 10; \
fi

RUN if [ $COMPILER_NAME = clang ] ; then \
    update-alternatives --install /usr/bin/clang clang /usr/bin/clang-$COMPILER_VERSION 10        \
&&  update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-$COMPILER_VERSION 10; \
fi

RUN useradd -m clion

# COPY conanfile.py "/tmp"
# RUN chmod 444 /tmp/conanfile.py
#
USER clion
RUN conan profile new "$HOME/.conan/profiles/default" --detect       \
&& conan profile update settings.compiler.libcxx=libstdc++11 default \
&& conan profile update settings.compiler.cppstd=17 default
#
# RUN cd /tmp \
# && conan install "conanfile.py"       \
#               --build outdated        \
#               -s build_type=Debug     \
# && conan install "conanfile.py"       \
#               --build outdated        \
#               -s build_type=Release

USER root
RUN mkdir -p /root/.ssh /home/clion/.ssh \
&& echo "$AUTHORIZED_KEYS" >> /root/.ssh/authorized_keys \
&& echo "$AUTHORIZED_KEYS" >> /home/clion/.ssh/authorized_keys \
&& chown -R clion:clion /home/clion \
&& chmod 750 /root/.ssh/ /home/clion/.ssh/   \
&& chmod 640 /root/.ssh/* /home/clion/.ssh/*

RUN service ssh start
EXPOSE 22

LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL ver=1.0.0
WORKDIR /data

CMD ["/usr/sbin/sshd","-D"]
