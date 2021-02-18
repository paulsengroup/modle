FROM julia:1.5.3 AS julia_base

ARG build_dir='/tmp/combine_bwig_stats'
ARG ver='1.0.0'

COPY src                 "$build_dir/src"
COPY dummy_data          "$build_dir/dummy_data"
COPY Manifest.toml       "$build_dir/Manifest.toml"
COPY Project.toml        "$build_dir/Project.toml"
COPY precompile_app.jl   "$build_dir/precompile_app.jl"
COPY build_app.sh        "$build_dir/build_app.sh"

# Update system repo and install required tools
RUN apt-get update            \
    && apt-get install -y gcc \
    && cd "$build_dir"        \
    && ./build_app.sh         \
    && mkdir -p /opt/         \
    && cp -r /tmp/combine_bwig_stats_app /opt \
    && apt-get remove -y gcc \
    && rm -rf /var/lib/apt/lists/* /root/.julia/*

ENV SHELL=/usr/bin/bash
ENV PATH='/usr/bin:/usr/bin/local:/opt/combine_bwig_stats_app/bin'

LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version="$ver"
WORKDIR /data
ENTRYPOINT ["compute_bwig_stats"]
