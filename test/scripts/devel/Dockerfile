# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


# The image generated from this Dockerfile can be used to generate the datasets
# used to test modle and modle_tools

FROM mambaorg/micromamba:1.1.0 as base

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG PIP_NO_CACHE_DIR=0

RUN micromamba install -y \
    -c conda-forge \
    -c bioconda \
    "python>=3.9" \
    r=4.2 \
    cooler=0.8.11 \
    "matplotlib<3.6" \
    "numpy<1.24" \
    pandas=1.5.2 \
    rpy2 \
    r-wcorr=1.9.5 \
    scikit-image=0.19.3 \
    scipy=1.9.3 \
    xz \
&& micromamba clean --all -y

RUN touch /opt/conda/lib/R/etc/.Rprofile

ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/usr/bin/bash"]
WORKDIR /data

# We have to explicitly set these R_* env variables in order for the
# container to work correctly when running using Apptainer
ENV R_HOME=/opt/conda/lib/R
ENV R_LIBS=/opt/conda/lib/R/lib
ENV R_ENVIRON=/opt/conda/lib/R/etc/Renviron
ENV R_HISTFILE=/tmp/.Rhistory

ENV R_HOME_USER='$R_HOME'
ENV R_LIBS_USER='$R_LIBS'
ENV R_ENVIRON_USER='$R_ENVIRON'
ENV R_PROFILE_USER=/opt/conda/lib/R/etc/.Rprofile

# https://github.com/opencontainers/image-spec/blob/main/annotations.md#pre-defined-annotation-keys
LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/modle'
LABEL org.opencontainers.image.documentation='https://github.com/paulsengroup/modle'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/modle'
LABEL org.opencontainers.image.licenses='MIT'
