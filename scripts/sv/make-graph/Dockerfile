FROM ubuntu:18.04
MAINTAINER vgteam
ARG vg_git_revision=8a4011928147a488e8d5cf3b593c156f78d2f266

# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true
ARG THREADS=4

# Install the base packages needed to let vg install packages.
# Make sure this runs after vg sources are imported so vg will always have an
# up to date package index to get its dependencies.
# We don't need to clean the package index since we don't ship this image and
# don't care about its size.
# We don't want to install too much stuff here, because we want to test vg's
# make get-deps to make sure it isn't missing something
RUN apt-get -qq -y update && \
    apt-get -qq -y upgrade && \
    apt-get -qq -y install \
    make \
    sudo \
    git

# fetch the desired git revision of vg
RUN git clone https://github.com/vgteam/vg.git /vg
WORKDIR /vg
RUN git fetch --tags origin && git checkout "$vg_git_revision" && git submodule update --init --recursive

# If we're trying to build from a non-recursively-cloned repo, go get the
# submodules.
RUN bash -c "[[ -e deps/sdsl-lite/CMakeLists.txt ]] || git submodule update --init --recursive"

# To increase portability of the docker image, set the target CPU architecture to
# Nehalem (2008) rather than auto-detecting the build machine's CPU.
# This has no AVX1, AVX2, or PCLMUL, but it does have SSE4.2.
# UCSC has a Nehalem machine that we want to support.
RUN sed -i s/march=native/march=nehalem/ deps/sdsl-lite/CMakeLists.txt
# Do the build. Trim down the resulting binary but make sure to include enough debug info for profiling.
RUN make get-deps && . ./source_me.sh && env && make include/vg_git_version.hpp && CXXFLAGS=" -march=nehalem " make -j $((THREADS < $(nproc) ? THREADS : $(nproc))) && make static && strip -d bin/vg

ENV PATH /vg/bin:$PATH


## python, snakemake and awscli
RUN apt-get update \
        && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        git screen wget curl gcc less nano \
        python3 \
        python3-pip \
        python3-setuptools \
        python3-dev \
        make \
        libncurses5-dev libncursesw5-dev \
        zlib1g-dev libbz2-dev liblzma-dev \
        && rm -rf /var/lib/apt/lists/*

RUN pip3 install --upgrade pip

RUN pip3 install --no-cache-dir requests awscli snakemake==5.8.2 biopython pyfaidx pyvcf pandas

RUN apt-get update \
        && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        bzip2 \
        && rm -rf /var/lib/apt/lists/*

## bcftools
RUN wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 && \
        tar -xjf bcftools-1.10.2.tar.bz2 && \
        cd bcftools-1.10.2 && \
        ./configure && make && make install && \
        cd .. && rm -rf bcftools-1.10.2 bcftools-1.10.2.tar.bz2

WORKDIR /home

COPY run.sh /home
