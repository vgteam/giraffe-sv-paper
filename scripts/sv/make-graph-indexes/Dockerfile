FROM ubuntu:18.04
MAINTAINER vgteam

# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

## python, snakemake and awscli
RUN apt-get update \
        && apt-get install -y --no-install-recommends \
        git screen wget curl gcc less nano \
        sudo \
        pkg-config \
        python3 \
        python3-pip \
        python3-setuptools \
        python3-dev \
        make \
        pigz tabix \
        bzip2 \
        libncurses5-dev libncursesw5-dev \
        zlib1g-dev libbz2-dev liblzma-dev \
        && rm -rf /var/lib/apt/lists/*

RUN pip3 install --upgrade pip

RUN pip3 install --no-cache-dir requests awscli snakemake==5.8.2 biopython pyfaidx pyvcf pandas boto3 numpy

## bcftools
RUN wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 && \
        tar -xjf bcftools-1.10.2.tar.bz2 && \
        cd bcftools-1.10.2 && \
        ./configure && make && make install && \
        cd .. && rm -rf bcftools-1.10.2 bcftools-1.10.2.tar.bz2

WORKDIR /bin

RUN wget --no-check-certificate https://github.com/vgteam/vg/releases/download/v1.28.0/vg && \
        chmod +x vg

WORKDIR /home
