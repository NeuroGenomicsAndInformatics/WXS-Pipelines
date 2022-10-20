##########################################
FROM ubuntu:20.04

ARG SAM_VER=1.10
ARG BCF_VER=1.15.1

RUN apt-get update && \
 	apt-get upgrade -y && \
	apt-get install -y --no-install-recommends \
  build-essential \
  git \
  wget \
	unzip \
	less \
  libz-dev \
	libbz2-dev \
  rsync \
  tabix \
	bedtools \
  libncurses5-dev \
  libncursesw5-dev \
  zlib1g \
  zlib1g-dev \
  liblzma-dev && \
  apt-get -y clean && \
	apt-get -y autoclean && \
	apt-get -y autoremove

RUN wget -q -P /tmp https://github.com/samtools/samtools/releases/download/${SAM_VER}/samtools-${SAM_VER}.tar.bz2 \
  && tar xjvf /tmp/samtools-${SAM_VER}.tar.bz2 -C /tmp \
  && cd /tmp/samtools-${SAM_VER}/ \
  && ./configure \
  && make \
  && cp /tmp/samtools-${SAM_VER}/samtools /usr/local/bin

RUN wget -q -P /tmp https://github.com/samtools/bcftools/releases/download/${BCF_VER}/bcftools-${BCF_VER}.tar.bz2 \
  && tar xjvf /tmp/bcftools-${BCF_VER}.tar.bz2 -C /tmp \
  && cd /tmp/bcftools-${BCF_VER}/ \
  && ./configure \
  && make \
  && cp /tmp/bcftools-${BCF_VER}/bcftools /usr/local/bin \
  && cp /tmp/bcftools-${BCF_VER}/misc/plot-vcfstats /usr/local/bin

COPY ./scripts/* /scripts/

CMD ["/bin/bash"]
