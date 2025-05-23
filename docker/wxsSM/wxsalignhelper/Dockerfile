##########################################
FROM rocker/r-ver:4.1.3

ARG SAM_VER=1.15.1
ENV JAVA_LIBRARY_PATH /usr/lib/jni
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/

RUN apt-get update && \
 	apt-get upgrade -y && \
	apt-get install -y --no-install-recommends \
  build-essential \
  git \
  wget \
  bedtools \
  libz-dev \
	libbz2-dev \
  rsync \
	python2 \
  libncurses5-dev \
  libncursesw5-dev \
  zlib1g \
  zlib1g-dev \
  liblzma-dev \
	openjdk-8-jdk \
  software-properties-common && \
  ln -s /usr/bin/python2.7 /usr/bin/python && \
  apt-get -y clean && \
	apt-get -y autoclean && \
	apt-get -y autoremove

RUN install2.r gfortran gsalib

RUN git clone --recursive https://github.com/bwa-mem2/bwa-mem2.git \
  && cd bwa-mem2 \
  && make \
  && cp bwa-mem2* /usr/local/bin

RUN wget -q -P /tmp https://github.com/samtools/samtools/releases/download/${SAM_VER}/samtools-${SAM_VER}.tar.bz2 \
  && tar xjvf /tmp/samtools-${SAM_VER}.tar.bz2 -C /tmp \
  && cd /tmp/samtools-${SAM_VER}/ \
  && ./configure \
  && make \
  && cp /tmp/samtools-${SAM_VER}/samtools /usr/local/bin

COPY ./scripts/* /scripts/

CMD ["/bin/bash"]
