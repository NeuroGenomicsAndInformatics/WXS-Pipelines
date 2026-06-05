##########################################
FROM rocker/r-ver:4.5.3

ARG SAM_VER=1.21.1
ENV JAVA_LIBRARY_PATH=/usr/lib/jni
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64/

RUN apt-get update && \
 	apt-get upgrade -y && \
	apt-get install -y --no-install-recommends \
  build-essential \
  wget \
  libz-dev \
	libbz2-dev \
  rsync \
  tabix \
	python3 \
  libncurses5-dev \
  libncursesw5-dev \
  zlib1g \
  zlib1g-dev \
  liblzma-dev \
	openjdk-17-jre \
  software-properties-common && \
  apt-get -y clean && \
	apt-get -y autoclean && \
	apt-get -y autoremove

RUN install2.r gfortran gsalib

RUN wget -q -P /tmp https://github.com/samtools/samtools/releases/download/${SAM_VER}/samtools-${SAM_VER}.tar.bz2 \
  && tar xjvf /tmp/samtools-${SAM_VER}.tar.bz2 -C /tmp \
  && cd /tmp/samtools-${SAM_VER}/ \
  && ./configure \
  && make \
  && cp /tmp/samtools-${SAM_VER}/samtools /usr/local/bin

RUN ln -s /usr/bin/python3 /usr/bin/python

COPY ./scripts/* /scripts/

CMD ["/bin/bash"]
