##########################################
FROM rocker/r-ver:4.1.3

ARG SAM_VER=1.15.1

RUN apt-get update && \
 	apt-get upgrade -y && \
	apt-get install -y --no-install-recommends \
  build-essential \
  wget \
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

RUN wget -q -P /tmp https://github.com/samtools/samtools/releases/download/${SAM_VER}/samtools-${SAM_VER}.tar.bz2 \
  && tar xjvf /tmp/samtools-${SAM_VER}.tar.bz2 -C /tmp \
  && cd /tmp/samtools-${SAM_VER}/ \
  && ./configure \
  && make \
  && cp /tmp/samtools-${SAM_VER}/samtools /usr/local/bin

COPY ./scripts/* /scripts/

CMD ["/bin/bash"]
