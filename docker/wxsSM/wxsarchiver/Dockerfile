##########################################
FROM rocker/r-ver:4.1.3

ARG BWA_VER=0.7.17
ARG SAM_VER=1.10
ARG BCF_VER=1.15.1
ENV JAVA_LIBRARY_PATH /usr/lib/jni
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/

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
	python2 \
  python3-matplotlib \
  texlive-latex-base \
  texlive-latex-recommended \
  texlive-latex-extra \
  texlive-fonts-recommended \
  bedtools \
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

RUN install2.r gfortran ggplot2 gsalib gplots reshape

RUN git clone --recursive https://github.com/bwa-mem2/bwa-mem2.git \
  && cd bwa-mem2 \
  && make \
  && cp bwa-mem2* /usr/local/bin

RUN wget --no-check-certificate -q -P /tmp https://github.com/lh3/bwa/releases/download/v${BWA_VER}/bwa-${BWA_VER}.tar.bz2 \
    && tar xjvf /tmp/bwa-${BWA_VER}.tar.bz2 -C /tmp \
    && cd /tmp/bwa-${BWA_VER}/ \
    && make \
    && cp /tmp/bwa-${BWA_VER}/bwa /usr/local/bin

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
