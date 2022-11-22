##########################################
FROM ubuntu:22.04

ARG BWA_VER=0.7.17
ARG SAM_VER=1.10
ARG BCF_VER=1.15.1
ENV JAVA_LIBRARY_PATH /usr/lib/jni
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/

RUN apt-get update && \
 	apt-get upgrade -y && \
	apt-get install -y --no-install-recommends \
  bc \
  libz-dev \
	libbz2-dev \
  rsync \
  tabix \
  bedtools \
  samtools \
  zlib1g \
  zlib1g-dev \
  liblzma-dev && \
  apt-get -y clean && \
	apt-get -y autoclean && \
	apt-get -y autoremove

COPY ./scripts/* /scripts/

CMD ["/bin/bash"]
