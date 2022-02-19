##########################################
FROM rocker/r-base:latest

ARG SAM_VER
ARG BED_VER

RUN apt-get update && \
 	apt-get upgrade -y && \
	apt-get install -y --no-install-recommends \
  build-essential \
  wget \
	unzip \
	less \
	bzip2 \
  rsync \
	python2 \
	openjdk-8-jdk \
	tabix \
	bwa \
  samtools \
  bedtools \
  software-properties-common && \
	apt-get -y clean && \
	apt-get -y autoclean && \
	apt-get -y autoremove

ENV JAVA_LIBRARY_PATH /usr/lib/jni
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/

COPY ./scripts/* /scripts/

CMD ["/bin/bash"]
