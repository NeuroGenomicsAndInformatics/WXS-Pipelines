##########################################
FROM rocker/r-ver:latest

ENV JAVA_LIBRARY_PATH /usr/lib/jni
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/

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
  ln -s /usr/bin/python2.7 /usr/bin/python && \
  apt-get -y clean && \
	apt-get -y autoclean && \
	apt-get -y autoremove

RUN install2.r gfortran ggplot2 gsalib gplots

COPY ./scripts/* /scripts/

CMD ["/bin/bash"]
